
from enum import Enum
import numpy as np
from scipy.optimize import root, minimize_scalar
from utils import swap_yz

ge = 9.80665  # standard gravity which is used to calculate exhaust velocity


def transition_mat(dt: float) -> np.ndarray:
    """State transition matrix

    Args:
        dt (float): The time interval

    Returns:
        np.ndarray: The state transition matrix(3x3)
    """
    return np.block(
        [[np.cos(dt) * np.eye(3), np.sin(dt) * np.eye(3)],
         [-np.sin(dt) * np.eye(3), np.cos(dt) * np.eye(3)]])


def gamma_mat(t: float) -> np.ndarray:
    """Gamma matrix is defined as\n
    [sin(t) -cos(t)
     cos(t)  sin(t)] \n
     as R^6x6. Note that the input is absolute time, not delta time.

    Args:
        t (float): The final Time

    Returns:
        np.ndarray: The corresponding gamma matrix(6 x 6)
    """
    return np.block([[np.sin(t) * np.eye(3), -np.cos(t) * np.eye(3)],
                     [np.cos(t) * np.eye(3), np.sin(t) * np.eye(3)]])


def lambda_t(t: float, t0: float, p_r_0: np.ndarray, p_v_0: np.ndarray):
    """Calculate the costate vectors at the given time.

    Args:
        t (float): The final time
        t0 (float): The initial time
        p_r_0 : The initial costate of r(3x1)
        p_v_0 : The initial costate of v(3x1)

    Returns:
        Vector p_r and p_v at the given time, respectively.
    """
    lambda_0 = np.vstack((p_r_0, p_v_0))
    temp = transition_mat(t - t0) @ lambda_0
    p_r_t, p_v_t = np.vsplit(temp, 2)
    return p_r_t, p_v_t


def K(t, p_v):
    Lambda = np.split(transition_mat(t), 2)[-1]
    norm = np.linalg.norm(p_v)
    _1_pv = p_v / norm
    temp = 1. / norm * (np.eye(3) - _1_pv @ _1_pv.T)
    return temp @ Lambda


class Status(Enum):
    Evaluate_t2 = 1
    PDI = 2
    Powered_descent = 3
    Terminal = 4
    Finished = 5


class UPG(object):

    def __init__(self, r_0: float, g_0: float, r_target, v_target,  # target
                 r_t, v_t, mass: float,  # current states
                 max_thrust: float, min_throttle: float, specific_impulse: float, k: float,
                 t_0: float, t_1: float,  # time
                 p_r_0, p_v_0, t_f: float, t_2_bound: list  # initial guesses
                 ):
        """The UPG class encapsules the calculation of the UPG algorithm.

        Args:
            r_0 (float): The equatorial radius of the celestial body.
            g_0 (float): The surface gravity of the body.
            r_target (3x1): Vector pointing at target position.
            v_target (3x1): The target velocity vector.
            r_t (3x1): Current position.
            v_t (3x1): Current velocity.
            mass (float): Spacecraft's mass.
            min_throttle (float): min throttle of the engine.
            specific_impulse (float): engine's vaccum specific_impulse.
            k (float): Coefficient for bolza landing.
            t_0 (float): The current time.
            t_1 (float): Duration of first burn.
            p_r_0 (3x1): Costate vector of r(guess)
            p_v_0 (3x1): Costate vector of v(guess)
            t_f (float): Total time duration(guess)
            t_2_bound (list): Upper and lower bound of t_2.
        """
        # constants
        self.g_0 = g_0

        # scalers to normalize variables for better numerical condition
        self.t_scaler = np.sqrt(r_0 / self.g_0)
        self.r_scaler = r_0
        self.v_scaler = np.sqrt(r_0 * self.g_0)

        # targets
        self.r_target = r_target / self.r_scaler
        self.v_target = v_target / self.v_scaler

        # states
        r_t = r_t / self.r_scaler
        v_t = v_t / self.v_scaler
        self.x = np.vstack((r_t, v_t))

        # other needed variables
        self.m = mass
        self.T_max = max_thrust
        self.T_min = max_thrust * min_throttle
        self.min_throttle = min_throttle
        self.isp = specific_impulse
        self.t_0 = 0.
        # t_1 is a normalized time interval
        self.t_1 = t_1 / self.t_scaler
        self.k = k  # coefficient in bolza problem
        self.t_2_bound =\
            [value / self.t_scaler for value in t_2_bound]  # boundary for t_2

        # guess
        t_f = t_f / self.t_scaler
        self.z = np.hstack((p_r_0.T, p_v_0.T, t_f)).T

        # store variables to avoid repeat calculation
        self.x_final = np.zeros((6, 1))
        self.m_final = 0.
        self.pv_final = np.zeros((3, 1))
        self.pr_final = np.zeros((3, 1))

        # variables to be used later
        self.t_2 = 0.
        self.t_solved = 0.

        # save the start time
        self.t_start = t_0 / self.t_scaler

        self.status = Status.Evaluate_t2
        self.convergence = False
        self.solver_status = 0
        self.last_err_msg = ""

    def mass_t(self, t: float, t0: float, m_0: float, thrust: float) -> float:
        """Calculate the mass at the given time t and initial mass m_0.

        Args:
            t (float): The final time
            t0 (float): The initial time
            m_0 (float): The initial mass
            thrust (float): Engine's thrust during this time duration

        Returns:
            float: The spacecraft's mass at time t.
        """
        return m_0 - (t - t0) * self.t_scaler * thrust / (self.isp * ge)

    def i(self, t: float, t_0: float, p_r_0, p_v_0, thrust: float, m_0: float):
        """The integrand of thrust integral. The variable correspond to time t 
            are first computed by t_0, then the integrand is evaluated.

        Args:
            t (float): Time at which the integrand is computed
            t_0 (float): The initial time.
            p_r_0 (ndarray): Initial costate vectors.
            p_v_0 (ndarray): Initial costate vectors.
            thrust (float): Engine's thrust during integration.
            m_0 (float): Initial mass.

        Returns:
            ndarray: The stacked integrand[i_c, i_s].T
        """
        _, p_v_t = lambda_t(t, t_0, p_r_0, p_v_0)
        p_v_t /= np.linalg.norm(p_v_t)
        m_t = self.mass_t(t, t_0, m_0, thrust)
        temp = p_v_t * thrust / (m_t * self.g_0)
        return np.vstack((temp * np.cos(t), temp * np.sin(t)))

    def j_i(self, t: float, t_0: float, p_r_0, p_v_0, thrust: float, m_0: float):
        _, p_v_t = lambda_t(t, t_0, p_r_0, p_v_0)
        K_t = K(t, p_v_t)
        m_t = self.mass_t(t, t_0, m_0, thrust)
        temp = thrust / (m_t * self.g_0)
        out = temp * np.block([[np.cos(t) * np.eye(3)],
                              [np.sin(t) * np.eye(3)]]) @ K_t
        return out

    def thrust_integral(self, t: float, t_0: float, p_r_0, p_v_0, thrust: float, m_0: float):
        """Calculate thrust integral with milne's rule

        Args:
            t (float): Upper bound of integral
            t_0 (float): Lower bound of integral
            p_r_0 (ndarray): Costate vectors at lower bound, 3x1
            p_v_0 (ndarray): Costate vectors at lower bound, 3x1
            thrust (float): Engine's thrust during integration
            m_0 (float): Initial mass at lower bound

        Returns:
            ndarray: a 6x1 vector representing the thrust integral, 
            the cos term on the top and the sin term on the bottom.
        """
        coefs = [7., 32., 12., 32., 7.]
        delta = (t - t_0) / 4.
        inte = np.asarray([coefs[i] * self.i(t_0 + i * delta, t_0, p_r_0, p_v_0, thrust, m_0)
                           for i in range(len(coefs))])
        jac = np.asarray([coefs[i] * self.j_i(t_0 + i * delta, t_0, p_r_0, p_v_0, thrust, m_0)
                          for i in range(len(coefs))])
        return ((t - t_0) / 90.) * np.sum(inte, axis=0), ((t - t_0) / 90.) * np.sum(jac, axis=0)

    def state_at_t(self, t: float, t0: float, x_0, p_r_0, p_v_0, thrust: float, m_0: float):
        """Calculate the state `x`, costate vectors `p_r` & `p_v` and the mass `m` 
        at the given time `t`.

        Args:
            t (float): The final time
            t0 (float): The initial time
            x_0 (_type_): State vector at the initial time
            p_r_0 (_type_): Initial costate vectors
            p_v_0 (_type_): Initial costate vectors
            thrust (float): Engine's thrust during this time
            m_0 (float): Initial mass of the spacecraft at initial time.

        Returns:
            The state `x`, costate vectors `p_r` & `p_v` and the mass `m` 
        at the given time `t`, respectively.
        """
        I_t, j_i = self.thrust_integral(t, t0, p_r_0, p_v_0, thrust, m_0)
        x_t = transition_mat(t - t0) @ x_0 + gamma_mat(t) @ I_t
        j_i = gamma_mat(t) @ j_i
        p_r_t, p_v_t = lambda_t(t, t0, p_r_0, p_v_0)
        m_t = self.mass_t(t, t0, m_0, thrust)
        return x_t, p_r_t, p_v_t, m_t, j_i

    def update_final_state(self, t_f: float, p_r_0, p_v_0):
        """Calculate the final state of the guidance. The final costate vectors
        and the final mass are updated togather.

        Args:
            t_f (float): The final time as interval from current time
        """
        # t_elapsed = self.t_0 - self.t_start
        t_1 = self.t_1 - self.t_0
        t_2 = self.t_1 + self.t_2 - self.t_0
        t_f = t_f - self.t_0
        if 0. < t_1:  # first burn arc
            # state at t1
            x_t, p_r_t, p_v_t, m_t, j_i = self.state_at_t(
                t_1, 0., self.x, p_r_0, p_v_0, self.T_max, self.m)
            dx1_dlam = j_i
            
            # state at t2
            x_t, p_r_t, p_v_t, m_t, j_i = self.state_at_t(
                t_2, t_1, x_t, p_r_t, p_v_t, self.T_min, m_t)
            dx2_dlam = transition_mat(t_2 - t_1) @ dx1_dlam + j_i
            
            # state at t_final
            self.x_final, self.pr_final, self.pv_final, self.m_final, j_i = self.state_at_t(
                t_f, t_2, x_t, p_r_t, p_v_t, self.T_max, m_t)
            
            dxf_dlam = transition_mat(t_f - t_2) @ dx2_dlam + j_i
            
            dxf_dtf = -gamma_mat(t_f - t_2) @ x_t +\
                transition_mat(t_f) @ self.thrust_integral(t_f, t_2, p_r_t, p_v_t, self.T_max, m_t)[0] +\
                gamma_mat(t_f) @ self.i(t_f, t_f, self.pr_final, self.pv_final, self.T_max, self.m_final)

        elif 0. < t_2:  # second burn arc, namely throttle back phase
            # state at t2
            x_t, p_r_t, p_v_t, m_t, j_i = self.state_at_t(
                t_2, 0., self.x, p_r_0, p_v_0, self.T_min, self.m)
            dx2_dlam = j_i
            
            # at t
            self.x_final, self.pr_final, self.pv_final, self.m_final, j_i = self.state_at_t(
                t_f, t_2, x_t, p_r_t, p_v_t, self.T_max, m_t)
            
            dxf_dlam = transition_mat(t_f - t_2) @ dx2_dlam + j_i
            
            dxf_dtf = -gamma_mat(t_f-t_2) @ x_t +\
                transition_mat(t_f) @ self.thrust_integral(t_f, t_2, p_r_t, p_v_t, self.T_max, m_t)[0] +\
                gamma_mat(t_f) @ self.i(t_f, t_f, self.pr_final, self.pv_final, self.T_max, self.m_final)

        else:  # final burn
            self.x_final, self.pr_final, self.pv_final, self.m_final, j_i = self.state_at_t(
                t_f, 0., self.x, p_r_0, p_v_0, self.T_max, self.m)
            
            dxf_dlam = j_i
            
            dxf_dtf = -gamma_mat(t_f - 0.) @ self.x +\
                transition_mat(t_f) @ self.thrust_integral(t_f, 0, p_r_0, p_v_0, self.T_max, self.m)[0] +\
                gamma_mat(t_f) @ self.i(t_f, t_f, self.pr_final, self.pv_final, self.T_max, self.m_final)

        return dxf_dlam, dxf_dtf

    def target_function_bl(self, z, k, use_jac=False):
        p_r_0, p_v_0, t_f = np.vsplit(z.reshape(7, 1), [3, 6])

        dxf_dlam, dxf_dtf = self.update_final_state(t_f, p_r_0, p_v_0)
        r_final, v_final = np.vsplit(self.x_final, 2)
        pvf_norm = np.linalg.norm(self.pv_final)

        # final position constraint
        s_1 = (r_final.T @ r_final - self.r_target.T @ self.r_target).flatten()

        # final velocity constraint
        s_2 = (v_final - self.v_target).flatten()
        # transversality condition
        s_3 = (self.pr_final.T @ v_final - self.pv_final.T @ r_final +
               pvf_norm * self.T_max / (self.m_final * self.g_0) - 1.).flatten()
        # reduced transversality conditions
        a_1 = a_2 = None
        if r_final[2, 0] != 0.:  # r_3 not 0
            a_1 = np.array(
                [[0,  0,  1.],
                 [0,  0,  0],
                 [-1., 0,  0]]
            )
            a_2 = np.array(
                [[0,  0,  0],
                 [0,  0,  1.],
                 [0, -1., 0]]
            )
        elif r_final[1, 0] != 0.:  # r_2 not 0
            a_1 = np.array(
                [[0,  1., 0],
                 [-1., 0,  0],
                 [0,  0,  0]]
            )
            a_2 = np.array(
                [[0,  0,  0],
                 [0,  0, -1.],
                 [0,  1., 0]]
            )
        else:  # r_1
            a_1 = np.array(
                [[0, -1., 0],
                 [1., 0,  0],
                 [0,  0,  0]]
            )
            a_2 = np.array(
                [[0,  0, -1.],
                 [0,  0,  0],
                 [1., 0,  0]]
            )
        s_4 = ((a_1 @ r_final).T @ (self.pr_final +
                                   2. * k * (r_final - self.r_target))).flatten()
        s_5 = ((a_2 @ r_final).T @ (self.pr_final +
                                   2. * k * (r_final - self.r_target))).flatten()
        # normalize s4 and s5 to conpensate for k
        scaler = np.sqrt(1. + k ** 2)
        s_4 /= scaler
        s_5 /= scaler
        s = np.concatenate([s_1, s_2, s_3, s_4, s_5])
        s = list(s)
        if use_jac is True:
            # calculate jacobian
            dr_dlam, dv_dlam = np.split(dxf_dlam, 2, axis=0)
            dr_dtf, dv_dtf = np.split(dxf_dtf, 2, axis=0)
            dpr_dlam, dpv_dlam = np.split(transition_mat(t_f), 2, axis=0)
            dpr_dtf, dpv_dtf = np.split(
                -gamma_mat(t_f) @ np.concatenate([p_r_0, p_v_0], axis=0), 2, axis=0)
            

            j_1 = 2. * np.block([r_final.T @ dr_dlam, r_final.T @ dr_dtf])
            j_2 = np.eye(3) @ np.block([dv_dlam, dv_dtf])
            j_3 = np.block([
                v_final.T @ dpr_dlam + self.pr_final.T @ dv_dlam -
                r_final.T @ dpv_dlam - self.pv_final.T @ dr_dlam +
                self.T_max / self.m_final / self.g_0 /
                pvf_norm * (self.pv_final.T @ dpv_dlam),
                
                v_final.T @ dpr_dtf + self.pr_final.T @ dv_dtf -
                r_final.T @ dpv_dtf - self.pv_final.T @ dr_dtf +
                self.T_max / self.m_final / self.g_0 / pvf_norm * (self.pv_final.T @ dpv_dtf) +
                self.T_max ** 2 * pvf_norm * self.t_scaler /
                self.m_final ** 2 / self.g_0 / (self.isp * ge)
            ])
            temp1 = (self.pr_final + 2. * k * (r_final - self.r_target)).T
            j_4 = np.block([
                temp1 @ a_1 @ dr_dlam +
                r_final.T @ a_1.T @ (dpr_dlam + 2. * k * dr_dlam),
                temp1 @ a_1 @ dr_dtf +
                r_final.T @ a_1.T @ (dpr_dtf + 2. * k * dr_dtf)
            ])
            j_5 = np.block([
                temp1 @ a_2 @ dr_dlam +
                r_final.T @ a_2.T @ (dpr_dlam + 2. * k * dr_dlam),
                temp1 @ a_2 @ dr_dtf +
                r_final.T @ a_2.T @ (dpr_dtf + 2. * k * dr_dtf)
            ])
            j_4 /= scaler
            j_5 /= scaler

            jac = np.concatenate([j_1, j_2, j_3, j_4, j_5])
            return (s, jac)
        return s

    def target_function_pl(self, z, use_jac=False):
        p_r_0, p_v_0, t_f = np.vsplit(z.reshape(7, 1), [3, 6])

        dxf_dlam, dxf_dtf = self.update_final_state(t_f, p_r_0, p_v_0)
        r_final, v_final = np.vsplit(self.x_final, 2)
        pvf_norm = np.linalg.norm(self.pv_final)

        # final position constraint
        scaler = 1e-3
        s_1 = scaler * (r_final - self.r_target).flatten() # scale down for better convergence

        # final velocity constraint
        s_2 = (v_final - self.v_target).flatten()
        # transversality condition
        s_3 = (self.pr_final.T @ v_final - self.pv_final.T @ r_final +
               pvf_norm * self.T_max / (self.m_final * self.g_0) - 1.).flatten()
        s = np.concatenate([s_1, s_2, s_3])
        s = list(s)
        
        if use_jac is True:
            # calculate jacobian
            dr_dlam, dv_dlam = np.split(dxf_dlam, 2, axis=0)
            dr_dtf, dv_dtf = np.split(dxf_dtf, 2, axis=0)
            dpr_dlam, dpv_dlam = np.split(transition_mat(t_f), 2, axis=0)
            dpr_dtf, dpv_dtf = np.split(
                -gamma_mat(t_f) @ np.concatenate([p_r_0, p_v_0], axis=0), 2, axis=0)
            
            j_1 = scaler * np.eye(3) @ np.block([dr_dlam, dr_dtf])
            j_2 = np.eye(3) @ np.block([dv_dlam, dv_dtf])
            j_3 = np.block([
                v_final.T @ dpr_dlam + self.pr_final.T @ dv_dlam -
                r_final.T @ dpv_dlam - self.pv_final.T @ dr_dlam +
                self.T_max / self.m_final / self.g_0 /
                pvf_norm * (self.pv_final.T @ dpv_dlam),
                v_final.T @ dpr_dtf + self.pr_final.T @ dv_dtf -
                r_final.T @ dpv_dtf - self.pv_final.T @ dr_dtf +
                self.T_max / self.m_final / self.g_0 / pvf_norm * (self.pv_final.T @ dpv_dtf) +
                self.T_max ** 2 * pvf_norm * self.t_scaler /
                self.m_final ** 2 / self.g_0 / (self.isp * ge)
            ])

            jac = np.concatenate([j_1, j_2, j_3])
            return (s, jac)
        return s
        
    def target_function_t_2(self, t_2):
        # m_f is calculated during root finding
        self.t_2 = t_2
        # while True:
        sol = root(self.target_function_bl, self.z, args=(0.,True), jac=True)
            # if sol.success:
            #     break
        self.z = sol.x.reshape(-1, 1)
        return self.m - self.m_final

    def solve(self, mode: str, jac=True):
        sol = None
        z = self.new_z()
        k = self.k if mode == 'bolza' else 0.
        if mode == 'pinpoint':
            sol = root(fun=self.target_function_pl, x0=z, args=(jac,), jac=jac)
        else:
            sol = root(fun=self.target_function_bl, x0=z, args=(k, jac), jac=jac)        

        if sol.success:
            z = sol.x
            self.z = z.reshape(7, 1)
            self.t_solved = self.t_0

        self.norm = np.linalg.norm(sol.fun)
        self.fun = sol.fun
        self.solver_status = sol.status
        self.convergence = sol.success
        self.last_err_msg = sol.message

    def solve_t_2(self):
        """Calculates the optimal t_2.
        """
        print("Terminate program if it stucks at solving t_2.")
        sol = minimize_scalar(
            self.target_function_t_2, bounds=self.t_2_bound, method='Bounded')
        # self.t_2 = sol.x
        # t_2 is updated during solving, no need to assign.
        if sol.success:
            self.status = Status.PDI
        else:
            print("Solving t_2 failed, retry...")
            self.solve_t_2()

    def update_state(self, r_t, v_t, mass, max_thrust, t_0):
        r_t = r_t / self.r_scaler
        v_t = v_t / self.v_scaler
        self.x = np.vstack((r_t, v_t))

        self.m = mass
        if self.T_max != max_thrust:
            self.T_max = max_thrust
            self.T_min = max_thrust * self.min_throttle

        self.t_0 = t_0 / self.t_scaler - self.t_start

    def update_start_time(self, t_0):
        self.t_start = t_0 / self.t_scaler

    def new_z(self):
        p_r_0, p_v_0, t_f = np.vsplit(self.z, [3, 6])
        p_r, p_v = lambda_t(self.t_0, self.t_solved, p_r_0, p_v_0)
        return np.concatenate([p_r, p_v, t_f])

    def set_target(self, r_f, v_f):
        self.r_target = r_f / self.r_scaler
        self.v_target = v_f / self.v_scaler

    @property
    def r_final(self):
        return self.x_final[0:3, :] * self.r_scaler

    @property
    def v_final(self):
        return self.x_final[3:6, :] * self.v_scaler

    @property
    def t_f(self):
        return (self.z[6, 0] - self.t_0).item() * self.t_scaler

    @property
    def get_t_2(self):
        return self.t_2.item() * self.t_scaler

    @property
    def t_1_go(self):
        return (self.t_1 - self.t_0).item() * self.t_scaler

    @property
    def t_2_go(self):
        return (self.t_2 + self.t_1 - self.t_0).item() * self.t_scaler

    @property
    def v_go(self):
        v_exh = self.isp * ge
        return v_exh * np.log(self.m / self.m_final).item()

    @property
    def last_convergence(self):
        return (self.t_0 - self.t_solved) * self.t_scaler

    @property
    def thrust_direction(self):
        z = self.new_z()
        pr, pv, tf = np.vsplit(z, [3, 6])
        pv = pv / np.linalg.norm(pv)
        # conver back to left handed
        pv = swap_yz(pv)
        return tuple(pv.flatten())

    @property
    def throttle(self):
        if self.t_1 < self.t_0 < self.t_1 + self.t_2:
            return 0.0001 if self.k != 0. else 0.
        else:
            return 1.
