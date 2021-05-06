import numpy as np
from scipy.optimize import root, minimize_scalar
from enum import Enum

ge = 9.80665  # standard gravity which is used to calculate exhaust velocity


def transition_mat(t, t0):
    """
    state vector transfer matrix
    :param t: final time
    :param t0: initial time
    :return: the corresponding transfer matrix(6 x 6)
    """
    return np.block([[np.cos(t - t0) * np.eye(3), np.sin(t - t0) * np.eye(3)],
                     [-np.sin(t - t0) * np.eye(3), np.cos(t - t0) * np.eye(3)]])


def gamma_mat(t, t0):
    """
    matrix gamma defined in paper, used for calculate state vector with thrust integral
    :param t: final time
    :param t0: initial time
    :return: the corresponding gamma matrix(6 x 6)
    """
    return np.block([[np.sin(t - t0) * np.eye(3), -np.cos(t - t0) * np.eye(3)],
                     [np.cos(t - t0) * np.eye(3), np.sin(t - t0) * np.eye(3)]])


def lambda_t(t, t0, pv0, pr0):
    """
    calculate the costate vectors at t(dimensionless)
    :param t: desired t
    :param t0: initial t
    :param pr0: initial Pr, 3x1
    :param pv0: initial Pv, 3x1
    :return: Pv,Pr at t
    """
    lambda_0 = np.vstack((pv0, -pr0))
    temp = transition_mat(t, t0) @ lambda_0
    return temp[0:3, :], -temp[3:6, :]


class Status(Enum):
    Evaluate_t2 = 1
    PDI = 2
    Powered_descent = 3


class UPG(object):

    def __init__(
            self, body,
            r_f, v_f,  # target
            r_t, v_t, mass,  # current states
            max_thrust, min_throttle, specific_impulse, k,
            t_0, t_1,  # time
            pv0, pr0, t_f  # initial guesses
    ):
        # constants
        self.r0 = body.equatorial_radius
        self.g0 = body.surface_gravity

        # multipliers to normalize variables
        self.time_multiplier = np.sqrt(self.r0 / self.g0)
        self.position_multiplier = self.r0
        self.velocity_multiplier = np.sqrt(self.r0 * self.g0)

        # targets
        self.r_f = r_f / self.position_multiplier
        self.v_f = v_f / self.velocity_multiplier

        # states
        r_t = r_t / self.position_multiplier
        v_t = v_t / self.velocity_multiplier
        self.x = np.vstack((r_t, v_t))

        # other needed variables
        self.m_0 = mass
        self.T_max = max_thrust
        self.T_min = max_thrust * min_throttle
        self.min_throttle = min_throttle
        self.isp = specific_impulse
        self.t_0 = t_0 / self.time_multiplier
        self.t_1 = t_1 / self.time_multiplier   # t_1 is a normalized time interval
        self.k = k  # coefficient in bolza problem

        # guess
        tf = t_f / self.time_multiplier + self.t_0
        self.z = np.hstack((pv0.T, pr0.T, tf)).T

        # store variables to avoid repeat calculation
        self.__x_f = np.zeros((6, 1))
        self.__m_f = 0.

        # variables to be used later
        self.__t_2 = 0.

        # fix t1 and t2 in main guidance loop

        self.t_1_fix = self.t_0 + self.t_1
        self.t_2_fix = self.t_0

        # to avoid access variables before solved
        self.t_2_solved: bool = False
        self.solved: bool = False
        self.status = Status.Evaluate_t2

    def mass_t(self, t, t0, isp, m0, thrust):
        """
        calculate the expected mass at dimensionless time t
        :param t0: initial time
        :param t: dimensionless time
        :param isp: engine's specific impulse
        :param m0: initial mass
        :param thrust: engine thrust
        :return: mass at t
        """
        return m0 - (t - t0) * self.time_multiplier * thrust / (isp * ge)

    def x_t(self, x0, t, t0, t1, t2, pv0, pr0, thrust, min_thrust, isp, m0):
        """
        calculate state vector at given time t,this method also update the final mass
        :param x0: initial state vector
        :param t: final time
        :param t0: initial time
        :param t1: first thrust switch time, which we turn thrust to T_min(thrust * min_t)
        :param t2: second thrust switch time, which we turn thrust back
        :param pv0: initial Pv
        :param pr0: initial Pr
        :param thrust: vessel's max thrust
        :param min_thrust: min thrust during descent; set to 0 for coast
        :param isp: specific impulse
        :param m0: initial mass
        :return: the corresponding satate vector
        """

        def i(t, t0, pv, pr, thrust, isp, m0):
            i_c = self.thrust_integral('c', t, t0, pv, pr, thrust, isp, m0)
            i_s = self.thrust_integral('s', t, t0, pv, pr, thrust, isp, m0)
            return np.vstack((i_c, i_s))

        def transfer(t, t0, x_0, pv0, pr0, thrust, isp, m0):
            i_t = i(t, t0, pv0, pr0, thrust, isp, m0)
            x_temp = transition_mat(t, t0) @ x_0 + gamma_mat(t, t0) @ i_t
            pv_t, pr_t = lambda_t(t, t0, pv0, pr0)
            m_t = self.mass_t(t, t0, isp, m0, thrust)
            return x_temp, pv_t, pr_t, m_t

        if t1 > t2:
            print("warning: t1>t2, t1=" + str(t1) + " t2=" + str(t2))
            t2 = t1
        if t0 < t1:
            # at t1
            x_t, pv_t, pr_t, m_t = transfer(t1, t0, x0, pv0, pr0, thrust, isp, m0)
            # at t2
            x_t, pv_t, pr_t, m_t = transfer(t2, t1, x_t, pv_t, pr_t, min_thrust, isp, m_t)
            # at t
            x_t, _, _, self.__m_f = transfer(t, t2, x_t, pv_t, pr_t, thrust, isp, m_t)
        else: # t0 >= t1:
            if t0 < t2:
                # at t2
                x_t, pv_t, pr_t, m_t = transfer(t2, t0, x0, pv0, pr0, min_thrust, isp, m0)
                # at t
                x_t, _, _, self.__m_f = transfer(t, t2, x_t, pv_t, pr_t, thrust, isp, m_t)
            else:  # t0 >= t2:
                x_t, _, _, self.__m_f = transfer(t, t0, x0, pv0, pr0, thrust, isp, m0)
        return x_t

    def thrust_integral(self, type_name, t, t0, pv0, pr0, thrust, isp, m0):
        """
        calculate thrust integral with milne's rule
        :param type_name: determins this is the cos integral or sin integral, should be 'c' or 's'
        :param t: end time
        :param t0: initial time
        :param pv0: initial pv
        :param pr0: initial pr
        :param thrust: engine thrust
        :param isp: engine isp
        :param m0: initial mass
        :return: the estimated thrust integral
        """

        def i(t):
            pv, _ = lambda_t(t, t0, pv0, pr0)
            pv = pv / np.linalg.norm(pv)
            temp = pv * thrust / (self.mass_t(t, t0, isp, m0, thrust) * self.g0)
            if type_name == 'c':
                return temp * np.cos(t - t0)
            elif type_name == 's':
                return temp * np.sin(t - t0)
            else:
                raise Exception()

        delta = (t - t0) / 4.
        return ((t - t0) / 90.) * (7 * i(t0) + 32 * i(t0 + delta) + 12 * i(t0 + 2 * delta) + 32 * i(t0 + 3 * delta) +
                                   7 * i(t0 + 4 * delta))

    def target_function_sl(self, z, r_f, v_f, x_0, t_0, t_1, t_2, thrust, min_thrust, isp, m0):
        """
        target function for soft landing
        :param z: solution to be optimized
        :param r_f: target position
        :param v_f: target velocity
        :param x_0: initial state vector
        :param t_0: initial time
        :param t_1: first thrust switch time
        :param t_2: second thrust switch time
        :param thrust: max thrust
        :param min_thrust: min thrust
        :param isp: specific impulse
        :param m0: initial mass
        :return: the 7 equations for solving z
        """
        return self.target_function_bl(z, 0, r_f, v_f, x_0, t_0, t_1, t_2, thrust, min_thrust, isp, m0)

    def target_function_bl(self, z, k, r_f, v_f, x_0, t_0, t_1, t_2, thrust, min_thrust, isp, m0):
        """
        target function for Bolza landing
        :param k: coefficient for performance
        :param z: solution to be optimized
        :param r_f: target position
        :param v_f: target velocity
        :param x_0: initial state vector
        :param t_0: initial time
        :param t_1: first thrust switch time
        :param t_2: second thrust switch time
        :param thrust: max thrust
        :param min_thrust: min thrust
        :param isp: specific impulse
        :param m0: initial mass
        :return: the 7 equations for solving z
        """
        pv = z[0:3].reshape(-1, 1)
        pr = z[3:6].reshape(-1, 1)
        tf = z[6]
        # final position constraint
        self.__x_f = self.x_t(x_0, tf, t_0, t_1, t_2, pv, pr, thrust, min_thrust, isp, m0)
        r_tf = self.__x_f[0:3, :]
        s1 = r_tf.T @ r_tf - r_f.T @ r_f
        # final velocity constraint
        v_tf = self.__x_f[3:6, :]
        s2 = v_tf - v_f  # s2 includes 3 equations
        # transversality condition
        pvf, prf = lambda_t(tf, t_0, pv, pr)

        s3 = prf.T @ v_tf - pvf.T @ r_tf + \
             np.linalg.norm(pvf) * thrust / (self.__m_f * self.g0) - 1
        s4 = r_tf[2, 0] * prf[0, 0] - r_tf[0, 0] * prf[2, 0] + k * (r_tf[0, 0] * r_f[2, 0] - r_tf[2, 0] * r_f[0, 0])
        s5 = r_tf[2, 0] * prf[1, 0] - r_tf[1, 0] * prf[2, 0] + k * (r_tf[1, 0] * r_f[2, 0] - r_tf[2, 0] * r_f[1, 0])
        return [s1[0, 0], s2[0, 0], s2[1, 0], s2[2, 0], s3[0, 0], s4, s5]

    def target_function_pl(self, z, r_f, v_f, x_0, t_0, t_1, t_2, thrust, min_thrust, isp, m0):
        # pinpoint landing is not robust. don't use it
        pv = z[0:3].reshape(-1, 1)
        pr = z[3:6].reshape(-1, 1)
        tf = z[6]
        # final position constraint
        self.__x_f = self.x_t(x_0, tf, t_0, t_1, t_2, pv, pr, thrust, min_thrust, isp, m0)
        r_tf = self.__x_f[0:3, :]
        s1 = r_tf - r_f  # s1 includes 3 equations
        # final velocity constraint
        v_tf = self.__x_f[3:6, :]
        s2 = v_tf - v_f  # s2 includes 3 equations
        # transversality condition
        pvf, prf = lambda_t(tf, t_0, pv, pr)
        s3 = prf.T @ v_tf - pvf.T @ r_tf + \
             np.linalg.norm(pvf) * thrust / (self.mass_t(tf, t_0, isp, m0, thrust) * self.g0) - 1
        return [s1[0, 0], s1[1, 0], s1[2, 0], s2[0, 0], s2[1, 0], s2[2, 0], s3[0, 0]]

    def target_function_t2_sl(self, t2, z0, r_f, v_f, x_0, t_0, t_1, thrust, min_thrust, isp, m0):
        # m_f is calculated during root finding
        root(
            self.target_function_sl, z0, args=(r_f, v_f, x_0, t_0, t_1, t2, thrust, min_thrust, isp, m0),
            jac=False)
        return m0 - self.__m_f

    def solve(self, mode: str):
        if not self.t_2_solved:
            raise Exception("t2 is not solved before solving guidance!")
        if mode == 'bolza':
            self.z = root(
                fun=self.target_function_bl, x0=self.z,
                args=(
                    self.k, self.r_f, self.v_f, self.x, self.t_0, self.t_1_fix, self.t_2_fix,
                    self.T_max, self.T_min, self.isp, self.m_0),
                jac=False).x.reshape(-1, 1)
        elif mode == 'soft':
            self.z = root(
                fun=self.target_function_sl, x0=self.z,
                args=(
                    self.r_f, self.v_f, self.x, self.t_0, self.t_1_fix, self.t_2_fix,
                    self.T_max, self.T_min, self.isp, self.m_0),
                jac=False).x.reshape(-1, 1)
        self.solved = True

    def solve_t_2(self):
        t_2 = minimize_scalar(  # bounds are set to between t1 and tf
            self.target_function_t2_sl, bounds=(self.t_1_fix, self.z[6, 0]),
            args=(self.z, self.r_f, self.v_f, self.x, self.t_0, self.t_1_fix, self.T_max, self.T_min,
                  self.isp, self.m_0), method='Bounded'
        ).x
        self.__t_2 = t_2 - self.t_0  # convert to time interval
        self.t_2_solved = True
        self.status = Status.PDI

    def update(self, r_t, v_t, mass, max_thrust, t_0):
        r_t = r_t / self.position_multiplier
        v_t = v_t / self.velocity_multiplier
        self.x = np.vstack((r_t, v_t))

        self.m_0 = mass
        if self.T_max != max_thrust:
            self.T_max = max_thrust
            self.T_min = max_thrust * self.min_throttle
        self.t_0 = t_0 / self.time_multiplier

    def update_time(self):
        self.t_1_fix = self.t_0 + self.t_1
        self.t_2_fix = self.t_0 + self.__t_2

    def set_target(self, r_f, v_f):
        self.r_f = r_f / self.position_multiplier
        self.v_f = v_f / self.velocity_multiplier

    @property
    def r_tf(self):
        if not self.solved:
            raise Exception("accessing variables before solved them")
        return self.__x_f[0:3, :] * self.position_multiplier

    @property
    def v_tf(self):
        if not self.solved:
            raise Exception("accessing variables before solved them")
        return self.__x_f[3:6, :] * self.velocity_multiplier

    @property
    def t_f(self):
        if not self.solved:
            raise Exception("accessing variables before solved them")
        return (self.z[6, 0] - self.t_0) * self.time_multiplier

    @property
    def t_2(self):
        if not self.t_2_solved:
            raise Exception("accessing variables before solved them")
        return self.__t_2 * self.time_multiplier

    @property
    def t_1_go(self):
        if not self.solved:
            raise Exception("accessing variables before solved them")
        return (self.t_1_fix - self.t_0) * self.time_multiplier

    @property
    def t_2_go(self):
        if not self.solved:
            raise Exception("accessing variables before solved them")
        return (self.t_2_fix - self.t_0) * self.time_multiplier

    @property
    def v_go(self):
        v_exh = self.isp * ge
        return v_exh * np.log(self.m_0 / self.__m_f)

    @property
    def thrust_direction(self):
        if not self.solved:
            raise Exception("accessing variables before solved them")
        pv = self.z[0:3, :]
        return pv / np.linalg.norm(pv)

    @property
    def throttle(self):
        if not self.solved:
            raise Exception("accessing variables before solved them")
        if self.t_1_fix < self.t_0 < self.t_2_fix:
            return 0.0001 * self.min_throttle
        else:
            return 1.
