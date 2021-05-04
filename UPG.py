import numpy as np
from scipy.optimize import root

ge = 9.80665  # standard gravity which is used to calculate exhaust velocity




class UPG(object):

    def __init__(self, conn, ):
        self.conn = conn
        self.space_center = conn.space_center
        self.vessel = self.space_center.active_vessel
        self.ap = self.vessel.auto_pilot
        self.body = self.vessel.orbit.body
        self.body_frame = self.body.reference_frame
        self.flight = self.vessel.flight(self.body_frame)
        self.r0 = self.body.equatorial_radius
        self.g0 = self.body.surface_gravity
        # multipliers to normalize variables
        self.time_multiplier = np.sqrt(self.r0 / self.g0)
        self.position_multiplier = self.r0
        self.velocity_multiplier = np.sqrt(self.r0 * self.g0)
        # stored in class to avoid duplicate calculation
        self.x_f = np.zeros((6, 1))
        self.m_f = 0

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

    def transfer_mat(self, t, t0):
        """
        state vector transfer matrix
        :param t: final time
        :param t0: initial time
        :return: the corresponding transfer matrix(6 x 6)
        """
        return np.block([[np.cos(t - t0) * np.eye(3), np.sin(t - t0) * np.eye(3)],
                         [-np.sin(t - t0) * np.eye(3), np.cos(t - t0) * np.eye(3)]])

    def gamma_mat(self, t, t0):
        """
        matrix gamma defined in paper, used for calculate state vector with thrust integral
        :param t: final time
        :param t0: initial time
        :return: the corresponding gamma matrix(6 x 6)
        """
        return np.block([[np.sin(t - t0) * np.eye(3), -np.cos(t - t0) * np.eye(3)],
                         [np.cos(t - t0) * np.eye(3), np.sin(t - t0) * np.eye(3)]])

    def lambda_t(self, t, t0, pv0, pr0):
        """
        calculate the costate vectors at t(dimensionless)
        :param t: desired t
        :param t0: initial t
        :param pr0: initial Pr, 3x1
        :param pv0: initial Pv, 3x1
        :return: Pv,Pr at t
        """
        lambda_0 = np.vstack((pv0, -pr0))
        temp = self.transfer_mat(t, t0) @ lambda_0
        return temp[0:3, :], temp[3:6, :]

    def x_t(self, x0, t, t0, t1, t2, pv0, pr0, thrust, min_t, isp, m0):
        """
        calculate state vector at given time t
        :param x0: initial state vector
        :param t: final time
        :param t0: initial time
        :param t1: first thrust switch time, which we turn thrust to T_min(thrust * min_t)
        :param t2: second thrust switch time, which we turn thrust back
        :param pv0: initial Pv
        :param pr0: initial Pr
        :param thrust: vessel's max thrust
        :param min_t: min throttle
        :param isp: specific impulse
        :param m0: initial mass
        :return: the corresponding satate vector
        """

        def i(t, t0, pv, pr, thrust, isp, m0):
            i_c = self.thrust_integral('c', t, t0, pv, pr, thrust, isp, m0)
            i_s = self.thrust_integral('s', t, t0, pv, pr, thrust, isp, m0)
            return np.vstack((i_c, i_s))

        if t1 > t2:
            print("warning: t1>t2, t1=" + str(t1) + " t2=" + str(t2))
            t2 = t1
        if t0 < t1:
            # x at t1
            i_t1 = i(t1, t0, pv0, pr0, thrust, isp, m0)
            x_temp = self.transfer_mat(t1, t0) @ x0 + self.gamma_mat(t1, t0) @ i_t1
            pv_t1, pr_t1 = self.lambda_t(t1, t0, pv0, pr0)
            m_t1 = self.mass_t(t1, t0, isp, m0, thrust)
            # x at t2
            i_t2 = i(t2, t1, pv_t1, pr_t1, thrust * min_t, isp, m_t1)
            x_temp = self.transfer_mat(t2, t1) @ x_temp + self.gamma_mat(t2, t1) @ i_t2
            pv_t2, pr_t2 = self.lambda_t(t2, t1, pv_t1, pr_t1)
            m_t2 = self.mass_t(t2, t1, isp, m_t1, thrust * min_t)
            # x at t
            i_t = i(t, t2, pv_t2, pr_t2, thrust, isp, m_t2)
            self.m_f = self.mass_t(t, t2, isp, m_t2, thrust)
            return self.transfer_mat(t, t2) @ x_temp + self.gamma_mat(t, t2) @ i_t
        elif t0 > t1:
            if t0 < t2:
                # x at t2
                i_t2 = i(t2, t0, pv0, pr0, thrust * min_t, isp, m0)
                pv_t2, pr_t2 = self.lambda_t(t2, t0, pv0, pr0)
                m_t2 = self.mass_t(t2, t0, isp, m0, thrust * min_t)
                x_temp = self.transfer_mat(t2, t0) @ x0 + self.gamma_mat(t2, t0) @ i_t2
                # x at t
                i_t = i(t, t2, pv_t2, pr_t2, thrust, isp, m_t2)
                self.m_f = self.mass_t(t, t2, isp, m_t2, thrust)
                return self.transfer_mat(t, t2) @ x_temp + self.gamma_mat(t, t2) @ i_t
            elif t0 > t2:
                i_t = i(t, t0, pv0, pr0, thrust, isp, m0)
                self.m_f = self.mass_t(t, t0, isp, m0, thrust)
                return self.transfer_mat(t, t0) @ x0 + self.gamma_mat(t, t0) @ i_t

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
            pv, _ = self.lambda_t(t, t0, pv0, pr0)
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

    def target_function_sl(self, z, r_f, v_f, x_0, t_0, t_1, t_2, thrust, min_t, isp, m0):
        pv = z[0:3].reshape(-1, 1)
        pr = z[3:6].reshape(-1, 1)
        tf = z[6]
        # final position constraint
        self.x_f = self.x_t(x_0, tf, t_0, t_1, t_2, pv, pr, thrust, min_t, isp, m0)
        r_tf = self.x_f[0:3, :]
        s1 = r_tf.T @ r_tf - r_f.T @ r_f
        # final velocity constraint
        v_tf = self.x_f[3:6, :]
        s2 = v_tf - v_f  # s2 includes 3 equations
        # transversality condition
        pvf, prf = self.lambda_t(tf, t_0, pv, pr)

        s3 = prf.T @ v_tf - pvf.T @ r_tf + \
             np.linalg.norm(pvf) * thrust / (self.mass_t(tf, t_0, isp, m0, thrust) * self.g0) - 1
        s4 = r_tf[2, 0] * prf[0, 0] - r_tf[0, 0] * prf[2, 0]
        s5 = r_tf[2, 0] * prf[1, 0] - r_tf[1, 0] * prf[2, 0]
        return [s1[0, 0], s2[0, 0], s2[1, 0], s2[2, 0], s3[0, 0], s4, s5]

    def target_function_t2_sl(self, t2, z0, r_f, v_f, x_0, t_0, t_1, thrust, min_t, isp, m0):
        # m_f is calculated during root finding
        root(
            self.target_function_sl, z0, args=(r_f, v_f, x_0, t_0, t_1, t2, thrust, min_t, isp, m0),
            method='lm', jac=False)
        '''z = sol.x.reshape(-1, 1)
        tf = z[6, 0]
        m_t1 = self.mass_t(t_1, t_0, isp, m0, thrust)
        m_t2 = self.mass_t(t2, t_1, isp, m_t1, thrust * min_t)
        m_tf = self.mass_t(tf, t2, isp, m_t2, thrust)'''
        return m0 - self.m_f

    def v_go(self, isp, m0):
        v_exh = isp * ge
        return v_exh*np.log(m0/self.m_f)

    def update(self, velocity_stream, position_stream, mass_stream, max_thrust_stream, ut_stream):
        v_0 = np.asarray(velocity_stream()).reshape(-1, 1) / self.velocity_multiplier
        r_0 = np.asarray(position_stream()).reshape(-1, 1) / self.position_multiplier
        x = np.vstack((r_0, v_0))
        m0 = mass_stream()
        thrust = max_thrust_stream()
        t_0 = ut_stream() / self.time_multiplier
        return x, m0, thrust, t_0
