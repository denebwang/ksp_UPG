import numpy as np

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
                return self.transfer_mat(t, t2) @ x_temp + self.gamma_mat(t, t2) @ i_t
            elif t0 > t2:
                i_t = i(t, t0, pv0, pr0, thrust, isp, m0)
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
