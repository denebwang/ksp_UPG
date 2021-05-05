import numpy as np


def fun_at(at, v0, sin_gamma, h0, g0):
    return at ** 2 + ((v0 ** 2 * sin_gamma) / (2 * h0)) * at - \
           (((v0 ** 2) * g0 * (1 + sin_gamma ** 2) / (4 * h0)) + g0 ** 2)


def t_togo(v0, sin_gamma, at, g0):
    return 1.15 * (v0 / 2) * (((1 + sin_gamma) / (at + g0)) + ((1 - sin_gamma) / (at - g0)))

class APDG(object):

    def __init__(self, r_f, v_f, a_f, t_f, t_0, m_0, r_0, v_0):
        self.r_f = r_f
        self.v_f = v_f
        self.t_f = t_f
        self.a_f = a_f

        self.t_0 = t_0
        # command values
        self.thrust_direction = np.zeros((3, 1))
        self.thrust = 0.
        # variables
        self.m_0 = m_0
        self.r_0 = r_0
        self.v_0 = v_0
        self.t_go = 0.

    def update(self, r_0, v_0, m_0, t_0):
        self.r_0 = r_0
        self.v_0 = v_0
        self.m_0 = m_0
        self.t_0 = t_0

    def compute(self):
        self.t_go = self.t_f - self.t_0
        a_t = \
            (-6 / self.t_go) * (self.v_f - self.v_0) + \
            (12 / self.t_go ** 2) * (self.r_f - self.r_0 - self.v_0 * self.t_go) + \
            self.a_f
        a = np.linalg.norm(a_t)
        self.thrust = a * self.m_0
        self.thrust_direction = a_t / a
