import numpy as np
import utils


def a_t(v_0: float, sin_gamma:float, h0:float, g:float):
    b = (v_0 ** 2) * sin_gamma / (2. * h0)
    c = (v_0 ** 2) * g * (1. + sin_gamma ** 2) / (4. * h0) - g ** 2
    delta = b ** 2 + 4. * c
    return (np.sqrt(delta) - b) / 2.


def t_togo(v_0, sin_gamma, a_t, g):
    return 1.2 * (v_0 / 2.) *\
    (((1. + sin_gamma) / (a_t + g)) + ((1. - sin_gamma) / (a_t - g)))


class APDG(object):

    def __init__(self, r_f, v_f, a_f_k, t_f, t_0, up ,m, r_0, v_0, g_0, k_r):
        self.r_f = r_f
        self.v_f = -v_f * up
        self.t_f = t_f
        self.g_0 = -g_0 * up
        self.up = up
        self.k_r = k_r
        
        k = (a_f_k - 1.) / (k_r / 6. - 1.) + 1.
        self.a_f = -k * self.g_0
        # command values
        self.a_t = np.zeros((3, 1))
        # variables
        self.m = m
        self.r_0 = r_0
        self.v_0 = v_0
        self.t_0 = t_0
        self.t_go = t_f - t_0

    def update(self, r_0, v_0, m_0, t_0):
        self.r_0 = r_0
        self.v_0 = v_0
        self.m = m_0
        self.t_go = self.t_f - t_0

    def compute(self):
        a_t = (2 * (1 - self.k_r / 3.) / self.t_go) * (self.v_f - self.v_0) +\
            (self.k_r / self.t_go ** 2) * (self.r_f - self.r_0 - self.v_0 * self.t_go) + \
            (self.k_r - 6.)/ 6. * self.a_f + (self.k_r - 12.)/ 6. * self.g_0
        self.a_t = a_t# - self.g_0

    @property
    def thrust_direction(self):
        unit_a_t = self.a_t / np.linalg.norm(self.a_t)
        return tuple(unit_a_t.flatten())

    @property
    def thrust(self):
        return np.linalg.norm(self.a_t) * self.m

