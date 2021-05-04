
def fun_at(at, v0, sin_gamma, h0, g0):
    return at ** 2 + ((v0 ** 2 * sin_gamma) / (2 * h0)) * at - \
           (((v0 ** 2) * g0 * (1 + sin_gamma ** 2) / (4 * h0)) + g0 ** 2)


def t_togo(v0, sin_gamma, at, g0):
    return 1.15 * (v0 / 2) * (((1 + sin_gamma) / (at + g0)) + ((1 - sin_gamma) / (at - g0)))
