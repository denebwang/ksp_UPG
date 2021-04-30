import time

import krpc
import json
import numpy as np
from scipy import optimize

conn = krpc.connect()
space_center = conn.space_center
vessel = space_center.active_vessel
ap = vessel.auto_pilot
body = vessel.orbit.body
body_frame = body.reference_frame
flight = vessel.flight(body_frame)
control = vessel.control
r0 = body.equatorial_radius
g0 = body.surface_gravity
ge = 9.80665  # standard gravity which is used to calculate exhaust velocity

# soft landing prob with full T
with open("target.json") as f:
    target = json.load(f)
# target vector
lon = target['lon']
lat = target['lat']
v = target['vf']
r_f = np.asarray(body.position_at_altitude(lat, lon,
                                           body.surface_height(lat, lon) + 100, body_frame)).reshape(-1, 1)

temp1 = space_center.ReferenceFrame.create_relative(body_frame,
                                                    rotation=(0., np.sin(-lon / 2. * np.pi / 180), 0.,
                                                              np.cos(-lon / 2. * np.pi / 180)))
temp2 = space_center.ReferenceFrame.create_relative(temp1,
                                                    rotation=(0., 0., np.sin(lat / 2. * np.pi / 180),
                                                              np.cos(lat / 2. * np.pi / 180)))
height = body.surface_height(lat, lon) + body.equatorial_radius
target_frame = space_center.ReferenceFrame.create_relative(temp2, position=(height, 0., 0.))
del temp1, temp2
'''line = conn.drawing.add_line((10,0,0),(0,0,0), target_frame)
line.color = (1,0,0)
line.thickness = 2'''
target_v = np.array([-v, 0., 0.])
v_f = np.asarray(space_center.transform_velocity((0., 0., 0.), target_v, target_frame, body_frame)).reshape(-1, 1)

# normalize to dimensionless
r_f = r_f / r0
v_f = v_f / np.sqrt(r0 * g0)


def mass_t(t, isp, m0, thrust):
    """
    calculate the expected mass at dimensionless time t
    :param t: dimensionless time
    :param isp: engine's specific impulse
    :param m0: initial mass
    :param thrust: engine thrust
    :return: mass at t
    """
    return m0 - t * np.sqrt(r0 / g0) * thrust / (isp * ge)


def transfer_mat(t, t0):
    return np.block([[np.cos(t - t0) * np.eye(3), np.sin(t - t0) * np.eye(3)],
                     [-np.sin(t - t0) * np.eye(3), np.cos(t - t0) * np.eye(3)]])


def gamma_mat(t, t0):
    return np.block([[np.sin(t - t0) * np.eye(3), -np.cos(t - t0) * np.eye(3)],
                     [np.cos(t - t0) * np.eye(3), np.sin(t - t0) * np.eye(3)]])


def lambda_t(t, t0, pv0, pr0):
    """
    calculates the costate vectors at t(dimensionless)
    :param t: desired t
    :param t0: initial t
    :param pr0: initial Pr, 3x1
    :param pv0: initial Pv, 3x1
    :return: Pv,Pr at t
    """
    lambda_0 = np.vstack((pv0, -pr0))
    return transfer_mat(t, t0) @ lambda_0


def x_t(x0, t, t0, pv0, pr0, thrust, isp, m0):
    i_c = thrust_integral('c', t, t0, pv0, pr0, thrust, isp, m0)
    i_s = thrust_integral('s', t, t0, pv0, pr0, thrust, isp, m0)
    i = np.vstack((i_c, i_s))
    #print(i)
    return transfer_mat(t, t0) @ x0 + gamma_mat(t, t0) @ i


def thrust_integral(type_name, t, t0, pv0, pr0, thrust, isp, m0):
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
        pv = lambda_t(t, t0, pv0, pr0)[0:3, :]
        pv = pv / np.linalg.norm(pv)
        temp = pv * thrust / (mass_t(t, isp, m0, thrust) * g0)
        if type_name == 'c':
            return temp * np.cos(t - t0)
        elif type_name == 's':
            return temp * np.sin(t - t0)
        else:
            raise Exception()

    delta = (t - t0) / 4.
    return ((t - t0) / 90.) * (7 * i(t0) + 32 * i(t0 + delta) + 12 * i(t0 + 2 * delta) + 32 * i(t0 + 3 * delta) +
                               7 * i(t0 + 4 * delta))


# get some useful values
specific_impulse = vessel.vacuum_specific_impulse
max_thrust = vessel.max_thrust
mass = conn.add_stream(getattr, vessel, 'mass')
velocity = conn.add_stream(vessel.velocity, body_frame)
position = conn.add_stream(vessel.position, body_frame)
m0 = mass()
v_0 = np.asarray(velocity()).reshape(-1, 1) / np.sqrt(r0 * g0)
r_0 = np.asarray(position()).reshape(-1, 1) / r0

# initial guess

pv = np.asarray(-v_0 / np.linalg.norm(v_0)).reshape(-1, 1)
pr = np.zeros((3, 1))
tf = np.array([300. / np.sqrt(r0 / g0)]).reshape(1, 1)
x = np.vstack((r_0, v_0))

z0 = np.hstack((pv.T, pr.T, tf)).T


def fun(z):
    pv = z[0:3].reshape(-1, 1)
    pr = z[3:6].reshape(-1, 1)
    tf = z[6]
    # final position constraint
    x_f = x_t(x, tf, 0, pv, pr, max_thrust, specific_impulse, m0)
    r_tf = x_f[0:3, :]
    s1 = r_tf.T @ r_tf - r_f.T @ r_f
    # final velocity constraint
    v_tf = x_f[3:6, :]
    s2 = v_tf - v_f  # s2 includes 3 equations
    # transversality condition
    lambdat = lambda_t(tf, 0, pv, pr)
    pvf = lambdat[0:3, :]
    prf = lambdat[3:6, :]
    s3 = prf.T @ v_tf - pvf.T @ r_tf + \
         np.linalg.norm(pvf) * max_thrust / (mass_t(tf, specific_impulse, m0, max_thrust) * g0) - 1
    s4 = r_tf[2, 0] * prf[0, 0] - r_tf[0, 0] * prf[2, 0]
    s5 = r_tf[2, 0] * prf[1, 0] - r_tf[1, 0] * prf[2, 0]
    return [s1[0, 0], s2[0, 0], s2[1, 0], s2[2, 0], s3[0, 0], s4, s5]


ap.reference_frame = body_frame
ap.engage()

while True:
    sol = optimize.root(fun, z0, method='lm', jac=False)
    z = sol.x.reshape(-1, 1)
    z0 = z
    pv = z[0:3, :]
    pr = z[3:6, :]
    tf = z[6, 0]
    tgo =tf * np.sqrt(r0 / g0)
    distance = np.linalg.norm((r_0 - r_f) * r0)
    print("tgo:" + str(tgo))
    #print("distance:" + str(distance))
    unit_t = pv / np.linalg.norm(pv)
    ap.target_direction = tuple(unit_t.reshape(1, -1)[0])
    if tgo < 0.5 and flight.surface_altitude < 2000:
        control.throttle = 0
        break
    time.sleep(0.5)
    # get new values for next loop
    v_0 = np.asarray(velocity()).reshape(-1, 1) / np.sqrt(r0 * g0)
    r_0 = np.asarray(position()).reshape(-1, 1) / r0
    x = np.vstack((r_0, v_0))
    m0 = mass()
ap.disengage()
