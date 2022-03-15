import krpc
import numpy as np
import sys


def vector(x):
    return np.asarray(x).reshape(-1, 1)


def target_reference_frame(space_center, body,  latitude: float, longitude: float):
    """Get the target frame at given latitude and longitude, with origin on the ground.

    Args:
        space_center (_type_): SpaceCenter object to create new frame
        body (_type_): The celestial body target is on
        latitude (float): Latitude of target
        longitude (float): Longitude of target

    Returns:
        _type_: The target reference frame with origin at the target point,
        x points up, y points to the north, and z points to the east
    """
    temp1 = space_center.ReferenceFrame.create_relative(
        body.reference_frame,
        rotation=(0., np.sin(-longitude / 2. * np.pi / 180), 0.,
                  np.cos(-longitude / 2. * np.pi / 180)))
    temp2 = space_center.ReferenceFrame.create_relative(
        temp1, rotation=(0., 0., np.sin(latitude / 2. * np.pi / 180),
                         np.cos(latitude / 2. * np.pi / 180)))
    height = body.surface_height(latitude, longitude) + body.equatorial_radius
    return space_center.ReferenceFrame.create_relative(
        temp2, position=(height, 0., 0.))


def thrust2throttle(thrust, max_thrust, min_throttle):
    return (thrust - max_thrust * min_throttle)\
        / (max_thrust * (1. - min_throttle))


def move_position2height(height, position, body, frame):
    position = tuple(position.flatten())
    lat = body.latitude_at_position(position, frame)
    lon = body.longitude_at_position(position, frame)
    return vector(body.position_at_altitude(
        lat, lon, body.surface_height(lat, lon) + height, frame))


def downrange(from_, to, body, frame):
    from_ = from_.flatten()
    to = to.flatten()
    theta = np.arccos(np.clip(np.dot(from_, to) / np.linalg.norm(from_) / np.linalg.norm(to), -1, 1))
    r = body.altitude_at_position(tuple(to.flatten()), frame) + body.equatorial_radius
    return r * theta
