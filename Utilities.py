import krpc
import numpy as np


def target_reference_frame(space_center, body,  latitude: int, longitude: int):
    """
    get the target frame at given latitude and longitude,on the ground.
    :param space_center: SpaceCenter object to create new frame
    :param body: the desired Celestial body
    :param latitude: latitude of target
    :param longitude: longitude of target
    :return: the target reference frame with origin at the target point,
    x points up, y points to the north, and z points to the east
    """
    temp1 = space_center.ReferenceFrame.create_relative(body.reference_frame,
                                                        rotation=(0., np.sin(-longitude / 2. * np.pi / 180), 0.,
                                                                  np.cos(-longitude / 2. * np.pi / 180)))
    temp2 = space_center.ReferenceFrame.create_relative(temp1,
                                                        rotation=(0., 0., np.sin(latitude / 2. * np.pi / 180),
                                                                  np.cos(latitude / 2. * np.pi / 180)))
    height = body.surface_height(latitude, longitude) + body.equatorial_radius
    return space_center.ReferenceFrame.create_relative(temp2, position=(height, 0., 0.))
