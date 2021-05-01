
class UPG(object):
    ge = 9.80665  # standard gravity which is used to calculate exhaust velocity

    def __init__(self, conn):
        self.conn = conn
        space_center = conn.space_center
        vessel = space_center.active_vessel
        ap = vessel.auto_pilot
        body = vessel.orbit.body
        body_frame = body.reference_frame
        flight = vessel.flight(body_frame)
        r0 = body.equatorial_radius
        g0 = body.surface_gravity

