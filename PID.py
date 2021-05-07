# -*- coding: utf-8 -*-
import time


class PID(object):

    def __init__(self, Kp=1.0, Ki=0.1, Kd=0.01):
        self.Kp = Kp    # P controls reaction to the instantaneous error
        self.Ki = Ki    # I controls reaction to the history of error
        self.Kd = Kd   # D prevents overshoot by considering rate of change
        self.P = 0.0
        self.I = 0.0
        self.D = 0.0
        self.SetPoint = 0.0  # Target value for controller
        self.ClampI = 20.0  # clamps i_term to prevent 'windup.'
        self.LastTime = time.time()
        self.LastError = 0.0

    def update(self, measure):
        now = time.time()
        dt = now - self.LastTime
        if not dt:
            dt = 1.0   # avoid potential divide by zero if PID just created.

        error = self.SetPoint - measure
        self.P = error
        self.I += error
        self.I = self.clamp_i(self.I)   # clamp to prevent windup lag
        self.D = (error - self.LastError) / dt

        self.LastMeasure = measure  # store data for next update
        self.lastTime = now

        return (self.Kp * self.P) + (self.Ki * self.I) + (self.Kd * self.D)

    def clamp_i(self, i):
        if i > self.ClampI:
            return self.ClampI
        elif i < -self.ClampI:
            return -self.ClampI
        else:
            return i

    def setpoint(self, value):
        self.SetPoint = value
        self.I = 0.0
