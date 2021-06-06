#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# constants and masses
Ms = 1.988500 * 10e30 #Solar mass [kg]
AU = 149597870700 #Astronomic Unit [m]
G = (6.67430 * 10e-11)*(1/(AU**3))*Ms*(60*60*24)**2 #Gravitational const. in AU^3 * SolarMass * Day^-2
Me = (5.9724 * 10e26)/Ms #Mass of the earth as a ration of Solar mass
Mo = 0/Ms # Mass of arbitrary object and Mo << Me,Ms [kg]


class Object:                   # define the Sun, Earth and the object
    def __init__(self, name, mass, r0, v0):
        self.name = name
        self.mass = mass
        self.r = [r0]
        self.v = [v0]

class System:
    def __init__(self, objects):
        self.objects = objects
        self.time = [0]

    def update(self):
        dt = 1
        for p in self.objects:
            acc = G*
            p.v += acc * dt
            p.xs.append(p.r[0])
            p.ys.append(p.r[1])
            p.plot.set_offsets(p.r[:2])
            p.line.set_xdata(p.xs)
            p.line.set_ydata(p.ys)
        return plots + lines
                





