#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os.path

filename = "out.txt"
if not os.path.isfile(filename):
    print('File does not exist.')
else:
    f = open(filename)
        


positions = [data.strip() for line in f.readlines() for data in line.split(',') if data.strip()]
rx1 = []
ry1 = []
rx2 = []
ry2 = []
rxp = []
ryp = []
i = 0
while (i < len(positions)):
    rx1.append(float(positions[i]))
    ry1.append(float(positions[i+1]))
    rx2.append(float(positions[i+2]))
    ry2.append(float(positions[i+3]))
    rxp.append(float(positions[i+4]))
    ryp.append(float(positions[i+5]))
    i += 6

fig1, ax = plt.subplots()
ax.set_box_aspect(1)
ax.scatter(rx1,ry1,s=100,c='#FD2E24')
ax.scatter(rx2,ry2,s=5,c='#2E1A81')
ax.scatter(rxp,ryp,s=1,c='black')
plt.xlabel('x')
plt.ylabel('y')
ax.legend(['sun','earth','particle'])

plt.show()

