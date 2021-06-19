#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os.path

filename1 = "out1.txt"
if not os.path.isfile(filename1):
    print('File does not exist.')
else:
    f1 = open(filename1)
        


positions1 = [data.strip() for line in f1.readlines() for data in line.split(',') if data.strip()]
rx1 = []
ry1 = []
rx2 = []
ry2 = []
rxp = []
ryp = []
i = 0
while (i < len(positions1)):
    rx1.append(float(positions1[i]))
    ry1.append(float(positions1[i+1]))
    rx2.append(float(positions1[i+2]))
    ry2.append(float(positions1[i+3]))
    rxp.append(float(positions1[i+4]))
    ryp.append(float(positions1[i+5]))
    i += 6

fig1, ax1 = plt.subplots()
ax1.set_box_aspect(1)
ax1.scatter(rx1,ry1,s=100,c='#FD2E24')
ax1.scatter(rx2,ry2,s=5,c='#2E1A81')
ax1.scatter(rxp,ryp,s=1,c='black')
plt.xlabel('x')
plt.ylabel('y')
ax1.legend(['sun','earth','particle'])

plt.show()


filename2 = "out2.txt"
if not os.path.isfile(filename2):
    print('File does not exist.')
else:
    f2 = open(filename2)
        


positions2 = [data.strip() for line in f2.readlines() for data in line.split(',') if data.strip()]
rxCM =0
ryCM =0
uxp = []
uyp = []
#print(positions2)
i = 0
while (i < len(positions2)):
    if i == 0:
        ux1 = float(positions2[i])
        uy1 = float(positions2[i+1])
        ux2 = float(positions2[i+2])
        uy2 = float(positions2[i+3])
        i += 4
    else:
        uxp.append(float(positions2[i]))
        uyp.append(float(positions2[i+1]))
        i += 2

fig2, ax2 = plt.subplots()
ax2.set_box_aspect(1)
ax2.scatter(rxCM,ryCM,s=10,c='black')
ax2.scatter(ux1,uy1,s=100,c='#FD2E24')
ax2.scatter(ux2,uy2,s=5,c='#2E1A81')
ax2.scatter(uxp,uyp,s=1,c='black')
plt.xlabel('x')
plt.ylabel('y')
ax2.legend(['center of mass','sun','earth','particle'])

plt.show()