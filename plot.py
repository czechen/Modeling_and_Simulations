#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import os.path

# Reading the first output file
filename1 = "out1.txt"
if not os.path.isfile(filename1):
    print('File does not exist.')
else:
    f1 = open(filename1)
        

# Putting all the lines in a vector
positions1 = [data.strip() for line in f1.readlines() for data in line.split(',') if data.strip()]
rx1 = []
ry1 = []
rx2 = []
ry2 = []
rxp = []
ryp = []
i = 0
# Reading all the values from the postion vector
while (i < len(positions1)):
    rx1.append(float(positions1[i]))
    ry1.append(float(positions1[i+1]))
    rx2.append(float(positions1[i+2]))
    ry2.append(float(positions1[i+3]))
    rxp.append(float(positions1[i+4]))
    ryp.append(float(positions1[i+5]))
    i += 6
# Plot
fig1, ax1 = plt.subplots()
ax1.set_box_aspect(1)
ax1.scatter(rx1,ry1,s=100,c='#FD2E24')
ax1.scatter(rx2,ry2,s=5,c='#2E1A81')
ax1.scatter(rxp,ryp,s=1,c='black')
plt.xlabel('x')
plt.ylabel('y')
ax1.legend(['sun','earth','particle'])

plt.show()

# For the second and the third output file it is the same process

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
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
ax2.legend(['center of mass','sun','earth','particle'])

plt.show()


filename3 = "out3.txt"
if not os.path.isfile(filename3):
    print('File does not exist.')
else:
    f3 = open(filename3)

potential3 = [data.strip() for line in f3.readlines() for data in line.split(',') if data.strip()]
omega = []
x = []
y = []
j = 0
while(j<len(potential3)):
    omega.append(float(potential3[j]))
    x.append(float(potential3[j+1]))
    y.append(float(potential3[j+2]))
    j +=3

fig3, ax3 = plt.subplots(subplot_kw={"projection": "3d"})

k = int(np.sqrt(len(x)))
X = np.reshape(x, (k, k))
Y = np.reshape(y, (k, k))
Omega = np.reshape(omega, (k, k))


#3D plot of the potential
surf = ax3.plot_surface(X,Y,Omega,cmap=plt.cm.jet,linewidth=0, antialiased=False)
plt.show()

#Contour plot of the potential
fig4,ax4=plt.subplots()
c = np.linspace(-7, -2, num=50)
plt.contourf(X,Y,Omega,c,cmap=plt.cm.jet)
ax4.scatter(rxCM,ryCM,s=10,c='black')
ax4.scatter(ux1,uy1,s=100,c='#FD2E24')
ax4.scatter(ux2,uy2,s=5,c='#2E1A81')
ax4.scatter(uxp,uyp,s=1,c='black')
ax4.set_box_aspect(1)
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()
plt.show()

