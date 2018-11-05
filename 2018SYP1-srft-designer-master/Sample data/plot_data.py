# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 10:54:22 2013

@author: rob
"""

from pylab import *

fwm,zrwm,ziwm=loadtxt("with_matching.txt",unpack=True)
fnm,zrnm,zinm=loadtxt("no_matching.txt",unpack=True)
fgp,zrgp,zigp=loadtxt("green_probe.txt",unpack=True)

znm=zrnm+1j*zinm
zm=zrwm+1j*ziwm

figure(1)
subplot(211)
plot(fwm,zrwm,label="With Matching")
plot(fnm,zrnm,label="No Matching")
plot(fgp,zrgp,label="Green Probe")
legend()
subplot(212)
plot(fwm,ziwm,label="With Matching")
plot(fnm,zinm,label="No Matching")
plot(fgp,zigp,label="Green Probe")
legend()


R0=50.0
tpg_nm=4*R0*znm/abs(R0+znm)**2
tpg_m=4*R0*zm/abs(R0+zm)**2

figure(2)
plot(fnm,tpg_nm, label="No Match")
plot(fwm,tpg_m, label="Match")
show()
