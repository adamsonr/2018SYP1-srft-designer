# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 10:54:22 2013

@author: rob
"""

#This script loads in an impedance file and displays
#the real and imaginary parts and the transducer power efficiency


from pylab import *
import easygui

R0=50

filename=easygui.fileopenbox(msg="Select original impedance file")
f,zr,zi=loadtxt(filename,unpack=True)
filename=easygui.fileopenbox(msg="Select matched impedance file")
f,zrm,zim=loadtxt(filename,unpack=True)


z=zr+1j*zi
zm=zrm+1j*zim

figure(1)
subplot(211)
plot(f/1e6,zr,label="Original Real")
plot(f/1e6,zi,label="Original Imag")
plot(f/1e6,zrm,label="Matched Real")
plot(f/1e6,zim,label="Matched Imag")
xlabel("Freq (MHz)")
ylabel("Z (Ohm)")
legend()
grid()
subplot(212)
#Transducer power gain
TPGm=4*R0*zrm/abs(R0+zm)**2  
TPG=4*R0*zr/abs(R0+z)**2  

plot(f,TPG,label="No Matching")
plot(f,TPGm,label="With Matching")
xlabel("Freq (MHz)")
ylabel("Transducer Power Gain")
legend()
grid()

show()
