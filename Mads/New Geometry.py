import openmc
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


import scipy.interpolate as scint


matplotlib.use('webagg')
def inches_to_cm(inches):
    return inches * 2.54

def pythagorasSide(a,c):
    return np.sqrt(c**2 - a**2)
def pythagorasHypotenuse(a,b):
    return np.sqrt(a**2 + b**2)


#All the balls.
fuel_Ball_outer = openmc.Sphere(r = inches_to_cm(16))
#fuel_Ball_outer = openmc.Sphere(r = inches_to_cm(16 + 5/16))
zirc_ball_inner = fuel_Ball_outer
zirc_Ball_outer = openmc.Sphere(r = inches_to_cm(16 + 5/16))

blanket_ball_inner = zirc_Ball_outer
blanket_Ball_outer = openmc.Sphere(r = inches_to_cm(30))

steelClad_Ball_inner = blanket_Ball_outer
steelClad_Ball_outer = openmc.Sphere(r = inches_to_cm(30+0.4))

carbSteel_Ball_inner = steelClad_Ball_outer
carbSteel_Ball_outer = openmc.Sphere(r = inches_to_cm(30.4+4.4))


#Cylindrical geometry, top
fuel_topCyl_outer = openmc.ZCylinder(0,inches_to_cm(12*5+10),r=inches_to_cm(4))

zirc_topCyl_inner = fuel_topCyl_outer
zirc_topCyl_outer = openmc.ZCylinder(inches_to_cm(16),inches_to_cm(12*5+10),r=inches_to_cm(4 + 5/16))

steelClad_topCyl_inner = zirc_topCyl_outer
steelClad_topCyl_outer = openmc.ZCylinder(inches_to_cm(30),inches_to_cm(12*5+10),r=inches_to_cm(4 + 5/16 + 0.4)) #Vacuum tag should be added

#Cylindrical geometry, bottom




plateSpacingsTopCone = [1+47/64,1+23/32,1+3/8,1+1/8]
plateThicknessesTopCone = [1/8,1/8,1/8,1/8] #Might have an off by one error, check this
plateSpacingsBottomCone = [2+21/32,2+51/64,2+3/8,1+13/16]


cenToTopCone = inches_to_cm((11+1/16))
topConeHeight = inches_to_cm(np.sum(plateThicknessesTopCone)+np.sum(plateSpacingsTopCone))
topConeRadius = inches_to_cm(pythagorasSide(11+1/16,16))


#inches_to_cm(-np.sum(plateThicknessesTopCone)-np.sum(plateSpacingsTopCone))
fuel_topCone_outer = openmc.ZCylinder(cenToTopCone,-inches_to_cm(pythagorasSide(11+1/16,16))-inches_to_cm((11+1/16)),r=inches_to_cm(pythagorasSide(11+1/16,16))) #These mesurements aint' right (the height is bad)

zirc_topCone_inner = fuel_topCone_outer
zirc_topCone_outer = openmc.ZCylinder(inches_to_cm(-(11+1/16)),,r=inches_to_cm(pythagorasSide(11+1/16,16)))

fuel_bottomCone_outer = openmc.ZCylinder(inches_to_cm(-(11+1/16)),,r=inches_to_cm(pythagorasSide(11+1/16,16))) #These mesurements aint' right (the height is bad)


#carbSteel_topCyl_inner = steelClad_topCyl_outer
PV_inner = openmc.Sphere(r = inches_to_cm(30))              
PV_outer = openmc.Sphere(r = inches_to_cm(30 + 4.4), boundary_type = 'vacuum')      #base er apparently at alle lag er transmissive som boundary men kan ikke have transmissive boundary som yderste, skal være periodic, reflective eller vacuum    
