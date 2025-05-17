import openmc
import numpy as np
import matplotlib.pyplot as plt

import matplotlib


import scipy.interpolate as scint


matplotlib.use('webagg')
def inches_to_cm(inches):
    return inches * 2.54



class InnerGeometry:
    def __init__(self,rBall,rOut,rIn,hTot,hCone1,theta1 = np.pi/2,theta2=np.pi/180*30,Inches = True): #Not gonna implement theta1 =/= np.pi/2
        if Inches: 
            rBall = inches_to_cm(rBall)
            rOut = inches_to_cm(rOut)
            rIn = inches_to_cm(rIn)
            hTot = inches_to_cm(hTot)
            hCone1 = inches_to_cm(hCone1)

        self.rBall = rBall
        self.rOut = rOut
        self.rIn = rIn
        self.hTot = hTot
        self.hCone1 = hCone1
        self.theta1 = theta1
        self.theta2 = theta2
    def getConeHeight(self,r,R,theta):
        return (R-r)/np.tan(theta)
    def getConeRadius(self,h,R,theta):
        return R-h*np.tan(theta)
    
    def totRegion(self):
        ball = +openmc.Sphere(r=self.rBall)
        topCylinder = -openmc.Cylinder(r=self.rOut, dz=self.hTot/2)
        bottomCylinder = -openmc.Cylinder(r=self.rIn, z0=-self.hTot/2, dz=self.hTot/2)

        cone1Z0 = self.rBall/np.sqrt(2)
        cone1R = self.rBall/np.sqrt(2)
        cone1h = cone1R/np.tan(self.theta1)
        cone1 = -openmc.Cone(r2=np.power(cone1R/cone1h,2), z0=-cone1Z0, dz=-cone1h)

        cone2Z0 = cone1Z0+self.getConeHeight(self.rIn,cone1R,self.theta1)
        cone2R = self.getConeRadius(self.hCone1,cone1R,self.theta2)
        cone2h = cone2R/np.tan(self.theta2)

        cone2 = -openmc.Cone(r2=np.power(cone2R/cone2h,2), z0=-cone2Z0, dz=-cone2h)

        return (ball | topCylinder) | (bottomCylinder | cone1) | cone2



        
"""
#Spherical Geometry

#Fuelball, zirc, blanket, steel clad, carbon steel, air, blast shield 
materialNames = ['fuel', 'zirc', 'blanket', 'steelClad', 'carbSteel', 'air', 'blastShield']
ballDiameters = inches_to_cm(np.asarray([32,32+10/16,60,60+0.4*2,60+0.4*2+4.4*2,74,74+2*(1+5/8)]))
ballRadii = ballDiameters/2
spheres = []
for i in range(len(ballDiameters)):
    spheres.append(openmc.Sphere(r=ballRadii[i]))


#Cylindrical Geometry

cylinderTopRadii = inches_to_cm(np.asarray([4,4+5/16,4+5/16+0.4]))
topCylinders = []
for i in range(len(cylinderTopRadii)):
    topCylinders.append(openmc.ZCylinder(r=cylinderTopRadii[i], z0=inches_to_cm(12*5+10)))


#fuel geometry
fuelRegion = openmc.Sphere()
"""

coreInlet = inches_to_cm(3.5/2)#Think they were referencing diameter
coreOutlet = inches_to_cm(4/2)
coreBallPlateThickness = inches_to_cm(5/16)
coreConePlateThickness = inches_to_cm(3/8) #Gonna need to assume this is the same as the ball plate thickness

hTot = inches_to_cm(11*12+6+3/16)
startPlate = inches_to_cm(7*12+9+1/4) - inches_to_cm(6*12+10+3/8)



plateSpacingsCone1 = np.asarray([1+47/64,1+23/32,1+3/8,1+1/8])
coneh1 = inches_to_cm(np.sum(plateSpacingsCone1)) + startPlate


fuelGeometry = InnerGeometry(16,coreOutlet,coreInlet,hTot,coneh1)
fuelRegion = fuelGeometry.totRegion()

D2O = openmc.Material(name='heavy water')
D2O.set_density('g/cm3', 1.1056)
# use nuclide H2 for deuterium and O16 for oxygen  
D2O.add_nuclide('H2', 2.0, 'ao')
D2O.add_nuclide('O16', 1.0, 'ao')
cell_fuel      = openmc.Cell(name='fuel', fill=D2O, region=fuelRegion)
root_univ = openmc.Universe(cells=[cell_fuel])



vox_plot = openmc.Plot()

vox_plot.type = 'voxel'
vox_plot.width = (100., 100., 50.)
vox_plot.pixels = (400, 400, 200)

plots = openmc.Plots([vox_plot])
plots.export_to_xml("testPlot.xml")

openmc.plot_geometry()
root_univ.plot()
plt.savefig("git/reaktorProject/Mads/testfig3")