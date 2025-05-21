import openmc
import numpy as np
import matplotlib.pyplot as plt

import matplotlib


import scipy.interpolate as scint


matplotlib.use('webagg')
def inches_to_cm(inches):
    return inches * 2.54



class InnerGeometry:
    def __init__(self,rBall,rCone1,rOut,rIn,hTot,hTop,theta1 = np.pi/2,theta2=np.pi/180*30,Inches = True): #Not gonna implement theta1 =/= np.pi/2
        if Inches: 
            rBall = inches_to_cm(rBall)
            rOut = inches_to_cm(rOut)
            rIn = inches_to_cm(rIn)
            hTot = inches_to_cm(hTot)
            hTop = inches_to_cm(hTop)
            rCone1 = inches_to_cm(rCone1)

        self.rBall = rBall
        self.rOut = rOut
        self.rIn = rIn
        self.hTot = hTot
        self.hTop = hTop
        self.rCone1 = rCone1
        self.theta1 = theta1/2 #These are the half angles of the cones
        self.theta2 = theta2/2
    @property
    def hBottom(self):
        return self.hTot - self.hTop
    def getConeTotH(self,R,theta):
        return R/np.tan(theta)
    def getConeHeight(self,r,R,theta):
        return (R-r)/np.tan(theta)
    def getConeRadius(self,h,R,theta):
        return R-h*np.tan(theta)
    @property
    def RCone2(self):
        return self.rCone1
    @property
    def RCone1(self):
        return self.rBall/np.sqrt(2)
    @property
    def hCone1(self):
        return self.getConeHeight(self.rCone1,self.RCone1,self.theta1)
    @property
    def HCone1(self):
        return self.getConeTotH(self.RCone1,self.theta1)
    @property
    def HCone2(self):
        return self.getConeTotH(self.RCone2,self.theta2)
    @property
    def z0Cone1(self):
        return self.rBall*np.sqrt(2)
    @property
    def hCone2(self):
        return self.rCone1/np.tan(self.theta2)
    @property
    def z0Cone2(self):
        return self.z0Cone1/2+self.hCone1+self.hCone2


    def totRegion(self,BoundType= "transmission",outside=False):

        
        ball = -openmc.Sphere(r=self.rBall,boundary_type=BoundType)
        topCylinder = -openmc.Cylinder(r=self.rOut, dz=self.hTot/2,boundary_type=BoundType) &-openmc.ZPlane(z0=self.hTop,boundary_type=BoundType) &+openmc.ZPlane(z0=0,boundary_type=BoundType)
        bottomCylinder = -openmc.Cylinder(r=self.rIn, z0=-self.hTot/2, dz=self.hTot/2,boundary_type=BoundType) &-openmc.ZPlane(z0=0,boundary_type=BoundType)&+openmc.ZPlane(z0=-self.hBottom,boundary_type=BoundType)
        
        #print(cone1h)
        cone1 = -openmc.Cone(r2=np.power(self.RCone1/self.HCone1,2), z0=-self.z0Cone1, dz=-1,boundary_type=BoundType) &+openmc.ZPlane(z0=-self.z0Cone1, boundary_type=BoundType)&-openmc.ZPlane(z0=-self.z0Cone1/2, boundary_type=BoundType)
        print(self.z0Cone1,self.hCone1,self.rCone1,self.RCone1)
  
        
        #cone2h = cone2R/np.tan(self.theta2)
        #print(cone2h,cone2R)
        #Why are the boundaries different? Idk
        cone2 = -openmc.Cone(r2=np.power(self.RCone2/self.HCone2,2), z0=-self.z0Cone2, dz=-1,boundary_type=BoundType)&+openmc.ZPlane(z0=-self.z0Cone2, boundary_type=BoundType)&-openmc.ZPlane(z0=-self.z0Cone1/2, boundary_type=BoundType)
            
        return cone1 | ball |cone2 | topCylinder | bottomCylinder



        
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
rBall = inches_to_cm(16)
coreInlet = inches_to_cm(3.5/2)#Think they were referencing diameter
coreOutlet = inches_to_cm(4/2)
coreBallPlateThickness = inches_to_cm(5/16)
coreConePlateThickness = inches_to_cm(3/8) #Gonna need to assume this is the same as the ball plate thickness

hTot = inches_to_cm(11*12+6+3/16)
hTop = inches_to_cm(6*12+10+3/8)
#startPlate = inches_to_cm(7*12+9+1/4) - inches_to_cm(6*12+10+3/8)
rCone1 = inches_to_cm(4+15/16)


#plateSpacingsCone1 = np.asarray([1+47/64,1+23/32,1+3/8,1+1/8])
#coneh1 = inches_to_cm(np.sum(plateSpacingsCone1)) + startPlate


zircThickness = inches_to_cm(5/16)

"""
All of this is just to plot it.
"""
fuelGeometry = InnerGeometry(rBall,rCone1,coreOutlet,coreInlet,hTot,hTop,Inches=False)
zircGeometry = InnerGeometry(rBall+zircThickness,rCone1+zircThickness/np.cos(np.pi/4),coreOutlet+zircThickness,coreInlet+zircThickness,hTot,hTop,Inches=False)
fuelRegion = fuelGeometry.totRegion(outside=False)
zircRegion = zircGeometry.totRegion(outside=False)&(~fuelRegion)

rBlanket = inches_to_cm(30)
blanketGeometry = -openmc.Sphere(r=rBlanket)
blanketRegion = blanketGeometry & (~fuelRegion) & (~zircRegion)

rSteelClad = inches_to_cm(30+0.4)
steelCladGeometry = -openmc.Sphere(r=rSteelClad)
steelCladRegion = steelCladGeometry & (~blanketRegion) & (~fuelRegion) & (~zircRegion)

rCarbSteel = inches_to_cm(30.4+4.4)
carbSteelGeometry = -openmc.Sphere(r=rCarbSteel)
carbSteelRegion = carbSteelGeometry & (~blanketRegion) & (~steelCladRegion) & (~fuelRegion) & (~zircRegion)

rAir = inches_to_cm(74/2)
airGeometry = -openmc.Sphere(r=rAir)
airRegion = airGeometry & (~blanketRegion) & (~steelCladRegion) & (~carbSteelRegion) & (~fuelRegion) & (~zircRegion)
rblastShield = rAir+inches_to_cm((1+5/8))
blastShieldGeometry = -openmc.Sphere(r=rblastShield)
blastShieldRegion = blastShieldGeometry & (~airRegion) & (~blanketRegion) & (~steelCladRegion) & (~carbSteelRegion) & (~fuelRegion) & (~zircRegion)

D2O = openmc.Material(name='heavy water')
D2O.set_density('g/cm3', 1.1056)
# use nuclide H2 for deuterium and O16 for oxygen  
D2O.add_nuclide('H2', 2.0, 'ao')
D2O.add_nuclide('O16', 1.0, 'ao')

SS304 = openmc.Material(name='SS304')
SS304.set_density('g/cc', 8.03)
SS304.add_element('Si', 0.0060, 'wo')
SS304.add_element('Cr', 0.1900, 'wo')
SS304.add_element('Mn', 0.0200, 'wo')
SS304.add_element('Fe', 0.6840, 'wo')
SS304.add_element('Ni', 0.1000, 'wo')

#zirconium = openmc.Material("zirconium")
#zirconium.add_element('Zr', 1.0, 'wo')
#zirconium.add_element('Sn', 0.015, 'wo')
#zirconium.set_density('g/cm3', 6.6)


mats = openmc.Materials([D2O,SS304])
#mats.export_to_xml()
# Export the geometry to XML



cell_fuel      = openmc.Cell(name='fuel', fill=D2O, region=fuelRegion)
cell_zirc      = openmc.Cell(name='zirc', fill=SS304, region=zircRegion) #TODO Not the right material, but zirconium doesn't wanna play
cell_blanket   = openmc.Cell(name='blanket', fill=D2O, region=blanketRegion) 
cell_steelClad = openmc.Cell(name='steelClad', fill=SS304, region=steelCladRegion) #TODO Not quite the right material
cell_carbSteel = openmc.Cell(name='carbSteel', fill=SS304, region=carbSteelRegion) #TODO Material is wrong
cell_air       = openmc.Cell(name='air', fill=D2O, region=airRegion)
cell_blastShield = openmc.Cell(name='blastShield', fill=SS304, region=blastShieldRegion) #TODO Material is wrong

root_univ = openmc.Universe(cells=[cell_fuel,cell_zirc,cell_blanket,cell_steelClad,cell_carbSteel,cell_air,cell_blastShield])
geometry = openmc.Geometry(root_univ)

mats.export_to_xml()
geometry.export_to_xml()

plot = openmc.Plot()
plot.filename = 'xy_slice'
plot.width = (hTot, hTot)     # Width of plot in cm
plot.pixels = (4000, 4000)        # Image resolution
plot.basis = 'xz'               # Slice plane

plots = openmc.Plots([plot])
plots.export_to_xml()
openmc.plot_geometry()
#root_univ.plot()
#plt.savefig("git/reaktorProject/Mads/testfig3")