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

cutPlane1 = openmc.YPlane(y0=0,boundary_type='periodic') #This is the plane that cuts the geometry in half, so we can see the inside

angle = np.pi/4

cutPlane2 = openmc.Plane(a=np.tan(angle),B=-1,D=0,boundary_type="periodic")
cutPlane1.periodic_surface = cutPlane2
cutPlane2.periodic_surface = cutPlane1
#slice = -cutPlane1 & -cutPlane2


fuelGeometry = InnerGeometry(rBall,rCone1,coreOutlet,coreInlet,hTot,hTop,Inches=False)
zircGeometry = InnerGeometry(rBall+zircThickness,rCone1+zircThickness/np.cos(np.pi/4),coreOutlet+zircThickness,coreInlet+zircThickness,hTot,hTop,Inches=False)
fuelRegion = fuelGeometry.totRegion(outside=False)
zircRegion = zircGeometry.totRegion(outside=False)&(~fuelRegion)#&slice

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
blastShieldGeometry = -openmc.Sphere(r=rblastShield,boundary_type='vacuum')
blastShieldRegion = blastShieldGeometry & (~airRegion) & (~blanketRegion) & (~steelCladRegion) & (~carbSteelRegion) & (~fuelRegion) & (~zircRegion)

vacuumBound = -openmc.ZCylinder(r=rblastShield+0.1,boundary_type='vacuum')&+openmc.ZPlane(z0=-fuelGeometry.hBottom-0.1,boundary_type='vacuum')&-openmc.ZPlane(z0=fuelGeometry.hTop+0.1,boundary_type='vacuum') #Big cylinder of vacuum around the whole thing
#vacuumRegion = vacuumBound & (~airRegion) & (~blanketRegion) & (~steelCladRegion) & (~carbSteelRegion) & (~fuelRegion) & (~zircRegion) #Not sure if this is needed or not.
slice = blastShieldGeometry

D2O = openmc.Material(name='heavy water')
D2O.set_density('g/cm3', 1.1056)
# use nuclide H2 for deuterium and O16 for oxygen  
D2O.add_nuclide('H2', 2.0, 'ao')
D2O.add_nuclide('O16', 1.0, 'ao')

#stainless steel
SS304 = openmc.Material(name='SS304')
SS304.set_density('g/cc', 8.03)
SS304.add_element('Si', 0.0060, 'wo')
SS304.add_element('Cr', 0.1900, 'wo')
SS304.add_element('Mn', 0.0200, 'wo')
SS304.add_element('Fe', 0.6840, 'wo')
SS304.add_element('Ni', 0.1000, 'wo')


#Carbon steel, this was available when the reactor was built according to chatgpt so don't trust it too much (it doesn't wanna give me a source though)
#https://www.azom.com/article.aspx?ArticleID=6772
S4340 = openmc.Material(name='S4340')
S4340.set_density('g/cc', 7.85)
S4340.add_element('C', 0.004, 'wo')
S4340.add_element('Si', 0.002, 'wo')
S4340.add_element('Mn', 0.008, 'wo')
S4340.add_element('P', 0.00035, 'wo')
S4340.add_element('S', 0.0004, 'wo')
S4340.add_element('Cr', 0.008, 'wo')
S4340.add_element('Mo', 0.0025, 'wo')
S4340.add_element('Ni', 0.02, 'wo')
S4340.add_element('Fe', 1-0.04525, 'wo') #just manually calculated this, should do something else


zirconium = openmc.Material(name="zirconium")
zirconium.set_density('g/cm3', 6.6)
zirconium.add_element('Zr', 1.0, 'wo')
zirconium.add_element('Sn', 0.015, 'wo')

air = openmc.Material(name='air')
air.set_density('g/cm3', 1.205e-3) #density at room temp
air.add_element('N', 0.7808, 'wo')
air.add_element('O', 0.2095, 'wo')
air.add_element('Ar', 0.0093, 'wo')
#Ignore the rest.
#air.add_elements_from_formula('CO2', 0.0004, 'wo')




vacuum = openmc.Material(name='vacuum')
vacuum.set_density('g/cm3', 1e-50)
vacuum.add_element('N', 0.7808, 'wo')
vacuum.add_element('O', 0.2095, 'wo')
vacuum.add_element('Ar', 0.0093, 'wo')


mats = openmc.Materials([D2O,SS304,zirconium,S4340,air,vacuum])
#mats.export_to_xml()
# Export the geometry to XML



cell_fuel      = openmc.Cell(name='fuel', fill=D2O, region=fuelRegion&slice)
cell_zirc      = openmc.Cell(name='zirc', fill=zirconium, region=zircRegion&slice)
cell_blanket   = openmc.Cell(name='blanket', fill=D2O, region=blanketRegion&slice) 
cell_steelClad = openmc.Cell(name='steelClad', fill=SS304, region=steelCladRegion&slice) 
cell_carbSteel = openmc.Cell(name='carbSteel', fill=S4340, region=carbSteelRegion&slice) #TODO Material is wrong
cell_air       = openmc.Cell(name='air', fill=air, region=airRegion&slice)
cell_blastShield = openmc.Cell(name='blastShield', fill=S4340, region=blastShieldRegion&slice) #TODO Material is wrong

#cell_vacuum    = openmc.Cell(name='vacuum', fill=vacuum, region=vacuumRegion&slice)

root_univ = openmc.Universe(cells=[cell_fuel,cell_zirc,cell_blanket,cell_steelClad,cell_carbSteel,cell_air,cell_blastShield])
geometry = openmc.Geometry(root_univ)

mats.export_to_xml()
geometry.export_to_xml()

"""
plot = openmc.Plot()
plot.type = "voxel"
plot.filename = 'voxel_plot'
plot.width = (hTot, hTot, hTot)
plot.pixels = (400, 400, 200)
"""

plot = openmc.Plot()
plot.filename = 'xy_slice'
plot.colors = {air:"white",D2O:"lightblue",SS304:"lightgrey",S4340:"grey",zirconium:"black"}
plot.width = (hTot, hTot)     # Width of plot in cm
plot.pixels = (4000, 4000)        # Image resolution
plot.basis = 'xz'               # Slice plane
plot.color_by = 'material'     # Color by material
plots = openmc.Plots([plot])
plots.export_to_xml()
openmc.plot_geometry()

#root_univ.plot()
#plt.savefig("git/reaktorProject/Mads/testfig3")