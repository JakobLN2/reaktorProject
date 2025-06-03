import openmc
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

"""
Fandt ud af at du bare kan gøre det her, så er openmc configered til at tage de korrekte cross section files
"""


import scipy.interpolate as scint


matplotlib.use('webagg')

OutName = "ShortPipes"


def inches_to_cm(inches):
    return inches * 2.54

"""
Define the geometry of the reactor
"""


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

        
        ball = -openmc.Sphere(r=self.rBall)
        topCylinder = -openmc.Cylinder(r=self.rOut, dz=self.hTot/2) &-openmc.ZPlane(z0=self.hTop) &+openmc.ZPlane(z0=0)
        bottomCylinder = -openmc.Cylinder(r=self.rIn, z0=-self.hTot/2, dz=self.hTot/2) &-openmc.ZPlane(z0=0)&+openmc.ZPlane(z0=-self.hBottom)
        
        #print(cone1h)
        cone1 = -openmc.Cone(r2=np.power(self.RCone1/self.HCone1,2), z0=-self.z0Cone1, dz=-1) &+openmc.ZPlane(z0=-self.z0Cone1)&-openmc.ZPlane(z0=-self.z0Cone1/2)
        print(self.z0Cone1,self.hCone1,self.rCone1,self.RCone1)
  
        
        #cone2h = cone2R/np.tan(self.theta2)
        #print(cone2h,cone2R)
        #Why are the boundaries different? Idk
        cone2 = -openmc.Cone(r2=np.power(self.RCone2/self.HCone2,2), z0=-self.z0Cone2, dz=-1)&+openmc.ZPlane(z0=-self.z0Cone2)&-openmc.ZPlane(z0=-self.z0Cone1/2)
            
        return  ball |cone2| cone1 | topCylinder | bottomCylinder



        

rBall = inches_to_cm(16)
coreInlet = inches_to_cm(3.5/2)#Think they were referencing diameter
coreOutlet = inches_to_cm(4/2)
#coreBallPlateThickness = inches_to_cm(5/16)
#coreConePlateThickness = inches_to_cm(3/8) #Gonna need to assume this is the same as the ball plate thickness

hTot = inches_to_cm(11*12+6+3/16)
hTop = inches_to_cm(6*12+10+3/8)
#startPlate = inches_to_cm(7*12+9+1/4) - inches_to_cm(6*12+10+3/8)
rCone1 = inches_to_cm(4+15/16)
zircThickness = inches_to_cm(5/16)


cutPlane1 = openmc.YPlane(y0=0,boundary_type='periodic') #This is the plane that cuts the geometry in half, so we can see the inside

angle = np.pi/4

cutPlane2 = openmc.Plane(a=np.tan(angle),B=-1,D=0,boundary_type="periodic")
cutPlane1.periodic_surface = cutPlane2
cutPlane2.periodic_surface = cutPlane1
#slice = -cutPlane1 & -cutPlane2

fuelGeometry = InnerGeometry(rBall,rCone1,coreOutlet,coreInlet,hTot,hTop,Inches=False)
zircGeometry = InnerGeometry(rBall+zircThickness,rCone1+zircThickness/np.cos(np.pi/4),coreOutlet+zircThickness,coreInlet+zircThickness,hTot,hTop,Inches=False)
fuelRegion = fuelGeometry.totRegion(outside=False)
zircRegion = zircGeometry.totRegion(outside=False)&(~fuelRegion)

#fuelRegion = -openmc.Sphere(r=rBall)
#zircGeometry = -openmc.Sphere(r=rBall+zircThickness)
#zircRegion = zircGeometry&(~fuelRegion) 


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
blastShieldGeometry = -openmc.Sphere(r=rblastShield,boundary_type='vacuum')#boundary_type='vacuum'
blastShieldRegion = blastShieldGeometry & (~airRegion) & (~blanketRegion) & (~steelCladRegion) & (~carbSteelRegion) & (~fuelRegion) & (~zircRegion)

vacuumBound = -openmc.ZCylinder(r=rblastShield+0.1,boundary_type='vacuum')&+openmc.ZPlane(z0=-fuelGeometry.hBottom-0.1,boundary_type='vacuum')&-openmc.ZPlane(z0=fuelGeometry.hTop+0.1,boundary_type='vacuum') #Big cylinder of vacuum around the whole thing
vacuumRegion = vacuumBound & (~airRegion) & (~blanketRegion) & (~steelCladRegion) & (~carbSteelRegion) & (~fuelRegion) & (~zircRegion)#Not sure if this is needed or not.
#slice = blastShieldGeometry
#slice = None




#Small table - Below 275 C
T = np.loadtxt(r"git/reaktorProject/russianWater.txt",max_rows=1)
P = np.loadtxt(r"git/reaktorProject/russianWater.txt",usecols=[0],skiprows=1)
data = np.loadtxt(r"git/reaktorProject/russianWater.txt",skiprows=1,usecols=range(1,len(T)+1))
f_1 = scint.RegularGridInterpolator((P,T), data, method="linear")

#Russian table - above or equal 275 C
T = np.loadtxt(r"git/reaktorProject/russianWater.txt",max_rows=1)
P = np.loadtxt(r"git/reaktorProject/russianWater.txt",usecols=[0],skiprows=1)*0.9807 # Convert to bar
data = np.loadtxt(r"git/reaktorProject/russianWater.txt",skiprows=1,usecols=range(1,len(T)+1))
f_2 = scint.RegularGridInterpolator((P,T), data, method="linear")
rho_2 = lambda P,T: 1/f_2((P,T))

rho = lambda P,T: f_1((P,T))/1000 if T < 275 else rho_2(P,T)




openmc.config["cross_sections"] = "jeff-3.3-hdf5/cross_sections.xml"
def inches_to_cm(inches):
    return inches * 2.54

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


zirconium = openmc.Material(name="zirconium",temperature=900)
zirconium.set_density('g/cm3', 6.6)
zirconium.add_element('Zr', 1.0, 'wo')
zirconium.add_element('Sn', 0.015, 'wo')

air = openmc.Material(name='air')
air.set_density('g/cm3', 1.205e-3) #density at room temp
air.add_element('N', 0.7808, 'wo')
air.add_element('O', 0.2095, 'wo')
#air.add_element('Ar', 0.0093, 'wo')



vacuum = openmc.Material(name='vacuum')
vacuum.set_density('g/cm3', 1e-50) #density at room temp
vacuum.add_element('N', 0.7808, 'wo')
vacuum.add_element('O', 0.2095, 'wo')
vacuum.add_element('Ar', 0.0093, 'wo')

# ——— Heavy water moderator (D₂O) ———  
# Density ≈1.1056 g/cm³ at 25 °C :contentReference[oaicite:0]{index=0}  #ved egentlig ikke hvad cooling water temperaturen er 

T= 280
P = 120

D2O = openmc.Material(name='heavy water')
D2O.set_density('g/cm3', rho(P,T))
# use nuclide H2 for deuterium and O16 for oxygen  
D2O.add_nuclide('H2', 2.0, 'ao')
D2O.add_nuclide('O16', 1.0, 'ao')


"""
tror du kan tilføje vandet mere clean via .add_elements_from_formula, ikke sikker
"""

"""

"""


"""
Assume density of water isn't impacted by uranium dissolved in it
"""

def makeFuel(rho,g_per_kgU,enrichment,fuel_id =None):
        
        fuel = openmc.Material(name='fuel solution',material_id =fuel_id)
        fuel.set_density('g/cm3', rho)
        U_weightFrac = g_per_kgU/1000
        #waterMass = 20
        #Umass = 235
        
        
        
        
        #Not sure what will happen, when we get mass of elements that aren't in fuel
        mU = openmc.data.atomic_weight("U") #Should take into account the enrichment, but this sould be an okay approximation
            
        mD2 = openmc.atomic_mass("H2")
        mO = openmc.atomic_mass("O16")
        mS = openmc.data.atomic_weight("S")
        mD2O = 2*mD2 + mO
        mUO2SO4 = mS + 6*mO+mU
            
            

            
        D2O_weightFrac = 1-mUO2SO4/mU*U_weightFrac
        O_weightFrac = 6*mO/mU*U_weightFrac+mO/mD2O*D2O_weightFrac #Oxygen from 
        H2_weightFrac = mD2/mD2O*D2O_weightFrac
        S_weightFrac = mS/mU*U_weightFrac


            
        fuel.add_nuclide("U238",percent=(1-enrichment/100)*U_weightFrac*100,percent_type="wo")
        fuel.add_nuclide("U235",percent=(enrichment/100)*U_weightFrac*100,percent_type="wo")
        fuel.add_element('S', percent=S_weightFrac*100,percent_type="wo")
        fuel.add_nuclide("O16",percent=O_weightFrac*100,percent_type="wo")

        fuel.add_nuclide("H2",percent=H2_weightFrac*100,percent_type="wo")
        return fuel
        #USO4_g_pr_kg =        

        
         


#fuel = makeFuel(rho(P,T),10.5,93) #g_per_kgU, enrichment

g_per_kg = 10.4                             #of uranium to heavy water
weight_frac = 10.4/1000
D2O_frac = 235 / (20) / weight_frac
# print(D2O_frac)

fuel = openmc.Material(name='fuel solution')
fuel.set_density('g/cm3', 0.8443)            # assume same density as pure D₂O
fuel.add_nuclide('H2', D2O_frac * 2, 'ao')
fuel.add_nuclide('O16', D2O_frac + 6, 'ao') #TODO This still feels scuffed
# uranium                       
fuel.add_nuclide('U235', 0.93, 'ao') 
fuel.add_nuclide('U238', 0.07, 'ao')
fuel.add_element('S', 1, 'ao')

print(fuel)
print(rho(P,T))

# fuel.add_elements_from_formula()





#mats = openmc.Materials([SS304, D2O, zirconium, fuel])


mats = openmc.Materials([fuel,D2O,SS304,zirconium,S4340,air,vacuum])
#mats.export_to_xml()
# Export the geometry to XML


#TODO Why the fuck doesn't this work :(

print("Here")
cell_fuel      = openmc.Cell(name='fuel', fill=fuel, region=fuelRegion                  &blastShieldGeometry)
cell_zirc      = openmc.Cell(name='zirc', fill=zirconium, region=zircRegion             &blastShieldGeometry)
cell_blanket   = openmc.Cell(name='blanket', fill=D2O, region=blanketRegion             &blastShieldGeometry) 
cell_steelClad = openmc.Cell(name='steelClad', fill=SS304, region=steelCladRegion       &blastShieldGeometry) 
cell_carbSteel = openmc.Cell(name='carbSteel', fill=S4340, region=carbSteelRegion       &blastShieldGeometry)
cell_air       = openmc.Cell(name='air', fill=air, region=airRegion                     &blastShieldGeometry)
cell_blastShield = openmc.Cell(name='blastShield', fill=S4340, region=blastShieldRegion &blastShieldGeometry)
#cell_vacuum    = openmc.Cell(name='vacuum', fill=vacuum, region=vacuumRegion           )
root_univ = openmc.Universe(cells=[cell_fuel,cell_zirc,cell_blanket,cell_steelClad,cell_carbSteel,cell_air,cell_blastShield])#,cell_carbSteel,cell_air,cell_blastShield,cell_vacuum

# root universe & geometry
#root_univ = openmc.Universe(cells=[cell_fuel,cell_zirc,cell_blanket,cell_steelClad,cell_carbSteel,cell_air,cell_blastShield])
geometry = openmc.Geometry(root_univ)                                                 #laver bare en klon af geometrien




geometry.export_to_xml()



plt.savefig("git/reaktorProject/Mads/testfig2")


settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.batches = 20
settings.inactive = 5
settings.particles = 20000
settings.temperature['default'] = 600           #close enough

# point‐source at center, isotropic, Watt fission spectrum
# source = openmc.Source()
# source.space = openmc.stats.Point((0, 0, 0))
# source.angle = openmc.stats.Isotropic()
# source.energy = openmc.stats.Watt(a=0.988e6, b=2.249e-6)  # typical U-235
# settings.source = source

#settings.export_to_xml()

# tally k-effective
tally = openmc.Tally(name='keff')
tally.estimator = 'analog'
tally.filters = []
tally.scores = ["flux"] #'kappa-fission'
# model = openmc.model.Model(geometry, settings, materials=[SS304, zirconium, D2O, fuel], 
#                            tallies=openmc.Tallies([tally]))
model = openmc.model.Model(
    geometry=geometry,
    settings=settings,
    materials=mats,           # your openmc.Materials object
    tallies=openmc.Tallies([tally])
)

sp_path = model.run(threads=8)                           #at den returnerer en path has to be retarded
sp = openmc.StatePoint(sp_path)

k_gen = sp.k_generation         #finder k_effective for hver generation

batches = np.arange(1, len(k_gen) + 1)


np.savez(OutName + "_Data",batches = batches, k_gen = k_gen)

fig,ax = plt.subplots(1,1)

sp.tallies[1]

ax.plot(batches, k_gen, marker='o', linestyle='-')
ax.set_xlabel('Batch number')
ax.set_ylabel('k-effective per generation')
ax.set_title('Convergence of k-effective (raw)')
ax.grid(True)
#plt.savefig("git/reaktorProject/Mads/convergence_advanced_geometry.png")

plot = openmc.Plot()
plot.filename = OutName
plot.colors = {air:"white",D2O:"lightblue",SS304:"lightgrey",S4340:"grey",fuel:"lightgreen",zirconium:"black"}
plot.width = (hTot, hTot)     # Width of plot in cm
plot.pixels = (4000, 4000)        # Image resolution
plot.basis = 'xz'               # Slice plane
plot.color_by = 'material'     # Color by material
plots = openmc.Plots([plot])
plots.export_to_xml()
openmc.plot_geometry()



#np.savez("git/reaktorProject/Mads/Keffs",enrichments = enrichments,concentrations = concentrations, keffs = kEffs,stdKeffs = stdKEffs) #Save the data, so we can make a nice looking plot without having to run the entire simulation again.

#before: 
"""
 k-effective (Collision)     = 1.09390 +/- 0.00257
 k-effective (Track-length)  = 1.09451 +/- 0.00248
 k-effective (Absorption)    = 1.09565 +/- 0.00191
 Combined k-effective        = 1.09877 +/- 0.00176
 Leakage Fraction            = 0.07126 +/- 0.00047
"""
"""
Total time for initialization     = 6.7664e+00 seconds
   Reading cross sections          = 6.7373e+00 seconds
 Total time in simulation          = 8.8938e+00 seconds
   Time in transport only          = 8.7828e+00 seconds
   Time in inactive batches        = 1.7757e+00 seconds
   Time in active batches          = 7.1181e+00 seconds
   Time synchronizing fission bank = 3.0948e-02 seconds
     Sampling source sites         = 2.6356e-02 seconds
     SEND/RECV source sites        = 3.3739e-03 seconds
   Time accumulating tallies       = 1.4165e-02 seconds
   Time writing statepoints        = 5.2321e-02 seconds
 Total time for finalization       = 1.5370e-04 seconds
 Total time elapsed                = 1.5683e+01 seconds
 Calculation Rate (inactive)       = 56316.9 particles/second
 Calculation Rate (active)         = 42146.1 particles/second
"""


plt.show()
