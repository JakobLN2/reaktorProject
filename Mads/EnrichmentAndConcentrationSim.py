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
def inches_to_cm(inches):
    return inches * 2.54

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

zirconium = openmc.Material(2, "zirconium", temperature=900)
zirconium.add_element('Zr', 1.0, 'wo')
zirconium.add_element('Sn', 0.015, 'wo')
zirconium.set_density('g/cm3', 6.6)


# ——— Heavy water moderator (D₂O) ———  
# Density ≈1.1056 g/cm³ at 25 °C :contentReference[oaicite:0]{index=0}  #ved egentlig ikke hvad cooling water temperaturen er 

T= 600
P = 150

D2O = openmc.Material(name='heavy water')
D2O.set_density('g/cm3', 1.1056)
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
def makeFuelOld(rho,g_per_kgU,enrichment,material_id =None):
                             #of uranium to heavy water
    weight_frac = g_per_kgU/1000



    D2O_frac = 235 / (20) / weight_frac #Idk how this is calculated, check later
    

    if material_id is None:
        fuel = openmc.Material(name='fuel solution')
    else:
        fuel = openmc.Material(material_id=material_id, name='fuel solution')
    fuel.set_density('g/cm3', rho)            # assume same density as pure D₂O
    # heavy‐water fraction ≈1–0.0099=0.9901 by mass  
    fuel.add_nuclide('H2', D2O_frac * 2, 'ao')
    fuel.add_nuclide('O16', D2O_frac + 6, 'ao')
    # uranium                       
    #fuel.add_nuclide('U235', 0.93, 'ao')       #Primitive enrichment scheme
    fuel.add_element("U",percent=1,percent_type="ao",enrichment=enrichment)
    #fuel.add_nuclide('U238', 0.07, 'ao')
    #Sulpher
    fuel.add_element('S', 1, 'ao')
    print(fuel.get_elements())
    return fuel

def makeFuel(rho,g_per_kgU,enrichment,fuel_id =None):
        
        fuel = openmc.Material(name='fuel solution',material_id =fuel_id)
        fuel.set_density('g/cm3', rho)
        weightProcU = g_per_kgU/1000
        #Not sure what will happen, when we get mass of elements that aren't in fuel
        mU = openmc.data.atomic_weight("U") #Should take into account the enrichment, but this sould be an okay approximation
        mD2O = 2*openmc.atomic_mass("H2") + openmc.atomic_mass("O16")
        mUSO4 = openmc.data.atomic_weight("S") + 4*openmc.atomic_mass("O16")+mU
        print("fuelcheck", mUSO4)
        R = (mU/weightProcU -mUSO4 )/mD2O #Ratio between number of D20 and USO4
        D2O_num_frac = R/(1+R)*100 #Mass fraction of D2O
        USO4_num_frac = 1/(1+R)*100 #Mass fraction of USO4
        
        fuel.add_nuclide("U238",percent=(1-enrichment/100)*USO4_num_frac/6,percent_type="ao")
        fuel.add_nuclide("U235",percent=(enrichment/100)*USO4_num_frac/6,percent_type="ao")
        fuel.add_element('S', percent=USO4_num_frac/6,percent_type="ao")
        fuel.add_element("O",percent=4*USO4_num_frac/6,percent_type="ao")

        fuel.add_nuclide("H2",percent=2*D2O_num_frac/3.,percent_type="ao")
        fuel.add_nuclide("O16",percent=D2O_num_frac/3,percent_type="ao")
        
        return fuel
        #USO4_g_pr_kg = g_per_kgU* fuel.get_mass("U") /(fuel.get_mass("U") + fuel.get_mass("S") + 4*fuel.get_mass("O"))
        
        

        
         



def modifyFuel(fuel, rho, g_per_kgU, enrichment):
    """
    Function to modify the fuel material with the given parameters
    """
    
    elements = fuel.get_elements()
    print("DEBUG!")
    
    for i in range(len(elements)):
        fuel.remove_element(elements[i])
    fuel.set_density('g/cm3', rho)
    weightProcU = g_per_kgU/1000
    #Not sure what will happen, when we get mass of elements that aren't in fuel
    mU = openmc.data.atomic_weight("U") #Should take into account the enrichment, but this sould be an okay approximation
    mD2O = 2*openmc.atomic_mass("H2") + openmc.atomic_mass("O16")
    mUSO4 = openmc.data.atomic_weight("S") + 4*openmc.atomic_mass("O16")+mU
    print("fuelcheck", mUSO4)
    R = (mU/weightProcU -mUSO4 )/mD2O #Ratio between number of D20 and USO4
    if np.isinf(R):
        D2O_num_frac = 100.
        USO4_num_frac = 0.
    else:
        D2O_num_frac = R/(1+R)*100 #Mass fraction of D2O
        USO4_num_frac = 1/(1+R)*100 #Mass fraction of USO4
    #fuel.add_element("U",percent=USO4_num_frac/6,percent_type="ao",enrichment=enrichment)
    fuel.add_nuclide("U238",percent=(1-enrichment/100)*USO4_num_frac/6,percent_type="ao")
    fuel.add_nuclide("U235",percent=(enrichment/100)*USO4_num_frac/6,percent_type="ao")
    fuel.add_element('S', percent=USO4_num_frac/6,percent_type="ao")

    fuel.add_element("O",percent=4*USO4_num_frac/6,percent_type="ao")
    fuel.add_nuclide("H2",percent=2*D2O_num_frac/3.,percent_type="ao")
    fuel.add_nuclide("O16",percent=D2O_num_frac/3,percent_type="ao")



T= 300
P = 300
"""
g_per_kg = 10.4                             #of uranium to heavy water
weight_frac = 10.4/1000
D2O_frac = 235 / (20) / weight_frac
print(D2O_frac)

fuel = openmc.Material(name='fuel solution')
fuel.set_density('g/cm3', rho(P,T))            # assume same density as pure D₂O
# heavy‐water fraction ≈1–0.0099=0.9901 by mass  
fuel.add_nuclide('H2', D2O_frac * 2, 'ao')
fuel.add_nuclide('O16', D2O_frac + 6, 'ao')
# uranium                       
fuel.add_nuclide('U235', 0.93, 'ao')       #Bruger highly enriched uranium
fuel.add_nuclide('U238', 0.07, 'ao')
#Sulpher
fuel.add_element('S', 1, 'ao')
"""

fuel = makeFuel(rho(P,T),10.4,93) #g_per_kgU, enrichment


# fuel.add_elements_from_formula()





mats = openmc.Materials([SS304, D2O, zirconium, fuel])

print(mats[3].nuclides[2])


clad_inner = openmc.Sphere(r = inches_to_cm(16))
clad_outer = openmc.Sphere(r = inches_to_cm(16 + 5/16))    
PV_inner = openmc.Sphere(r = inches_to_cm(30))              
PV_outer = openmc.Sphere(r = inches_to_cm(30 + 4.4), boundary_type = 'vacuum')      #base er apparently at alle lag er transmissive som boundary men kan ikke have transmissive boundary som yderste, skal være periodic, reflective eller vacuum    


fuel_region = -clad_inner                                  #~10 g U per kg D2O, UO_2SO_4 in heavy water?
core_vessel = +clad_inner & -clad_outer                    #zirconium alloy
moderator_region = +clad_outer & -PV_inner                 #Heavy water
blast_shield = +PV_inner & -PV_outer                       #Stainless steal 304



cell_fuel      = openmc.Cell(name='fuel', fill=fuel, region=fuel_region)
cell_clad      = openmc.Cell(name='clad (Zr)', fill=zirconium, region=core_vessel)
cell_moderator = openmc.Cell(name='moderator', fill=D2O, region=moderator_region)
cell_shield    = openmc.Cell(name='blast shield', fill=SS304, region=blast_shield)

# root universe & geometry
root_univ = openmc.Universe(cells=[cell_fuel, cell_clad, cell_moderator, cell_shield])  #laver en celle af de celler vi har defineret som man så kan gentage, kinda unødvendigt no?
geometry = openmc.Geometry(root_univ)                                                 #laver bare en klon af geometrien


def changeUraniumInMat(baseMaterial,concentration,Enrichment):
    """
    Function to change the USO4 concentration and enrichment in material

    Arguments:
        baseMaterial: the material to change the enrichment/concentration in

        Enrichment: then enrichment percent.

        concentration: The concentration in g pr L urainum.
    """
    #Clone the base material, so we don't do anything stupid
    ourMat = baseMaterial
    #Add our new enrichment
    ourMat.remove_element("U")
    ourMat.add_nuclide('U235', concentration/100., 'ao')       #Bruger highly enriched uranium
    ourMat.add_nuclide('U238', 1-concentration/100., 'ao')
    #ourMat.add_element("U",percent=percentFuel,percent_type="ao",enrichment=enrichment)
    return ourMat

def changeEnrichmentOfMat(enrichment,baseMaterial):
    """
    Function to change the enrichment of uranium in the given material

    Arguments:
        baseMaterial: the material to change the enrichment in

        enrichment: The enrichment in weight percent.
    """
    #Clone the base material, so we don't do anything stupid
    ourMat = baseMaterial
    #Add our new enrichment
    ourMat.remove_element("U")
    ourMat.add_nuclide('U235', enrichment/100., 'ao')       #Bruger highly enriched uranium
    ourMat.add_nuclide('U238', 1-enrichment/100., 'ao')
    #ourMat.add_element("U",percent=percentFuel,percent_type="ao",enrichment=enrichment)
    return ourMat


def kEffAtEnrichment(model,enrichments,fuelId):
    """
    Arguments:
        Enrichments:   Liste af enrichments man vil finde k eff af

        model:    Base model vi bruger; OBS: funktionen ændre på den
    Returns:
        keffs for all the given **enrichments**
    """
    kEffs = np.zeros(len(enrichments))
    kEffstds = np.zeros(len(enrichments))


    for i,enrich in enumerate(enrichments):
        print("Current step: ", i)

        model.materials[fuelId] = changeEnrichmentOfMat(enrich,model.materials[fuelId])
        #print(model.materials[fuelId])
        
        sp_path = model.run()                        #at den returnerer en path has to be retarded
        
        with openmc.StatePoint(sp_path) as sp:
            k_eff = sp.keff
            
            # k_avg = sp.keff[0].mean

            keff_var = sp.keff

            # Extract nominal value and uncertainty
            kEffs[i] = keff_var.nominal_value   # or keff_var.n
            kEffstds[i]  = keff_var.std_dev          # or keff_var.s
    return kEffs, kEffstds





def kEffAtConcentrationAndEnrichments(model,concentrations,enrichments,rho,fuelId):
    """
    Arguments:
        concentrations:   Liste af enrichments
    """
    kEffs = np.zeros((len(concentrations),len(enrichments)))
    kEffstds = np.zeros((len(concentrations),len(enrichments)))


    for i,conc in enumerate(concentrations):
        for j,enrich in enumerate(enrichments):
            print("Current step: ", i)
            modifyFuel(model.materials[fuelId],rho,conc,enrich) #This just needs to be garbage
            print(model.materials[fuelId])
            #model.materials[fuelId] = modifyFuel(model.materials[fuelId],rho,conc,enrich)
            #print(model.materials[fuelId])
            
            sp_path = model.run()                        #at den returnerer en
            with openmc.StatePoint(sp_path) as sp:
                k_eff = sp.keff
                
                # k_avg = sp.keff[0].mean

                keff_var = sp.keff

                # Extract nominal value and uncertainty
                kEffs[i,j] = keff_var.nominal_value   # or keff_var.n
                kEffstds[i,j]  = keff_var.std_dev          # or keff_var.s
    return kEffs, kEffstds



geometry.export_to_xml()

root_univ.plot()

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

settings.export_to_xml()

# tally k-effective
tally = openmc.Tally(name='keff')
tally.estimator = 'analog'
tally.filters = []
tally.scores = ['kappa-fission']
# model = openmc.model.Model(geometry, settings, materials=[SS304, zirconium, D2O, fuel], 
#                            tallies=openmc.Tallies([tally]))
model = openmc.model.Model(
    geometry=geometry,
    settings=settings,
    materials=mats,           # your openmc.Materials object
    tallies=openmc.Tallies([tally])
)

enrichments = np.linspace(1,100,50) #Enrichment can't go above 99.2 for some reason. Idk, have to investigate
concentrations = np.linspace(1,100,50) #g_pr_kg




kEffs, stdKEffs = kEffAtConcentrationAndEnrichments(model,concentrations,enrichments,rho(120, 280),3)
fig,ax = plt.subplots(1,1)
ax.contourf(enrichments,concentrations,kEffs)
ax.set_xlabel("Enrichment")
ax.set_ylabel("Concentration")
ax.set_title("K eff vs Enrichment and Concentration")
ax.set_xlim(0,100)
ax.set_ylim(0,100)
ax.set_aspect('equal')

np.savez("git/reaktorProject/Mads/Keffs",enrichments = enrichments,concentrations = concentrations, keffs = kEffs,stdKeffs = stdKEffs) #Save the data, so we can make a nice looking plot without having to run the entire simulation again.

fig.savefig("git/reaktorProject/Mads/keffvsEnrichmentDebug")


"""
# run and capture summary   
sp_path = model.run()                           #at den returnerer en path has to be retarded
sp = openmc.StatePoint(sp_path)
# k_avg = sp.keff[0].mean

keff_var = sp.keff

# Extract nominal value and uncertainty
mean_keff = keff_var.nominal_value   # or keff_var.n
std_keff  = keff_var.std_dev          # or keff_var.s

print(f"Estimated k-effective = {mean_keff:.5f} ± {std_keff:.5f}")

k_gen = sp.k_generation         #finder k_effective for hver generation
batches = np.arange(1, len(k_gen) + 1)


plt.plot(batches, k_gen, marker='o', linestyle='-')
plt.xlabel('Batch number')
plt.ylabel('k-effective per generation')
plt.title('Convergence of k-effective (raw)')
plt.grid(True)


cum_sum = np.cumsum(k_gen)
# running mean
k_cumavg = cum_sum / batches

plt.figure()
plt.plot(batches, k_cumavg, marker='o', linestyle='-')
plt.xlabel('Batch number')
plt.ylabel('Cumulative average k-effective')
plt.title('Convergence of k-effective (cumulative average)')
plt.grid(True)

plt.savefig("git/reaktorProject/Mads/testfig1")
"""
plt.show()
