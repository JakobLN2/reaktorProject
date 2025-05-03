import openmc
import os
import numpy as np
import matplotlib.pyplot as plt


"""
Fandt ud af at du bare kan gøre det her, så er openmc configered til at tage de korrekte cross section files
"""
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
D2O = openmc.Material(name='heavy water')
D2O.set_density('g/cm3', 1.1056)
# use nuclide H2 for deuterium and O16 for oxygen  
D2O.add_nuclide('H2', 2.0, 'ao')
D2O.add_nuclide('O16', 1.0, 'ao')


g_per_kg = 10.4                             #of uranium to heavy water
weight_frac = 10.4/1000
D2O_frac = 235 / (20) / weight_frac
print(D2O_frac)

fuel = openmc.Material(name='fuel solution')
fuel.set_density('g/cm3', 1.1056)            # assume same density as pure D₂O
# heavy‐water fraction ≈1–0.0099=0.9901 by mass  
fuel.add_nuclide('H2', D2O_frac * 2, 'ao')
fuel.add_nuclide('O16', D2O_frac + 6, 'ao')
# uranium                       
fuel.add_nuclide('U235', 0.93, 'ao')       #Bruger highly enriched uranium
fuel.add_nuclide('U238', 0.07, 'ao')
fuel.add_element('S', 1, 'ao')

# fuel.add_elements_from_formula()





mats = openmc.Materials([SS304, D2O, zirconium, fuel])
#mats.cross_sections = r'/home/candifloos/Reaktorfysik/Python/RFP/Data/jeff-3.3-hdf5/cross_sections.xml'


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
geometry = openmc.Geometry(root_univ)                                                   #laver bare en klon af geometrien
geometry.export_to_xml()






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
plt.show()