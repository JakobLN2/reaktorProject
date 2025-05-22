import openmc
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.interpolate as scint
plt.rc("axes", labelsize=30, titlesize=32)   # skriftstørrelse af xlabel, ylabel og title
plt.rc("xtick", labelsize=26, top=True, direction="in")  # skriftstørrelse af ticks, vis også ticks øverst og vend ticks indad
plt.rc("ytick", labelsize=26, right=True, direction="in") # samme som ovenstående
plt.rc("legend", fontsize=26) # skriftstørrelse af figurers legends
plt.rcParams["font.size"] = "20"
plt.rcParams["figure.figsize"] = (16,8)

savepath = r'/home/candifloos/Reaktorfysik/Figurer//'

# matplotlib.use('webagg')
def inches_to_cm(inches):
    return inches * 2.54

#Small table - Below 275 C
path = r'/home/candifloos/Reaktorfysik/reaktorProject/Jakob//'
T = np.loadtxt(path + 'badWater.txt',max_rows=1)
P = np.loadtxt(path + 'badWater.txt',usecols=[0],skiprows=1)
data = np.loadtxt(path + 'badWater.txt',skiprows=1,usecols=range(1,len(T)+1))
f_1 = scint.RegularGridInterpolator((P,T), data, method="linear")

#Russian table - above or equal 275 C
T = np.loadtxt(path + 'russianWater.txt',max_rows=1)
P = np.loadtxt(path + 'russianWater.txt',usecols=[0],skiprows=1)
data = np.loadtxt(path + 'russianWater.txt',skiprows=1,usecols=range(1,len(T)+1))
f_2 = scint.RegularGridInterpolator((P,T), data, method="linear")
rho_2 = lambda P,T: 1/f_2((P,T))

rho = lambda P,T: f_1((P,T))/1000 if T < 275 else rho_2(P,T)


SS304 = openmc.Material(name='SS304')
SS304.set_density('g/cc', 8.03)
SS304.add_element('Si', 0.0060, 'wo')
SS304.add_element('Cr', 0.1900, 'wo')
SS304.add_element('Mn', 0.0200, 'wo')
SS304.add_element('Fe', 0.6840, 'wo')
SS304.add_element('Ni', 0.1000, 'wo')

zirconium = openmc.Material(2, "zirconium", temperature=900)
zirconium.add_element('Zr', 0.985, 'wo')
zirconium.add_element('Sn', 0.015, 'wo')
zirconium.set_density('g/cm3', 6.6)


# ——— Heavy water moderator (D₂O) ———  
# Density ≈1.1056 g/cm³ at 25 °C :contentReference[oaicite:0]{index=0}  #actually at 300C
D2O = openmc.Material(name='heavy water')
D2O.set_density('g/cm3', rho(120, 280))
# use nuclide H2 for deuterium and O16 for oxygen  
D2O.add_nuclide('H2', 2.0, 'ao')
D2O.add_nuclide('O16', 1.0, 'ao')


g_per_kg = 10.4                             #of uranium to heavy water
weight_frac = 10.4/1000
D2O_frac = 235 / (20) / weight_frac
print(D2O_frac)

fuel = openmc.Material(name='fuel solution')
fuel.set_density('g/cm3', 0.8443 )#+ weight_frac)            # assume same density as pure D₂O
fuel.add_nuclide('H2', D2O_frac * 2, 'ao')
fuel.add_nuclide('O16', D2O_frac + 6, 'ao')
# uranium                       
fuel.add_nuclide('U235', 0.93, 'ao') 
fuel.add_nuclide('U238', 0.07, 'ao')
fuel.add_element('S', 1, 'ao')


clad_inner = openmc.Sphere(r = inches_to_cm(16))
clad_outer = openmc.Sphere(r = inches_to_cm(16 + 5/16))    
PV_inner = openmc.Sphere(r = inches_to_cm(30))              
PV_outer = openmc.Sphere(r = inches_to_cm(30 + 4.4), boundary_type = 'vacuum')      #base er apparently at alle lag er transmissive som boundary men kan ikke have transmissive boundary som yderste, skal være periodic, reflective eller vacuum    
                                                                                    #Hvis man bruger reflective ligger de for some reason meget tættere?


fuel_regions = [0]*4
fuel_cells = [0] * 4
split = openmc.YPlane(y0=0, boundary_type="transmission")
split_b1 = openmc.YPlane(y0=-inches_to_cm(8), boundary_type="transmission")
split_b2 = openmc.YPlane(y0=-inches_to_cm(16), boundary_type="transmission")
split_t1 = openmc.YPlane(y0=inches_to_cm(8), boundary_type="transmission")
split_t2 = openmc.YPlane(y0=inches_to_cm(16), boundary_type="transmission")



N = 1#128 #number of subdivisions in the core (temperature gradient)

b = -inches_to_cm(16); t = inches_to_cm(16)
T_grad = lambda y: (y - b) * (300 - 256)/(t - b) + 256

ys = np.linspace(b,t,N+1, endpoint=True)
splits = [openmc.YPlane(y,boundary_type="transmission") for y in ys]
fuel_regions = [0] * N; fuel_cells = [0] * N; fuels = [0] * N
for i in range(N):
    t_i = T_grad(0.5*(t - b)/N + ys[i])
    fuels[i] = fuel.clone()
    fuels[i].set_density('g/cm3', rho(100,t_i) + weight_frac)
    fuel_regions[i] = -clad_inner & +splits[i] & -splits[i + 1]
    fuel_cells[i] = openmc.Cell(name=f'fuel{i}', fill=fuels[i], region=fuel_regions[i])
    
mats = openmc.Materials([SS304, D2O, zirconium, *fuels])
mats.cross_sections = r'/home/candifloos/Reaktorfysik/Python/RFP/Data/jeff-3.3-hdf5/cross_sections.xml'


# fuel_region = -clad_inner                                  #~10 g U per kg D2O, UO_2SO_4 in heavy water?
# cell_fuel      = openmc.Cell(name='fuel', fill=fuel, region=fuel_region)

core_vessel = +clad_inner & -clad_outer                    #zirconium alloy
moderator_region = +clad_outer & -PV_inner                 #Heavy water
pressure_vessel = +PV_inner & -PV_outer                       #Stainless steal 304


cell_clad      = openmc.Cell(name='clad (Zr)', fill=zirconium, region=core_vessel)
cell_moderator = openmc.Cell(name='moderator', fill=D2O, region=moderator_region)
cell_shield    = openmc.Cell(name='blast shield', fill=SS304, region=pressure_vessel)

# root universe & geometry
# root_univ = openmc.Universe(cells=[cell_fuel, cell_clad, cell_moderator, cell_shield])  #laver en celle af de celler vi har defineret som man så kan gentage, kinda unødvendigt no?
root_univ = openmc.Universe(cells=[*fuel_cells, cell_clad, cell_moderator, cell_shield]) 
geometry = openmc.Geometry(root_univ)                                                   #laver bare en klon af geometrien
geometry.export_to_xml()

root_univ.plot()
plt.savefig(savepath + 'Spherical reactor.png', bbox_inches='tight')

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.batches = 50
settings.inactive = 4
settings.particles = 20000
settings.temperature['default'] = 600           #close enough

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


plt.figure()
plt.plot(batches, k_gen, marker='o', linestyle='-')
plt.xlabel('Batch number')
plt.ylabel('k-effective per generation')
plt.title('Convergence of k-effective (raw)')
plt.grid(True)


cum_sum = np.cumsum(k_gen)
# running mean
k_cumavg = cum_sum / batches
plt.savefig(savepath + 'keff over time.png', bbox_inches='tight')


plt.figure()
plt.plot(batches, k_cumavg, marker='o', linestyle='-')
plt.xlabel('Batch number')
plt.ylabel('Cumulative average k-effective')
plt.title('Convergence of k-effective (cumulative average)')
plt.grid(True)

plt.savefig(savepath + 'keff convergence.png', bbox_inches='tight')
sp.close()










fuel_regions = [0]*4
fuel_cells = [0] * 4
split = openmc.YPlane(y0=0, boundary_type="transmission")
split_b1 = openmc.YPlane(y0=-inches_to_cm(8), boundary_type="transmission")
split_b2 = openmc.YPlane(y0=-inches_to_cm(16), boundary_type="transmission")
split_t1 = openmc.YPlane(y0=inches_to_cm(8), boundary_type="transmission")
split_t2 = openmc.YPlane(y0=inches_to_cm(16), boundary_type="transmission")



N = 128 #number of subdivisions in the core (temperature gradient)

b = -inches_to_cm(16); t = inches_to_cm(16)
T_grad = lambda y: (y - b) * (300 - 256)/(t - b) + 256

ys = np.linspace(b,t,N+1, endpoint=True)
splits = [openmc.YPlane(y,boundary_type="transmission") for y in ys]
fuel_regions = [0] * N; fuel_cells = [0] * N; fuels = [0] * N
for i in range(N):
    t_i = T_grad(0.5*(t - b)/N + ys[i])
    fuels[i] = fuel.clone()
    fuels[i].set_density('g/cm3', rho(100,t_i) + weight_frac)
    fuel_regions[i] = -clad_inner & +splits[i] & -splits[i + 1]
    fuel_cells[i] = openmc.Cell(name=f'fuel{i}', fill=fuels[i], region=fuel_regions[i])
    
mats = openmc.Materials([SS304, D2O, zirconium, *fuels])
mats.cross_sections = r'/home/candifloos/Reaktorfysik/Python/RFP/Data/jeff-3.3-hdf5/cross_sections.xml'
    



# fuel_region = -clad_inner                                  #~10 g U per kg D2O, UO_2SO_4 in heavy water?
# cell_fuel      = openmc.Cell(name='fuel', fill=fuel, region=fuel_region)

core_vessel = +clad_inner & -clad_outer                    #zirconium alloy
moderator_region = +clad_outer & -PV_inner                 #Heavy water
pressure_vessel = +PV_inner & -PV_outer                       #Stainless steal 304


cell_clad      = openmc.Cell(name='clad (Zr)', fill=zirconium, region=core_vessel)
cell_moderator = openmc.Cell(name='moderator', fill=D2O, region=moderator_region)
cell_shield    = openmc.Cell(name='blast shield', fill=SS304, region=pressure_vessel)

# root universe & geometry
# root_univ = openmc.Universe(cells=[cell_fuel, cell_clad, cell_moderator, cell_shield])  #laver en celle af de celler vi har defineret som man så kan gentage, kinda unødvendigt no?
root_univ = openmc.Universe(cells=[*fuel_cells, cell_clad, cell_moderator, cell_shield]) 
geometry = openmc.Geometry(root_univ)                                                   #laver bare en klon af geometrien
geometry.export_to_xml()

root_univ.plot()
plt.savefig(savepath + 'Spherical reactor, grad.png', bbox_inches='tight')

settings = openmc.Settings()
settings.run_mode = 'eigenvalue'
settings.batches = 50
settings.inactive = 4
settings.particles = 20000
settings.temperature['default'] = 600           #close enough

settings.export_to_xml()

# tally k-effective
tally = openmc.Tally(name='keff')
tally.estimator = 'analog'
tally.filters = []
tally.scores = ['kappa-fission']

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

keff_var_grad = sp.keff

# Extract nominal value and uncertainty
mean_keff = keff_var_grad.nominal_value   # or keff_var.n
std_keff  = keff_var_grad.std_dev          # or keff_var.s

print(f"Estimated k-effective = {mean_keff:.5f} ± {std_keff:.5f}")

k_gen_grad = sp.k_generation         #finder k_effective for hver generation
batches_grad = np.arange(1, len(k_gen_grad) + 1)


plt.figure()
plt.plot(batches, k_gen, marker='o', linestyle='-', label = 'No gradient')
plt.plot(batches, k_gen_grad, marker='o', linestyle='-', label = 'Gradient')
plt.xlabel('Batch number')
plt.ylabel('k-effective per generation')
plt.title('Convergence of k-effective (raw)')
plt.grid(True)
plt.legend()


cum_sum_grad = np.cumsum(k_gen_grad)
# running mean
k_cumavg_grad = cum_sum_grad / batches_grad
plt.savefig(savepath + 'keff over time, grad.png', bbox_inches='tight')


plt.figure()
plt.plot(batches_grad, k_cumavg, marker='o', linestyle='-', label = 'No gradient')
plt.plot(batches_grad, k_cumavg_grad, marker='o', linestyle='-', label = 'Gradient')
plt.xlabel('Batch number')
plt.ylabel('Cumulative average k-effective')
plt.title('Convergence of k-effective (cumulative average)')
plt.grid(True)
plt.legend()
plt.savefig(savepath + 'keff convergence, grad.png', bbox_inches='tight')



plt.show()