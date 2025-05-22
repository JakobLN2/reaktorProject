import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as scint


matplotlib.use('webagg')
plt.rc("axes", labelsize=30, titlesize=32)   # skriftstørrelse af xlabel, ylabel og title
plt.rc("xtick", labelsize=26, top=True, direction="in")  # skriftstørrelse af ticks, vis også ticks øverst og vend ticks indad
plt.rc("ytick", labelsize=26, right=True, direction="in") # samme som ovenstående
plt.rc("legend", fontsize=26) # skriftstørrelse af figurers legends
plt.rcParams["font.size"] = "20"
plt.rcParams["figure.figsize"] = (16,8)

def inches_to_cm(inches):
    return inches * 2.54

path = r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/"
#Small table - Below 275 C
T = np.loadtxt(path + "badWater.txt",max_rows=1)
P = np.loadtxt(path + "badWater.txt",usecols=[0],skiprows=1)
data = np.loadtxt(path + "badWater.txt",skiprows=1,usecols=range(1,len(T)+1))
f_1 = scint.RegularGridInterpolator((P,T), data, method="linear")

#Russian table - above or equal 275 C
T = np.loadtxt(path + "russianWater.txt",max_rows=1)
P = np.loadtxt(path + "russianWater.txt",usecols=[0],skiprows=1)
data = np.loadtxt(path + "russianWater.txt",skiprows=1,usecols=range(1,len(T)+1))
f_2 = scint.RegularGridInterpolator((P,T), data, method="linear")

rho = lambda P,T: f_1((P,T))/1000 if T < 275 else 1/f_2((P,T))
rho_ThO2 = 10 #g/cm^3



SS304 = openmc.Material(1, name='SS304')
SS304.set_density('g/cc', 8.03)
SS304.add_element('Si', 0.0060, 'wo')
SS304.add_element('Cr', 0.1900, 'wo')
SS304.add_element('Mn', 0.0200, 'wo')
SS304.add_element('Fe', 0.6840, 'wo')
SS304.add_element('Ni', 0.1000, 'wo')

zirconium = openmc.Material(2, "zirconium", temperature=600)
zirconium.add_element('Zr', 0.985, 'wo')
zirconium.add_element('Sn', 0.015, 'wo')
zirconium.set_density('g/cm3', 6.6)


def makeGeometry(P, T):
        
    # ——— Heavy water moderator (D₂O) ———  
    th_weight = 0.6 #kg Th / kg D2O
    thO2_weight = (1 + 32/232) * th_weight #kg ThO2 / kg D2O
    rho_tot = lambda p,t: (1 + thO2_weight)/(1/rho(p,t) + thO2_weight/rho_ThO2) #assumes mixing with V = V1 + V2 (yikes)

    D2O = openmc.Material(3, name='heavy water')
    D2O.set_density('g/cm3', rho_tot(120, 280))
    D2O.volume = 4.0/3.0 * np.pi * (inches_to_cm(30)**3 - inches_to_cm(16 + 5/16)**3)
    D2O.add_nuclide('H2', 4/20, 'wo')
    D2O.add_nuclide('O16', 16/20 + 32/232*th_weight, 'wo')
    D2O.add_element('Th', th_weight, 'wo')


    # ——— Uranium fuel ———
    weight_frac = 10.4/1000
    D2O_frac = 235 / (20) / weight_frac

    fuel = openmc.Material(4, name='fuel solution')
    fuel.set_density('g/cm3', rho(P, T))
    fuel.volume = 4.0/3.0 * np.pi * inches_to_cm(16)**3
    fuel.add_nuclide('H2', D2O_frac * 2, 'ao')
    fuel.add_nuclide('O16', D2O_frac + 6, 'ao')

    fuel.add_nuclide('U235', 0.93, 'ao') 
    fuel.add_nuclide('U238', 0.07, 'ao')
    fuel.add_element('S', 1, 'ao')


    mats = openmc.Materials([SS304, zirconium, D2O, fuel])
    mats.cross_sections = r'/home/jakobln/devel/projects/reaktorfysik/data/jeff-3.3-hdf5/cross_sections.xml'


    clad_inner = openmc.Sphere(r = inches_to_cm(16))
    clad_outer = openmc.Sphere(r = inches_to_cm(16 + 5/16))    
    PV_inner = openmc.Sphere(r = inches_to_cm(30))              
    PV_outer = openmc.Sphere(r = inches_to_cm(30 + 4.4), boundary_type = 'vacuum')


    pressure_vessel = +PV_inner & -PV_outer                       #Stainless steal 304
    core_vessel = +clad_inner & -clad_outer                    #zirconium alloy
    moderator_region = +clad_outer & -PV_inner                 #Heavy water
    fuel_region = -clad_inner

    cell_shield    = openmc.Cell(name='blast shield', fill=SS304, region=pressure_vessel)
    cell_clad      = openmc.Cell(name='clad (Zr)', fill=zirconium, region=core_vessel)
    cell_moderator = openmc.Cell(name='moderator', fill=D2O, region=moderator_region)
    cell_fuel      = openmc.Cell(name='fuel', fill=fuel, region=fuel_region)

    root_univ = openmc.Universe(cells=[cell_shield, cell_clad, cell_moderator, cell_fuel]) 

    return openmc.Geometry(root_univ), mats

def makeModel(P, T):
    geometry, mats = makeGeometry(P, T)
    geometry.determine_paths()
    geometry.export_to_xml()

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = 5
    settings.inactive = 2
    settings.particles = 2000
    # settings.batches = 20
    # settings.inactive = 5
    # settings.particles = 20000
    settings.temperature['default'] = 600

    settings.export_to_xml()

    tally = openmc.Tally(name='keff')
    tally.estimator = 'analog'
    tally.filters = []
    tally.scores = ['kappa-fission']
    model = openmc.model.Model(
        geometry=geometry,
        settings=settings,
        materials=mats,
        tallies=openmc.Tallies([tally])
    )
    return model

def keff_T(P, T):
    model = makeModel(P, T)

    sp_path = model.run()
    sp = openmc.StatePoint(sp_path)
    keff_var = sp.keff

    mean_keff = keff_var.nominal_value
    std_keff  = keff_var.std_dev
    k_gen = sp.k_generation

    print(f"Estimated k-effective = {mean_keff:.5f} ± {std_keff:.5f}")
    return mean_keff, std_keff, k_gen



path_deplete = r"/home/jakobln/devel/projects/reaktorfysik/data/chain_casl_pwr.xml"

#operator = openmc.deplete.Operator(geom, settings, '/home/rfp/Depletion_chains/chain_casl_pwr.xml')
model = makeModel(100, 280)
operator = openmc.deplete.CoupledOperator(model, path_deplete)

powdens = 27.39726 #W/gHM; this power density corresponds to 10MWd/kgHM burnup over one year

# burnup_steps = np.array([5,5,5,5,5,5,5,5])
burnup_steps = np.array([5,10,25])

integrator = openmc.deplete.PredictorIntegrator(operator, timesteps=burnup_steps,
                                                power_density=powdens,timestep_units='MWd/kg')
integrator.integrate()

results = openmc.deplete.Results.from_hdf5("./depletion_results.h5")
time, k = results.get_keff()
time /= (24 * 60 * 60)  #seconds to days

fig, ax = plt.subplots()
ax.errorbar(time/365.25, k[:, 0], yerr=k[:, 1])
ax.set(xlabel = 'Time [yr]', ylabel = '$k_{eff}$')

_time, u5 = results.get_atoms("3", "U235",nuc_units='atom/b-cm') #we call it _time, because we already have a time variable in the correct day units which we intend to use
_time, pu239 = results.get_atoms("3", "Pu239",nuc_units='atom/b-cm')
_time, cs137 = results.get_atoms("3", "Cs137",nuc_units='atom/b-cm')
_time, xe135 = results.get_atoms("3", "Xe135",nuc_units='atom/b-cm')

fig, ax = plt.subplots()
ax.plot(time/365.25, u5, label="U235")
ax.plot(time/365.25, pu239, label="Pu239")
ax.plot(time/365.25, cs137, label="Cs137")
ax.plot(time/365.25, xe135, label="Xe135")
ax.set(xlabel = "Time [yr]", ylabel = "Atom concentration (atom/b-cm)")
ax.legend(loc = "upper right")


# print(keff_T(150,300))
# ts = np.linspace(120, 400, 29, endpoint=True)
# ps = np.array([50,100,150,170])

# k_effs = np.zeros((len(ps),len(ts)))
# with open(path + "PT dependence data.txt", "w") as file:
#     file.write("\t" + "\t".join(map(str, ts)))
#     for j,p in enumerate(ps):
#         file.write(f"\n{p}\t")
#         for i,t in enumerate(ts):
#             if(rho(p,t) < 0): 
#                 k_effs[j,i] = 0
#             else: 
#                 k_effs[j,i] = keff_T(p, t)[0]
#             file.write(f"{k_effs[j,i]:.6f}\t")
#             file.flush()


# k_effs = np.array([keff_T(t)[0] for t in ts])

# fig, ax = plt.subplots()
# ax.set(xlabel = 'Core temperature [$^o$C]', ylabel = 'k-effective', title = 'k-effective temperature dependence')
# ax.plot(ts, k_effs, marker='o', linestyle='-')
# ax.grid(True)




# keff, keff_err, k_gen = keff_T(300)
# batches = np.arange(1, len(k_gen) + 1)

# fig, ax = plt.subplots()
# ax.set(xlabel = 'Batch number', ylabel = 'k-effective per generation', title = 'Convergence of k-effective (raw)')
# ax.plot(batches, k_gen, marker='o', linestyle='-')
# ax.grid(True)


# cum_sum = np.cumsum(k_gen)
# k_cumavg = cum_sum / batches

# fig, ax = plt.subplots()
# ax.plot(batches, k_cumavg, marker='o', linestyle='-')
# ax.set(xlabel = 'Batch number', ylabel = 'Cumulative average k-effective', title = 'Convergence of k-effective (cumulative average)')
# ax.grid(True)


plt.show()