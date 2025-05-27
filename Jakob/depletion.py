import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as scint


matplotlib.use('TkAgg')
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



# ——— Heavy water moderator (D₂O) ———  
D2O = openmc.Material(3, name='moderator slurry')
D2O.set_density('g/cm3', rho(120, 280))
D2O.volume = 4.0/3.0 * np.pi * (inches_to_cm(30)**3 - inches_to_cm(16 + 5/16)**3)
D2O.add_nuclide('H2', 2, 'wo')
D2O.add_nuclide('O16', 1, 'wo')

def makeGeometry(P, T):
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
    fuel.depletable = True


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

print(keff_T(150,300))


path_deplete = r"/home/jakobln/devel/projects/reaktorfysik/data/chain_casl_pwr.xml"

#operator = openmc.deplete.Operator(geom, settings, '/home/rfp/Depletion_chains/chain_casl_pwr.xml')
model = makeModel(120, 280)
model.export_to_model_xml(path + "model.xml")
operator = openmc.deplete.CoupledOperator(model, path_deplete)

pow = 1e6 #W


burnup_steps = np.array([1e-4, 0.001, 0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.2,0.2,0.2,0.2,0.2,0.2,0.2])
# burnup_steps = np.array([0.1, 0.1, 0.1])

integrator = openmc.deplete.PredictorIntegrator(operator, timesteps=burnup_steps,
                                                power=pow,timestep_units='a')
integrator.integrate()
results = openmc.deplete.Results.from_hdf5("./depletion_results.h5")
time, k = results.get_keff()
time /= (24 * 60 * 60* 365.25)  #seconds to years

fig, ax = plt.subplots()
ax.errorbar(time, k[:, 0], yerr=k[:, 1])
ax.set(xlabel = 'Time [yr]', ylabel = '$k_{eff}$')
plt.savefig(path + "Depletion keff depletion.png", bbox_inches="tight")

_time, u235 = results.get_atoms("4", "U235",nuc_units='atom/b-cm') #we call it _time, because we already have a time variable in the correct day units which we intend to use
_time, cs137 = results.get_atoms("4", "Cs137",nuc_units='atom/b-cm')
_time, xe135 = results.get_atoms("4", "Xe135",nuc_units='atom/b-cm')

fig, ax = plt.subplots()
ax.plot(time, u235, label="U235")
ax.plot(time, cs137, label="Cs137")
ax.plot(time, xe135, label="Xe135")
ax.set(xlabel = "Time [yr]", ylabel = "Atom concentration (atom/b-cm)", title="Fuel populations", yscale="linear")
ax.legend(loc = "upper right")
plt.savefig(path + "Depletion fuel populations.png", bbox_inches="tight")

fig, ax = plt.subplots()
ax.plot(time[:-1], np.diff(u235)/np.diff(time), label="U235")
ax.set(xlabel = "Time [yr]", ylabel = "Change in atom concentration (atom/b-cm/yr)", title="Change in U235 concentration", yscale="linear")
ax.legend(loc = "upper right")
plt.savefig(path + "depletion rate.png", bbox_inches="tight")

plt.show()
