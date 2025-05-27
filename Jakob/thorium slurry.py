import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as scint
import json

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
rho_ThO2 = 10 #g/cm^3

mass_ThO2 = 0

geom_dic = {"rho_ThO2": rho_ThO2}



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


th_weight = 0.6 #kg Th / kg D2O
thO2_weight = (1 + 32/232) * th_weight #kg ThO2 / kg D2O
rho_tot = lambda p,t: (1 + thO2_weight)/(1/rho(p,t) + thO2_weight/rho_ThO2) #assumes mixing with V = V1 + V2 (yikes)
geom_dic["rho_tot"] = rho_tot(120,280)

# ——— Heavy water moderator (D₂O) ———  
D2O = openmc.Material(3, name='moderator slurry')
D2O.set_density('g/cm3', rho_tot(120, 280))
D2O.volume = 4.0/3.0 * np.pi * (inches_to_cm(30)**3 - inches_to_cm(16 + 5/16)**3)
D2O.add_nuclide('H2', 4/20, 'wo')
D2O.add_nuclide('O16', 16/20 + 32/232*th_weight, 'wo')
D2O.add_element('Th', th_weight, 'wo')
geom_dic["blanket_volume"] = D2O.volume

mass_moderator = D2O.volume * rho_tot(120, 280)
mass_ThO2 = mass_moderator * th_weight/(th_weight + 32/232 + 1)
geom_dic["m_ThO2"] = mass_ThO2
geom_dic["m_blanket"] = mass_moderator

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
    fuel.depletable = False
    geom_dic["fuel_volume"] = fuel.volume


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

# print(keff_T(150,300)) #For some reason this is required when the relevant xml files arent created? integrator apparantly wants the files before it starts


path_deplete = r"/home/jakobln/devel/projects/reaktorfysik/data/chain_casl_pwr.xml"

#operator = openmc.deplete.Operator(geom, settings, '/home/rfp/Depletion_chains/chain_casl_pwr.xml')
model = makeModel(120, 280)
model.export_to_model_xml(path + "model.xml")
operator = openmc.deplete.CoupledOperator(model, path_deplete)

pow = 5e6 #W
geom_dic["power"] = pow/1e6


burnup_steps = np.array([1e-4, 0.001, 0.001, 0.001, 0.001, 0.01, 0.01, 0.01, 0.01, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.5, 0.5, 0.5, 0.5, 2, 2, 2])
# burnup_steps = np.array([0.1, 0.1, 0.1])

integrator = openmc.deplete.PredictorIntegrator(operator, timesteps=burnup_steps,
                                                power=pow,timestep_units='a')
integrator.integrate()
results = openmc.deplete.Results.from_hdf5("./depletion_results.h5")
time, k = results.get_keff()
time /= (24 * 60 * 60 * 365.25)  #seconds to years


# _time, th232 = results.get_atoms("3", "Th232",nuc_units='atom/b-cm') 
# _time, u235 = results.get_atoms("3", "U235",nuc_units='atom/b-cm') #we call it _time, because we already have a time variable in the correct day units which we intend to use
# _time, u233 = results.get_atoms("3", "U233",nuc_units='atom/b-cm') 
# _time, pu239 = results.get_atoms("3", "Pu239",nuc_units='atom/b-cm')
# _time, cs137 = results.get_atoms("3", "Cs137",nuc_units='atom/b-cm')
# _time, xe135 = results.get_atoms("3", "Xe135",nuc_units='atom/b-cm')

depletion_results = {"t": list(time), "k": list(k[:,0]), "kerr": list(k[:,1])}
nucs = ["Th232", "U235", "U233", "Pu239", "Cs137", "Xe135"]
for nuc in nucs:
    depletion_results[nuc.lower()] = list(results.get_atoms("3", nuc,nuc_units='atoms')[1])
depletion_results.update(geom_dic)
with open(path + "thorium breeding results.txt", "w") as file:
    json.dump(depletion_results, file, indent=3, ensure_ascii=False)


# print(f"mass of thorium in blanket: {mass_ThO2:.2f} g")
# gamma = th232[0]/mass_ThO2
# ms_u233 = u233/gamma*1e-3

# fig, ax = plt.subplots()
# ax.errorbar(time, k[:, 0], yerr=k[:, 1])
# ax.set(xlabel = 'Time [yr]', ylabel = '$k_{eff}$')
# plt.savefig(path + "keff depletion.png", bbox_inches="tight")

# fig, ax = plt.subplots()
# ax.plot(time, th232, label="Th232")
# ax.plot(time, u235, label="U235")
# ax.plot(time, u233, label="U233")
# ax.plot(time, pu239, label="Pu239")
# ax.plot(time, cs137, label="Cs137")
# ax.plot(time, xe135, label="Xe135")
# ax.set(xlabel = "Time [yr]", ylabel = "Atom concentration (atom/b-cm)", title="Blanket population", yscale="linear")
# ax.legend(loc = "upper right")
# plt.savefig(path + "breeding populations.png", bbox_inches="tight")


# fig, ax = plt.subplots()
# ax.plot(time, ms_u233, label="U233")
# ax.set(xlabel = "Time [yr]", ylabel = "Mass [kg]", title="Bred mass", yscale="linear")
# ax.legend(loc = "lower right")
# plt.savefig(path + "bred mass.png", bbox_inches="tight")


# _time, u235 = results.get_atoms("4", "U235",nuc_units='atom/b-cm') #we call it _time, because we already have a time variable in the correct day units which we intend to use
# _time, cs137 = results.get_atoms("4", "Cs137",nuc_units='atom/b-cm')
# _time, xe135 = results.get_atoms("4", "Xe135",nuc_units='atom/b-cm')

# fig, ax = plt.subplots()
# ax.plot(time, u235, label="U235")
# ax.plot(time, cs137, label="Cs137")
# ax.plot(time, xe135, label="Xe135")
# ax.set(xlabel = "Time [yr]", ylabel = "Atom concentration (atom/b-cm)", title="Blanket population", yscale="linear")
# ax.legend(loc = "upper right")
# plt.savefig(path + "fuel populations.png", bbox_inches="tight")


# plt.show()
