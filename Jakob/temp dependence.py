import openmc
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



SS304 = openmc.Material(name='SS304')
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
    D2O = openmc.Material(name='heavy water')
    D2O.set_density('g/cm3', rho(120, 280))
    D2O.add_nuclide('H2', 2.0, 'ao')
    D2O.add_nuclide('O16', 1.0, 'ao')



    # ——— Uranium fuel ———
    weight_frac = 10.4/1000
    D2O_frac = 235 / (20) / weight_frac

    fuel = openmc.Material(name='fuel solution')
    fuel.set_density('g/cm3', rho(P, T))
    fuel.add_nuclide('H2', D2O_frac * 2, 'ao')
    fuel.add_nuclide('O16', D2O_frac + 6, 'ao')

    fuel.add_nuclide('U235', 0.93, 'ao') 
    fuel.add_nuclide('U238', 0.07, 'ao')
    fuel.add_element('S', 1, 'ao')


    mats = openmc.Materials([SS304, D2O, zirconium, fuel])
    mats.cross_sections = r'/home/jakobln/devel/projects/reaktorfysik/data/jeff-3.3-hdf5/cross_sections.xml'


    clad_inner = openmc.Sphere(r = inches_to_cm(16))
    clad_outer = openmc.Sphere(r = inches_to_cm(16 + 5/16))    
    PV_inner = openmc.Sphere(r = inches_to_cm(30))              
    PV_outer = openmc.Sphere(r = inches_to_cm(30 + 4.4), boundary_type = 'vacuum')


    fuel_region = -clad_inner
    core_vessel = +clad_inner & -clad_outer                    #zirconium alloy
    moderator_region = +clad_outer & -PV_inner                 #Heavy water
    pressure_vessel = +PV_inner & -PV_outer                       #Stainless steal 304

    cell_fuel      = openmc.Cell(name='fuel', fill=fuel, region=fuel_region)
    cell_clad      = openmc.Cell(name='clad (Zr)', fill=zirconium, region=core_vessel)
    cell_moderator = openmc.Cell(name='moderator', fill=D2O, region=moderator_region)
    cell_shield    = openmc.Cell(name='blast shield', fill=SS304, region=pressure_vessel)

    root_univ = openmc.Universe(cells=[cell_fuel, cell_clad, cell_moderator, cell_shield]) 

    return openmc.Geometry(root_univ), mats

def makeModel(P, T):
    geometry, mats = makeGeometry(P, T)
    geometry.export_to_xml()

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    settings.batches = 20
    settings.inactive = 5
    settings.particles = 20000
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




ts = np.linspace(120, 400, 29, endpoint=True)
ps = np.array([50,100,150,170])

k_effs = np.zeros((len(ps),len(ts)))
with open(path + "PT dependence data.txt", "w") as file:
    file.write("\t" + "\t".join(map(str, ts)))
    for j,p in enumerate(ps):
        file.write(f"\n{p}\t")
        for i,t in enumerate(ts):
            if(rho(p,t) < 0): 
                k_effs[j,i] = 0
            else: 
                k_effs[j,i] = keff_T(p, t)[0]
            file.write(f"{k_effs[j,i]:.6f}\t")
            file.flush()


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