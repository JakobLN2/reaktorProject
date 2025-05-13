import openmc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as scint


plt.rc("axes", labelsize=30, titlesize=32)   # skriftstørrelse af xlabel, ylabel og title
plt.rc("xtick", labelsize=26, top=True, direction="in")  # skriftstørrelse af ticks, vis også ticks øverst og vend ticks indad
plt.rc("ytick", labelsize=26, right=True, direction="in") # samme som ovenstående
plt.rc("legend", fontsize=26) # skriftstørrelse af figurers legends
plt.rcParams["font.size"] = "20"
plt.rcParams["figure.figsize"] = (16,8)
def spherical_cap_volume(h, R):
    """
    Compute volume of a spherical cap of height h on a sphere of radius R.
    V = (π/3) * h^2 * (3R - h)
    """
    return (np.pi / 3) * h**2 * (3*R - h)

def derivative_spherical_cap_volume(h, R):
    """
    Derivative of cap volume with respect to h:
    dV/dh = (π/3) * (6Rh - 3h^2) = π * h * (2R - h)
    """
    return np.pi * h * (2*R - h)

def find_fill_height(V, R, tol=1e-6, max_iter=100):
    """
    Find the fill height h for a given volume V in a sphere of radius R
    using Newton-Raphson.
    """
    # Initial guess: approximate from V ≈ π R h^2 => h0 = sqrt(V / (π R))
    h = min(2*R, np.sqrt(V / (np.pi * R)))
    
    for i in range(max_iter):
        f = spherical_cap_volume(h, R) - V
        fp = derivative_spherical_cap_volume(h, R)
        
        # Avoid division by zero
        if fp == 0:
            break
        
        h_new = h - f / fp
        
        if abs(h_new - h) < tol:
            return h_new
        h = h_new
    
    # If not converged, return the last estimate
    return h

def inches_to_cm(inches):
    return inches * 2.54

path = r'/home/candifloos/Reaktorfysik/reaktorProject/Jakob/'
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


def makeGeometry(weight_frac_old, weight_frac_new, dt, avgdensity = False, T = 280):
        
    # ——— Heavy water moderator (D₂O) ———  
    D2O = openmc.Material(name='heavy water')
    D2O.set_density('g/cm3', rho(120, 280))
    D2O.add_nuclide('H2', 2.0, 'ao')
    D2O.add_nuclide('O16', 1.0, 'ao')



    # ——— Uranium fuel ———
    new_frac = (FR * dt)/V
    if new_frac > 1:
        new_frac = 1
    
    D2O_frac_new = (235 / (20) / weight_frac_new) * new_frac
    D2O_frac_old = (235 / (20) / weight_frac_old) * (1 - new_frac)
    D2O_frac = D2O_frac_new + D2O_frac_old

    fuel = openmc.Material(name='fuel solution')
    fuel.set_density('g/cm3', rho(100, T))
    fuel.add_nuclide('H2', D2O_frac * 2, 'ao')
    fuel.add_nuclide('O16', D2O_frac + 6, 'ao')

    fuel.add_nuclide('U235', 0.93, 'ao') 
    fuel.add_nuclide('U238', 0.07, 'ao')
    fuel.add_element('S', 1, 'ao')


    fuel_old = openmc.Material(name='fuel solution old')
    fuel_old.set_density('g/cm3', rho(100, T))
    fuel_old.add_nuclide('H2', D2O_frac_old * 2, 'ao')
    fuel_old.add_nuclide('O16', D2O_frac_old + 6, 'ao')

    fuel_old.add_nuclide('U235', 0.93, 'ao') 
    fuel_old.add_nuclide('U238', 0.07, 'ao')
    fuel_old.add_element('S', 1, 'ao')
    

    fuel_new = openmc.Material(name='fuel solution new')
    fuel_new.set_density('g/cm3', rho(100, T))
    fuel_new.add_nuclide('H2', D2O_frac_new * 2, 'ao')
    fuel_new.add_nuclide('O16', D2O_frac_new + 6, 'ao')

    fuel_new.add_nuclide('U235', 0.93, 'ao') 
    fuel_new.add_nuclide('U238', 0.07, 'ao')
    fuel_new.add_element('S', 1, 'ao')


    """Det viser sig at det ikke er sjovt at løse tredjegrads ligninger analytisk"""
    # b = 3 * inches_to_cm(16) 
    # d = 3 * V / np.pi
    # print(np.sqrt(4 * b**3 * d + 27 *d**2))

    # e = (-(3 * np.sqrt(3) * np.sqrt(4 * b**3 * d + 27 *d**2) - 2*b**3 - 27 * d) / 2)**(1/3)
    # print(e)
    # h = 1/3 * e + e**(-1) * b**2 + b

    # y = -inches_to_cm(16) + h
    # print(y)
    h = find_fill_height(V = FR * dt, R = inches_to_cm(16))
    y = -inches_to_cm(16) + h

    mats = openmc.Materials([SS304, D2O, zirconium, fuel_old, fuel_new, fuel])
    mats.cross_sections = r'/home/candifloos/Reaktorfysik/Python/RFP/Data/jeff-3.3-hdf5/cross_sections.xml'


    clad_inner = openmc.Sphere(r = inches_to_cm(16))
    clad_outer = openmc.Sphere(r = inches_to_cm(16 + 5/16))    
    PV_inner = openmc.Sphere(r = inches_to_cm(30))              
    PV_outer = openmc.Sphere(r = inches_to_cm(30 + 4.4), boundary_type = 'vacuum')
    split = openmc.YPlane(y)

    fuel_region_old = -clad_inner & -split
    fuel_region_new = -clad_inner & +split
    core_vessel = +clad_inner & -clad_outer                    #zirconium alloy
    moderator_region = +clad_outer & -PV_inner                 #Heavy water
    pressure_vessel = +PV_inner & -PV_outer                       #Stainless steal 304

    if avgdensity:
        cell_fuel_old  = openmc.Cell(name='fuel_old', fill=fuel, region=fuel_region_old)
        cell_fuel_new  = openmc.Cell(name='fuel_new', fill=fuel, region=fuel_region_new)   
    else:
        cell_fuel_old  = openmc.Cell(name='fuel_old', fill=fuel_old, region=fuel_region_old)
        cell_fuel_new  = openmc.Cell(name='fuel_new', fill=fuel_new, region=fuel_region_new)
    
    cell_clad      = openmc.Cell(name='clad (Zr)', fill=zirconium, region=core_vessel)
    cell_moderator = openmc.Cell(name='moderator', fill=D2O, region=moderator_region)
    cell_shield    = openmc.Cell(name='blast shield', fill=SS304, region=pressure_vessel)

    root_univ = openmc.Universe(cells=[cell_fuel_old, cell_fuel_new, cell_clad, cell_moderator, cell_shield]) 

    root_univ.plot()

    return openmc.Geometry(root_univ), mats

def makeModel(weight_frac_old, weight_frac_new, dt, avgdensity = False, T = 280):
    geometry, mats = makeGeometry(weight_frac_old, weight_frac_new, dt, avgdensity, T)
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

def keff_T(weight_frac_old, weight_frac_new, dt, avgdensity = False, T = 280):
    model = makeModel(weight_frac_old, weight_frac_new, dt, avgdensity, T)

    sp_path = model.run()
    sp = openmc.StatePoint(sp_path)
    keff_var = sp.keff

    mean_keff = keff_var.nominal_value
    std_keff  = keff_var.std_dev
    k_gen = sp.k_generation

    print(f"Estimated k-effective = {mean_keff:.5f} ± {std_keff:.5f}")
    return mean_keff, std_keff, k_gen




V = 4 / 3 * np.pi * (inches_to_cm(16))**3     #cm^3
FR = 25236                      #cm^3/s, 400 gps
old = 10.4 /1000
new = 10 /1000


dts = np.linspace(1, V / FR, 5)
k_effs = np.array([keff_T(weight_frac_old = old, weight_frac_new = new, dt = dt, avgdensity = False)[0] for dt in dts])
# k_effs_2 = np.array([keff_T(weight_frac_old = old, weight_frac_new = new, dt = dt, avgdensity = True)[0] for dt in dts])

fig, ax = plt.subplots()
ax.set(xlabel = 'Time [s]', ylabel = 'k-effective', title = f'k-effective time dependence, {round(old * 100, 2)} -> {round(new * 100, 2)} % uranium concentration')
ax.plot(dts, k_effs, marker='o', linestyle='-')
# ax.plot(dts, k_effs_2, marker='o', linestyle='-')
ax.grid(True)




plt.show()