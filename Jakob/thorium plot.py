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

with open(path + "thorium breeding results.txt", "r") as file:
    data = json.load(file)
    time = np.array(data["t"])
    k = data["k"]; kerr = data["kerr"]
    V_fuel = data["fuel_volume"]; pow = data["power"]
    mass_ThO2 = data["m_ThO2"]
    # u235 = np.array(data["u235"]); cs137 = np.array(data["cs137"]); xe135 = np.array(data["xe135"])

# print(f"mass of thorium in blanket: {mass_ThO2:.2f} g")
# gamma = data["th232"][0]/mass_ThO2
# ms_u233 = np.array(data["u233"])/gamma*1e-3
ms_u233 = np.array(data["u233"])*233/6.022e23*1e-3 #kg

# if(1.0 in time):
#     idx = np.argwhere(time == 1.0)[0]
#     print(f"mass of U233 after 1 year @ {data["power"]:.1f} MW : {:.2f} kg")
# else: 
ms_u233_yr = scint.interp1d(time, ms_u233, kind="cubic")(1.0)
print(f"mass of U233 after 1 year @ {data["power"]:.1f} MW : {ms_u233_yr:.2f} kg")
#     print("bad")


fig, ax = plt.subplots()
ax.errorbar(time, data["k"], yerr=data["kerr"])
ax.set(xlabel = 'Time [yr]', ylabel = '$k_{eff}$', yscale="linear")
ax.grid()
plt.savefig(path + "keff depletion.png", bbox_inches="tight")

fig, ax = plt.subplots()
ax.plot(time, np.array(data["th232"])/data["blanket_volume"], label="Th232")
ax.plot(time, np.array(data["u235"])/data["blanket_volume"], label="U235")
ax.plot(time, np.array(data["u233"])/data["blanket_volume"], label="U233")
ax.plot(time, np.array(data["pu239"])/data["blanket_volume"], label="Pu239")
ax.plot(time, np.array(data["cs137"])/data["blanket_volume"], label="Cs137")
ax.plot(time, np.array(data["xe135"])/data["blanket_volume"], label="Xe135")
ax.set(xlabel = "Time [yr]", ylabel = "Atom concentration (atom/b-cm)", title="Blanket population", yscale="log", ylim=(1e7), xlim=(-0.1,2))
ax.legend(loc = "lower right")
ax.grid()
plt.savefig(path + "breeding populations.png", bbox_inches="tight")


fig, ax = plt.subplots()
ax.plot(time, ms_u233, label="U233", color="tab:green")
ax.set(xlabel = "Time [yr]", ylabel = "Mass [kg]", title="Bred mass", yscale="linear", xlim=(0,2))
ax.legend(loc = "lower right")
ax.grid()
plt.savefig(path + "bred mass.png", bbox_inches="tight")


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


plt.show()
