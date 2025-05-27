import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as scint
import json
import scipy.optimize as sc


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
linear = lambda x,*p: p[0]*x + p[1]

with open(path + "depletion results.txt", "r") as file:
    data = json.load(file)
    time = np.array(data["t"])
    k = data["k"]; kerr = data["kerr"]
    V_fuel = data["fuel volume"]; pow = data["power"]
    u235 = np.array(data["u235"]); cs137 = np.array(data["cs137"]); xe135 = np.array(data["xe135"])

idx = np.argwhere(time >= 0.5)[0][0]
p_opt, _ = sc.curve_fit(linear, time[:idx], u235[:idx] * 235/6.022e23, p0=[0,0])
u_rate = p_opt[0]
print(f"Rate of U235 consumption @ {pow} MW: {-u_rate*1e-3:.3f} kg/yr")


fig, ax = plt.subplots()
ax.errorbar(time, k, yerr=kerr, marker="o", capsize=3)
ax.set(xlabel = 'Time [yr]', ylabel = '$k_{eff}$', title="$k_{eff}$ over time", xlim=(0,1))
ax.grid()
plt.savefig(path + "Depletion keff depletion.png", bbox_inches="tight")

fig, ax = plt.subplots()
ax.plot(time, u235/V_fuel, label="U235")
ax.plot(time, cs137/V_fuel, label="Cs137")
ax.plot(time, xe135/V_fuel, label="Xe135")
ax.plot(time[:idx + 3], linear(time[:idx + 3], *p_opt)*6.022e23/235/V_fuel, color="k", linestyle="--", label="linear fit")
ax.set(xlabel = "Time [yr]", ylabel = "Atom concentration (atom/cm$^3$)", title="Fuel region populations", yscale="linear", xlim=(0,2))
ax.legend(loc = "upper right")
ax.grid()
plt.savefig(path + "Depletion fuel populations.png", bbox_inches="tight")

fig, ax = plt.subplots()
ax.plot(time[:-1], np.diff(u235)/np.diff(time)/V_fuel, label="U235")
ax.set(xlabel = "Time [yr]", ylabel = "Change in atom concentration (atom/cm$^3$/yr)", title="Change in U235 concentration", yscale="linear")
ax.legend(loc = "upper right")
plt.savefig(path + "depletion rate.png", bbox_inches="tight")

plt.show()