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


path = r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/"

ts = np.loadtxt(path + "PT dependence data.txt", max_rows=1)
ps = np.loadtxt(path + "PT dependence data.txt", usecols=[0],skiprows=1)
ks = np.loadtxt(path + "PT dependence data.txt", skiprows=1, usecols=range(1,len(ts) + 1))

fig, ax = plt.subplots()
ax.set(xlabel="Temperature [C]", ylabel="$K_{eff}$", title="Temperature dependence on k-effective value")
for i,k_t in enumerate(ks):
    filt = np.logical_and(0 < k_t, k_t < 1.2)
    ax.plot(ts[filt], k_t[filt],marker="",label=f"{ps[i]} bar")

ax.legend(loc = "lower left")
ax.grid()

fig, ax = plt.subplots()
ax.set(xlabel="Temperature [C]", ylabel="$K_{eff}$", title="Temperature dependence on k-effective value @ 170 bar",
       xlim=(ts[0], 350), ylim=(0.8,1.2))
ax.plot(ts, ks[-1],marker="",label=f"{ps[-1]} bar", color="red")
ax.legend(loc = "lower left")
ax.grid()
plt.show()
