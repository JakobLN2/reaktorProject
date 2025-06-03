import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as scint
import scipy.optimize as sc


matplotlib.use('webagg')
plt.rc("axes", labelsize=30, titlesize=32)   # skriftstørrelse af xlabel, ylabel og title
plt.rc("xtick", labelsize=26, top=True, direction="in")  # skriftstørrelse af ticks, vis også ticks øverst og vend ticks indad
plt.rc("ytick", labelsize=26, right=True, direction="in") # samme som ovenstående
plt.rc("legend", fontsize=26) # skriftstørrelse af figurers legends
plt.rcParams["font.size"] = "20"
plt.rcParams["figure.figsize"] = (16,8)


path = r"git/reaktorProject/Mads/Geometry Investigation"

dataNames = ["BasicBall","LongPipes","ShortPipes"]


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

f_lin = lambda x, *p: p[0] * x + p[1]
filt = np.logical_and(280 <= ts, ts <= 330)
# p_opt, p_cov = sc.curve_fit(f_lin, ts[filt], ks[-1][filt], sigma=ks_err[-1][filt], absolute_sigma=False, p0 = [0,0])
p_opt, p_cov = sc.curve_fit(f_lin, ts[filt], ks[-1][filt], absolute_sigma=False, p0 = [0,0])
ax.plot(ts, f_lin(ts, *p_opt), label="linear fit around 300$^o$C", color="k")
print(f"temperature coefficient rho = {p_opt[0]:.6f} +/- {np.sqrt(p_cov[0][0]):.6f}")

ax.legend(loc = "lower left")
ax.grid()

savepath = r'/home/candifloos/Reaktorfysik/Figurer//'
plt.savefig(savepath + 'Temperature coefficient.png', bbox_inches='tight')

plt.show()
