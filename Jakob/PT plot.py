"""Make a contour plot of the density of heavy water"""
import openmc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.interpolate as scint
matplotlib.use('webagg')

#Small table - Below 275 C
T = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/badWater.txt",max_rows=1)
P = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/badWater.txt",usecols=[0],skiprows=1)
data = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/badWater.txt",skiprows=1,usecols=range(1,len(T)+1))
f_1 = scint.RegularGridInterpolator((P,T), data, method="linear")

#Russian table - above or equal 275 C
T = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/russianWater.txt",max_rows=1)
P = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/russianWater.txt",usecols=[0],skiprows=1) * 0.9807
data = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/russianWater.txt",skiprows=1,usecols=range(1,len(T)+1))
f_2 = scint.RegularGridInterpolator((P,T), data, method="linear")
rho_2 = lambda P,T: 1/f_2((P,T))

rho = lambda P,T: f_1((P,T))/1000 if T < 275 else rho_2(P,T)

def rho_plot(P,T):
    try: 
        val = rho(P,T)
        # return val if val >= 0 else 0
        return val
    except ValueError: return 0

T_plot = np.linspace(3.82,425, 500)
P_plot = np.linspace(1,340*0.9807, 500) 

Z_plot = np.zeros((len(T_plot), len(P_plot)))
for i,t in enumerate(T_plot):
    for j,p in enumerate(P_plot):
        Z_plot[i,j] = rho_plot(p,t)

fig, ax = plt.subplots()
ax.set(xlabel="Temperature [C]", ylabel="Pressure [bar]",title="Heavy water density against pressure and temperature")
t,p = np.meshgrid(T_plot, P_plot, indexing="ij")

levels = np.linspace(0.0, 1.3, 40)
cb = ax.contourf(t,p,Z_plot,levels)
fig.colorbar(cb)
plt.gca().invert_yaxis()
# fig.savefig("../D2O PT graph.png")
plt.show()