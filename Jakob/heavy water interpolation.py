import numpy as np
import scipy.interpolate as scint


#Small table - Below 275 C
T = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/badWater.txt",max_rows=1)
P = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/badWater.txt",usecols=[0],skiprows=1)
data = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/badWater.txt",skiprows=1,usecols=range(1,len(T)+1))
f_1 = scint.RegularGridInterpolator((P,T), data, method="cubic")

#Russian table - above or equal 275 C
T = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/russianWater.txt",max_rows=1)
P = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/russianWater.txt",usecols=[0],skiprows=1)
data = np.loadtxt(r"/home/jakobln/devel/projects/reaktorfysik/reaktorProject/Jakob/russianWater.txt",skiprows=1,usecols=range(1,len(T)+1))
f_2 = scint.RegularGridInterpolator((P,T), data, method="linear")
rho_2 = lambda P,T: 1/f_2((P,T))

rho = lambda P,T: f_1((P,T))/1000 if T < 275 else rho_2(P,T)


num = 0
for i,p in enumerate(P):
    for j,t in enumerate(T):
        if(abs(1/rho(p,t) - data[i,j]) > 1e-4):
            num += 1
            print(rho(p,t), data[i,j])
            # exit()
print(np.prod(np.shape(data)))
print(f"{num = }")

# print(rho(100,280),rho_2(100,280))
# print(rho(100, 250))