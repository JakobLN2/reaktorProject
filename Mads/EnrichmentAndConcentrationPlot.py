import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize as opt

matplotlib.use('webagg')


plt.rc("axes", labelsize=30, titlesize=32)   # skriftstørrelse af xlabel, ylabel og title
plt.rc("xtick", labelsize=26, top=True, direction="in")  # skriftstørrelse af ticks, vis også ticks øverst og vend ticks indad
plt.rc("ytick", labelsize=26, right=True, direction="in") # samme som ovenstående
plt.rc("legend", fontsize=26) # skriftstørrelse af figurers legends
plt.rcParams["font.size"] = "20"
plt.rcParams["figure.figsize"] = (16,8)



data = np.load('git/reaktorProject/Mads/KeffsEnrichmentAndConcentration2.npz')
enrichs = data['enrichments']
cons = data['concentrations']
kEffs = data['keffs']
stdKeffs = data['stdKeffs']

# Plotting
fig, ax = plt.subplots(1,1)
contour = ax.contourf(enrichs,cons,kEffs)
ax.set_xlabel("Enrichment [%U235]")
ax.set_ylabel("Concentration\n[g U pr kg D2O]")
ax.set_title(r"Fits to $\frac{\alpha}{x}$")
ax.set_xlim(0,100)
ax.set_ylim(0,100)
ax.set_aspect('equal')


Cons = kEffs*0
Cons[:] = cons
Cons = np.transpose(Cons)

Enrichs = kEffs*0
Enrichs[:] = enrichs

plotKeffs = [1,0.8,0.6]
for i in range(0,len(plotKeffs)):
    mask = (np.abs(kEffs-plotKeffs[i]) <=stdKeffs*2)
    ax.scatter(Enrichs[mask],Cons[mask],color = "C3",marker=".")
    extraMask = Enrichs[mask]<97.5
    popt, pcov = opt.curve_fit(lambda x,c: c/x,Enrichs[mask][extraMask],Cons[mask][extraMask],p0 = [1480],absolute_sigma=True,sigma=Cons[mask][extraMask]*0+1) #TODO Idk if we should 

    xs = np.linspace(0.1,100,1000)
    print(popt)
    ax.plot(xs,popt[0]/xs,color = "C1", ls = "--")
    print(Enrichs[mask])


fig.colorbar(contour, ax=ax, label='K eff')
fig.savefig("git/reaktorProject/Mads/keffvsEnrichmentDebug")

fig, ax = plt.subplots(1,1)
contour = ax.contourf(enrichs,cons,kEffs)
ax.set_xlabel("Enrichment [%U235]")
ax.set_ylabel("Concentration\n[g U pr kg Water]")
ax.set_title(r"Fits to $\frac{\alpha}{x}+\beta$")
ax.set_xlim(0,100)
ax.set_ylim(0,100)
ax.set_aspect('equal')


Cons = kEffs*0
Cons[:] = cons
Cons = np.transpose(Cons)

Enrichs = kEffs*0
Enrichs[:] = enrichs

#plotKeffs = [,1,0.8,0.6]
for i in range(0,len(plotKeffs)):
    mask = (np.abs(kEffs-plotKeffs[i]) <=stdKeffs)
    ax.scatter(Enrichs[mask],Cons[mask],color = "C3",marker=".")
    extraMask = Cons[mask]<101
    popt, pcov = opt.curve_fit(lambda x,c,d: c/x+d,Enrichs[mask][extraMask],Cons[mask][extraMask],p0 = [5,0],absolute_sigma=True,sigma=Cons[mask][extraMask]*0+1) #TODO Idk if we should 

    xs = np.linspace(0.1,100,1000)
    print(popt)
    ax.plot(xs,popt[0]/xs+popt[1],color = "C1", ls = "--")
    print(Enrichs[mask])


fig.colorbar(contour, ax=ax, label='K eff')


plt.show()
