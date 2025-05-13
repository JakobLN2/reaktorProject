import numpy as np
import matplotlib.pyplot as plt

data = np.load('git/reaktorProject/Mads/KeffsEnrichmentAndConcentration.npz')
enrichs = data['enrichments']
cons = data['concentrations']
kEffs = data['keffs']
stdKeffs = data['stdKeffs']

# Plotting
fig, ax = plt.subplots(1,1)
contour = ax.contourf(enrichs,cons,kEffs)
ax.set_xlabel("Enrichment")
ax.set_ylabel("Concentration")
ax.set_title("K eff vs Enrichment and Concentration")
ax.set_xlim(0,100)
ax.set_ylim(0,100)
ax.set_aspect('equal')
fig.colorbar(contour, ax=ax, label='K eff')
fig.savefig("git/reaktorProject/Mads/keffvsEnrichmentDebug")
plt.show()
