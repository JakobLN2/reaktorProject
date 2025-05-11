"""
Plot kEff vs Enrichment
This script reads the "Keffs" file containing kEff and enrichment data, processes it, and generates a plot.
"""
import numpy as np
import matplotlib.pyplot as plt

# Load the data from the file
data = np.loadtxt("git/reaktorProject/Mads/Keffs")
enrichment = data[0]
kEff = data[1]
stdKeff = data[2]
# Create a figure and axis
fig, ax = plt.subplots()
# Plot the data with error bars
ax.errorbar(enrichment, kEff, yerr=stdKeff, fmt='.', label='kEff vs Enrichment')
ax.axhline(y=1, color='r', linestyle='--', label='kEff = 1')
# Set the labels and title
ax.set_xlabel('Enrichment (% atommic occorences)')
ax.set_ylabel('kEff')
ax.set_title('kEff vs Enrichment')
# Add a grid
ax.grid()
# Add a legend
ax.legend()
# Save the figure
fig.savefig("git/reaktorProject/Mads/keffvsEnrichment")
# Show the plot
plt.show()