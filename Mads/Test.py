import numpy as np

x = np.linspace(0, 10, 100)
y = np.linspace(0, 10, 100)

xx,yy = np.meshgrid(x, y)

z = np.sin(xx) * np.cos(yy)

np.savez('data.npz', x=x, y=y, z=z)

npzfile = np.load('data.npz')
x = npzfile['x']
y = npzfile['y']
z = npzfile['z']
print(x)
print(y)
print(z)
