# import matplotlib
# import matplotlib.pyplot as plt
# matplotlib.use('TkAgg')
# print(plt.get_backend())

# fig, ax = plt.subplots()
# ax.plot([0,1], [1,0])
# plt.show()

# arr = [0] * 2
# print(arr)

import numpy as np
N = 4
b = -5; t = 5
ls = np.linspace(b,t,N+1, endpoint=True)
T_grad = lambda y: (y - b) * (300 - 256)/(t - b) + 256
print(ls)
Ts_mid = [0]*N
for i in range(N):
    Ts_mid[i] = T_grad(1/2*(t - b)/N+ls[i])
print(Ts_mid)