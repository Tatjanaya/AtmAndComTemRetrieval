import numpy as np

t_8_n = 300
t_8_o = 295
t_9_n = 301
t_9_o = 297

c1 = 3.7404e8   
c2 = 14387

band1 = 10.852
band2 = 12

# l_t_8_n = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / t_8_n) - 1))
# l_t_8_o = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / t_8_o) - 1))
# l_t_9_n = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / t_9_n) - 1))
# l_t_9_o = c1 / (np.pi * pow(band2, 5) * (np.exp(c2 / band2 / t_9_o) - 1))

# print(l_t_8_n)
# print(l_t_8_o)
# print(l_t_9_n)
# print(l_t_9_o)

# L_s8_0 = c1 / (np.pi * pow(band1, 5) * (np.exp(c2 / band1 / T_ang_1) - 1))
# c2 / band1 / np.log((c1 / (np.pi * pow(band1, 5) * t_8_n))+1)

di_n_8 = 1.0174 * t_8_n + 0.6045 * (t_8_n - t_9_n) + 0.2652 * (t_8_n - t_9_n)**2 + -3.7037
di_o_8 = 1.0223 * t_8_o + 0.9217 * (t_8_o - t_9_o) + 0.2249 * (t_8_o - t_9_o)**2 + -4.3483
di_n_9 = 1.0292 * t_8_n + -0.0667 * (t_8_n - t_9_n) + 0.4339 * (t_8_n - t_9_n)**2 + -6.7047
di_o_9 = 1.035 * t_8_o + 0.3792 * (t_8_o - t_9_o) + 0.3414 * (t_8_o - t_9_o)**2 + -7.5637

# di_n_8 = 1.0269 * t_8_n + 0.6116 * (t_8_n - t_9_n) + 0.2294 * (t_8_n - t_9_n)**2 + -6.2953
# di_o_8 = 1.0296 * t_8_o + 1.0276 * (t_8_o - t_9_o) + 0.1855 * (t_8_o - t_9_o)**2 + -6.3340
# di_n_9 = 1.0489 * t_8_n + -0.1371 * (t_8_n - t_9_n) + 0.3874 * (t_8_n - t_9_n)**2 + -12.1453
# di_o_9 = 1.0534 * t_8_o + 0.4243 * (t_8_o - t_9_o) + 0.2982 * (t_8_o - t_9_o)**2 + -12.5983

print(di_n_8)
print(di_o_8)
print(di_n_9)
print(di_o_9)