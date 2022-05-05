import numpy as np

def Bel(array, angle_min, angle_int, Ban):
    a = array.shape[0]
    b = array.shape[1]
    for i in range(a):
        for j in range(b):
            if abs(angle_min + i * angle_int) <= 5:
                if Ban[j] < 450:
                    array[i][j] = 0.373
                if 450 <= Ban[j] <= 499:
                    array[i][j] = 0.310
                if 500 <= Ban[j] <= 549:
                    array[i][j] = 0.263
                if 550 <= Ban[j] <= 599:
                    array[i][j] = 0.224
                if 600 <= Ban[j] <= 649:
                    array[i][j] = 0.198
                if 650 <= Ban[j] <= 699:
                    array[i][j] = 0.175
                if 700 <= Ban[j] <= 749:
                    array[i][j] = 0.155
                if 750 <= Ban[j] <= 799:
                    array[i][j] = 0.139
                if Ban[j] >= 800:
                    array[i][j] = 0.125
            if 5 < abs(angle_min + i * angle_int) <= 15:
                if Ban[j] < 450:
                    array[i][j] = 0.376
                if 450 <= Ban[j] <= 499:
                    array[i][j] = 0.313
                if 500 <= Ban[j] <= 549:
                    array[i][j] = 0.265
                if 550 <= Ban[j] <= 599:
                    array[i][j] = 0.226
                if 600 <= Ban[j] <= 649:
                    array[i][j] = 0.200
                if 650 <= Ban[j] <= 699:
                    array[i][j] = 0.177
                if 700 <= Ban[j] <= 749:
                    array[i][j] = 0.157
                if 750 <= Ban[j] <= 799:
                    array[i][j] = 0.140
                if Ban[j] >= 800:
                    array[i][j] = 0.127
            if 15 < abs(angle_min + i * angle_int) <= 25:
                if Ban[j] < 450:
                    array[i][j] = 0.387
                if 450 <= Ban[j] <= 499:
                    array[i][j] = 0.322
                if 500 <= Ban[j] <= 549:
                    array[i][j] = 0.273
                if 550 <= Ban[j] <= 599:
                    array[i][j] = 0.232
                if 600 <= Ban[j] <= 649:
                    array[i][j] = 0.205
                if 650 <= Ban[j] <= 699:
                    array[i][j] = 0.182
                if 700 <= Ban[j] <= 749:
                    array[i][j] = 0.161
                if 750 <= Ban[j] <= 799:
                    array[i][j] = 0.144
                if Ban[j] >= 800:
                    array[i][j] = 0.130
            if 25 < abs(angle_min + i * angle_int) <= 35:
                if Ban[j] < 450:
                    array[i][j] = 0.405
                if 450 <= Ban[j] <= 499:
                    array[i][j] = 0.337
                if 500 <= Ban[j] <= 549:
                    array[i][j] = 0.286
                if 550 <= Ban[j] <= 599:
                    array[i][j] = 0.244
                if 600 <= Ban[j] <= 649:
                    array[i][j] = 0.216
                if 650 <= Ban[j] <= 699:
                    array[i][j] = 0.192
                if 700 <= Ban[j] <= 749:
                    array[i][j] = 0.170
                if 750 <= Ban[j] <= 799:
                    array[i][j] = 0.152
                if Ban[j] >= 800:
                    array[i][j] = 0.137
            if 35 < abs(angle_min + i * angle_int) <= 45:
                if Ban[j] < 450:
                    array[i][j] = 0.435
                if 450 <= Ban[j] <= 499:
                    array[i][j] = 0.363
                if 500 <= Ban[j] <= 549:
                    array[i][j] = 0.308
                if 550 <= Ban[j] <= 599:
                    array[i][j] = 0.264
                if 600 <= Ban[j] <= 649:
                    array[i][j] = 0.233
                if 650 <= Ban[j] <= 699:
                    array[i][j] = 0.207
                if 700 <= Ban[j] <= 749:
                    array[i][j] = 0.184
                if 750 <= Ban[j] <= 799:
                    array[i][j] = 0.163
                if Ban[j] >= 800:
                    array[i][j] = 0.149
            if 45 < abs(angle_min + i * angle_int) <= 55:
                if Ban[j] < 450:
                    array[i][j] = 0.482
                if 450 <= Ban[j] <= 499:
                    array[i][j] = 0.404
                if 500 <= Ban[j] <= 549:
                    array[i][j] = 0.343
                if 550 <= Ban[j] <= 599:
                    array[i][j] = 0.294
                if 600 <= Ban[j] <= 649:
                    array[i][j] = 0.260
                if 650 <= Ban[j] <= 699:
                    array[i][j] = 0.231
                if 700 <= Ban[j] <= 749:
                    array[i][j] = 0.205
                if 750 <= Ban[j] <= 799:
                    array[i][j] = 0.184
                if Ban[j] >= 800:
                    array[i][j] = 0.167
            if abs(angle_min + i * angle_int) > 55:
                if Ban[j] < 450:
                    array[i][j] = 0.559
                if 450 <= Ban[j] <= 499:
                    array[i][j] = 0.471
                if 500 <= Ban[j] <= 549:
                    array[i][j] = 0.401
                if 550 <= Ban[j] <= 599:
                    array[i][j] = 0.345
                if 600 <= Ban[j] <= 649:
                    array[i][j] = 0.306
                if 650 <= Ban[j] <= 699:
                    array[i][j] = 0.272
                if 700 <= Ban[j] <= 749:
                    array[i][j] = 0.242
                if 750 <= Ban[j] <= 799:
                    array[i][j] = 0.216
                if Ban[j] >= 800:
                    array[i][j] = 0.196
    return array

# Ban = [555, 645, 859]
# belta = np.zeros((1, 3))
# Bel(belta, 20, 1, Ban)
# print(belta)
