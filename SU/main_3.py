from math import cos, sin, pi

import numpy as np

d = 30.0/2                         # Базовый диаметр кулака (мм)
h = 12.0                         # Максимальное перемещение толкателя (мм)
z = 0.25                         # Тепловой зазор (мм)
f_pod = 35.0 / 180 * pi          # Фаза подъёма (град)
f_v = 5.0 / 180 * pi             # Фаза выдержки (град)
f_op = 70.0 / 180 * pi           # Фаза опускания (град)
f_z = 2 / 180 * pi               # Фаза теплового зазора (град)

fi_0 = 0
fi_1 = fi_0 + f_z / 2
fi_2 = fi_1 + f_z / 2
fi_3 = fi_2 + f_pod
fi_4 = fi_3 + f_pod
fi_5 = fi_4 + f_z

# Участки
#1 - fi_01
#2 - fi_12
#3 - fi_23
#4 - fi_34
#5 - fi_45

# Граничные условия для участков

delta1_Y = {
    'fi_0':0,
    'fi_5':0,
}
fun_Y = {
    'fi_0':d - z,
    'fi_1':d - z / 2,
    'fi_2':d,
    'fi_3':d + h / 2,
    'fi_4':d + h,
    'fi_5':d - z,
}

def delta2_fun(x, fi_1, fi_0, C, D):
    return -C * cos(pi * (x - fi_0) / (fi_1 - fi_0)) + D

def delta1_fun(x, fi_1, fi_0, C, D, const_1):
    return -C * sin(pi * (x - fi_0) / (fi_1 - fi_0)) / (pi / (fi_1 - fi_0)) + D * x + const_1

def const_1(fi_1, fi_0, C, D, Y):
    x = fi_0
    return Y - (-C * sin(pi * (x - fi_0) / (fi_1 - fi_0)) / (pi / (fi_1 - fi_0)) + D * x)

def fun(x, fi_1, fi_0, C, D, const_1, const_2):
    return C * cos(pi * (x - fi_0) / (fi_1 - fi_0)) / ((pi / (fi_1 - fi_0)) ** 2) + (1 / 2) * D * x**2 + const_1 * x + const_2

def const_2(fi_1, fi_0, C, D, const_1, Y):
    x = fi_0
    return Y - (C * cos(pi * (x - fi_0) / (fi_1 - fi_0)) / ((pi / (fi_1 - fi_0)) ** 2) + (1 / 2) * D * x**2 + const_1 * x)

def C_fun(Y_1, Y_0):
    return 0.5 * (Y_1 - Y_0)

def D_fun(Y_1, Y_0):
    return 0.5 * (Y_1 + Y_0)

'''
def find_line(Y_1, Y_0, fun_Y_0, fun_Y_2, delta1_Y_0, delta2_Y_0, fi_0, fi_1):
    fun_Y_1 = (fun_Y_2 + fun_Y_0) / 2
    while True:
'''

delta2_Y = {
    'fi_0':0,
    'fi_1':100,
    'fi_2':-10,
    'fi_3':40,
    'fi_4':-50,
    'fi_5':0,
}
C = {
    'fi_01':C_fun(delta2_Y['fi_1'], delta2_Y['fi_0']),
    'fi_12':C_fun(delta2_Y['fi_2'], delta2_Y['fi_1']),
    'fi_23':C_fun(delta2_Y['fi_3'], delta2_Y['fi_2']),
    'fi_34':C_fun(delta2_Y['fi_4'], delta2_Y['fi_3']),
    'fi_45':C_fun(delta2_Y['fi_5'], delta2_Y['fi_4']),
}
D = {
    'fi_01':D_fun(delta2_Y['fi_1'], delta2_Y['fi_0']),
    'fi_12':D_fun(delta2_Y['fi_2'], delta2_Y['fi_1']),
    'fi_23':D_fun(delta2_Y['fi_3'], delta2_Y['fi_2']),
    'fi_34':D_fun(delta2_Y['fi_4'], delta2_Y['fi_3']),
    'fi_45':D_fun(delta2_Y['fi_5'], delta2_Y['fi_4']),
}

# Участок 1
const_dict = {
    'fi_01':[],
    'fi_12':[],
    'fi_23':[],
    'fi_34':[],
    'fi_45':[],
}
const_dict['fi_01'].append(const_1(fi_1, fi_0, C['fi_01'], D['fi_01'], delta1_Y['fi_0']))
const_dict['fi_01'].append(const_2(fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], fun_Y['fi_0']))
delta1_Y['fi_1'] = delta1_fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0])
print(fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1]))
fun_Y['fi_1'] = fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1])
const_dict['fi_12'].append(const_1(fi_2, fi_1, C['fi_12'], D['fi_12'], delta1_Y['fi_1']))
const_dict['fi_12'].append(const_2(fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], fun_Y['fi_1']))
delta1_Y['fi_2'] = delta1_fun(fi_2, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0])
print(fun(fi_2, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1]))

const_dict['fi_23'].append(const_1(fi_3, fi_2, C['fi_23'], D['fi_23'], delta1_Y['fi_2']))
const_dict['fi_23'].append(const_2(fi_3, fi_2, C['fi_23'], D['fi_23'], const_dict['fi_01'][0], fun_Y['fi_2']))

delta1_Y['fi_3'] = delta1_fun(fi_3, fi_3, fi_2, C['fi_23'], D['fi_23'], const_dict['fi_23'][0])
fun_Y['fi_3'] = fun(fi_3, fi_3, fi_2, C['fi_23'], D['fi_23'], const_dict['fi_23'][0], const_dict['fi_23'][1])
print(fun(fi_3, fi_3, fi_2, C['fi_23'], D['fi_23'], const_dict['fi_23'][0], const_dict['fi_23'][1]))
const_dict['fi_34'].append(const_1(fi_4, fi_3, C['fi_34'], D['fi_34'], delta1_Y['fi_3']))
const_dict['fi_34'].append(const_2(fi_4, fi_3, C['fi_34'], D['fi_34'], const_dict['fi_34'][0], fun_Y['fi_3']))
delta1_Y['fi_4'] = delta1_fun(fi_4, fi_4, fi_3, C['fi_34'], D['fi_34'], const_dict['fi_34'][0])
print(fun(fi_4, fi_4, fi_3, C['fi_34'], D['fi_34'], const_dict['fi_34'][0], const_dict['fi_34'][1]))

X_1 = np.linspace(0, fi_1, 100)
Y_1 = []
for i in X_1:
    Y_1.append(fun(i, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1]))

X_2 = np.linspace(fi_1, fi_2, 100)
Y_2 = []
for i in X_2:
    Y_2.append(fun(i, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1]))

X_3 = np.linspace(fi_2, fi_3, 100)
Y_3 = []
for i in X_3:
    Y_3.append(fun(i, fi_3, fi_2, C['fi_23'], D['fi_23'], const_dict['fi_23'][0], const_dict['fi_23'][1]))

X_4 = np.linspace(fi_3, fi_4, 100)
Y_4 = []
for i in X_4:
    Y_4.append(fun(i, fi_4, fi_3, C['fi_34'], D['fi_34'], const_dict['fi_34'][0], const_dict['fi_34'][1]))

import matplotlib.pyplot as plt
plt.plot(X_1, Y_1)
plt.plot(X_2, Y_2)
plt.ylim(29.75, 30)
plt.show()

plt.plot(X_3, Y_3)
plt.plot(X_4, Y_4)
plt.show()