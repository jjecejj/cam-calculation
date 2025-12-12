from math import cos, sin, pi

import numpy as np

d = 30.0                         # Базовый диаметр кулака (мм)
h = 12.0                         # Максимальное перемещение толкателя (мм)
z = 0.25                         # Тепловой зазор (мм)
f_pod = 75.0 / 180 * pi          # Фаза подъёма (град)
f_v = 5.0 / 180 * pi             # Фаза выдержки (град)
f_op = 70.0 / 180 * pi           # Фаза опускания (град)
f_z = 2 / 180 * pi               # Фаза теплового зазора (град)

fi_0 = 0
fi_1 = fi_0 + f_z
fi_2 = fi_1 + f_pod
fi_3 = fi_2 + f_v
fi_4 = fi_3 + f_op
fi_5 = fi_4 + f_z

# Участки
#1 - fi_01
#2 - fi_12
#3 - fi_23
#4 - fi_34
#5 - fi_45

# Граничные условия для участков
delta3_Y = {
    'fi_0':0,
    'fi_5':0,
}
delta2_Y = {
    'fi_0':0,
    'fi_5':0,
}
delta1_Y = {
    'fi_0':0,
    'fi_5':0,
}
fun_Y = {
    'fi_0':d - z,
    'fi_1':d,
    'fi_2':d + h,
    'fi_3':d + h,
    'fi_4':d,
    'fi_5':d - z,
}

def delta4_fun(x, fi_1, fi_0, C, D):
    'Четвёртая производная'
    return -C * cos(pi * (x - fi_0) / (fi_1 - fi_0)) + D

def delta3_fun(x, fi_1, fi_0, C, D, const_1):
    'Третья производная (рывок)'
    return -C * sin(pi * (x - fi_0) / (fi_1 - fi_0)) / (pi / (fi_1 - fi_0)) + D * x + const_1

def const_1(fi_1, fi_0, C, D, Y):
    x = fi_0
    return Y - (-C * sin(pi * (x - fi_0) / (fi_1 - fi_0)) / (pi / (fi_1 - fi_0)) + D * x)

def delta2_fun(x, fi_1, fi_0, C, D, const_1, const_2):
    'Вторая производная (ускорение)'
    return C * cos(pi * (x - fi_0) / (fi_1 - fi_0)) / ((pi / (fi_1 - fi_0)) ** 2) + (1 / 2) * D * x**2 + const_1 * x + const_2

def const_2(fi_1, fi_0, C, D, const_1, Y):
    x = fi_0
    return Y - (C * cos(pi * (x - fi_0) / (fi_1 - fi_0)) / ((pi / (fi_1 - fi_0)) ** 2) + (1 / 2) * D * x**2 + const_1 * x)

def delta1_fun(x, fi_1, fi_0, C, D, const_1, const_2, const_3):
    'Первая производная (скорость)'
    return C * sin(pi * (x - fi_0) / (fi_1 - fi_0)) / ((pi / (fi_1 - fi_0)) ** 3) + (1 / 6) * D * x**3 + (1 / 2) * const_1 * x**2 + const_2 * x + const_3

def const_3(fi_1, fi_0, C, D, const_1, const_2, Y):
    x = fi_0
    return Y - (C * sin(pi * (x - fi_0) / (fi_1 - fi_0)) / ((pi / (fi_1 - fi_0)) ** 3) + (1 / 6) * D * x**3 + (1 / 2) * const_1 * x**2 + const_2 * x)

def fun(x, fi_1, fi_0, C, D, const_1, const_2, const_3, const_4):
    'Функция точки контакта кулачка и толкателя'
    return -C * cos(pi * (x - fi_0) / (fi_1 - fi_0)) / ((pi / (fi_1 - fi_0)) ** 4) + (1 / 24) * D * x ** 4 + (1 / 6) * const_1 * x ** 3 + (1 / 2) * const_2 * x**2 + const_3 * x + const_4

def const_4(fi_1, fi_0, C, D, const_1, const_2, const_3, Y):
    x = fi_0
    return Y - (-C * cos(pi * (x - fi_0) / (fi_1 - fi_0)) / ((pi / (fi_1 - fi_0)) ** 4) + (1 / 24) * D * x ** 4 + (1 / 6) * const_1 * x ** 3 + (1 / 2) * const_2 * x**2 + const_3 * x)

def C_fun(Y_1, Y_0):
    return 0.5 * (Y_1 - Y_0)

def D_fun(Y_1, Y_0):
    return 0.5 * (Y_1 + Y_0)

# Значение Четвёртой производной при соответсвующих углах
'''
min = [1e8, 0]
for i in np.linspace(-264677552.5871959 - 100, -264677552.5871959+ 100, 100000):
    print("I:", i)
    delta4_Y = {
        'fi_0':0,
        'fi_1':40000000,
        'fi_2':i,
        'fi_3':0,
        'fi_4':5,
        'fi_5':0,
    }
    C = {
        'fi_01':C_fun(delta4_Y['fi_1'], delta4_Y['fi_0']),
        'fi_12':C_fun(delta4_Y['fi_2'], delta4_Y['fi_1']),
        'fi_23':None,
        'fi_34':C_fun(delta4_Y['fi_4'], delta4_Y['fi_3']),
        'fi_45':C_fun(delta4_Y['fi_5'], delta4_Y['fi_4']),
    }
    D = {
        'fi_01':D_fun(delta4_Y['fi_1'], delta4_Y['fi_0']),
        'fi_12':D_fun(delta4_Y['fi_2'], delta4_Y['fi_1']),
        'fi_23':None,
        'fi_34':D_fun(delta4_Y['fi_4'], delta4_Y['fi_3']),
        'fi_45':D_fun(delta4_Y['fi_5'], delta4_Y['fi_4']),
    }

    # Участок 1
    const_dict = {
        'fi_01':[],
        'fi_12':[],
        'fi_23':None,
        'fi_34':[],
        'fi_45':[],
    }
    const_dict['fi_01'].append(const_1(fi_1, fi_0, C['fi_01'], D['fi_01'], delta3_Y['fi_0']))
    const_dict['fi_01'].append(const_2(fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], delta2_Y['fi_0']))
    const_dict['fi_01'].append(const_3(fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], delta1_Y['fi_0']))
    const_dict['fi_01'].append(const_4(fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], const_dict['fi_01'][2], fun_Y['fi_0']))
    #print(delta4_fun(0, fi_1, fi_0, C['fi_01'], D['fi_01']))
    delta3_Y['fi_1'] = delta3_fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0])
    delta2_Y['fi_1'] = delta2_fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1])
    delta1_Y['fi_1'] = delta1_fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], const_dict['fi_01'][2])
    #print(fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], const_dict['fi_01'][2], const_dict['fi_01'][3]))

    const_dict['fi_12'].append(const_1(fi_2, fi_1, C['fi_12'], D['fi_12'], delta3_Y['fi_1']))
    const_dict['fi_12'].append(const_2(fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], delta2_Y['fi_1']))
    const_dict['fi_12'].append(const_3(fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], delta1_Y['fi_1']))
    const_dict['fi_12'].append(const_4(fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], const_dict['fi_12'][2], fun_Y['fi_1']))
    #print(delta4_fun(fi_1, fi_2, fi_1, C['fi_12'], D['fi_12']))
    #print(delta3_fun(fi_1, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0]))
    #print(delta2_fun(fi_1, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1]))
    #print(delta1_fun(fi_1, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], const_dict['fi_12'][2]))
    temp = fun(fi_2, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], const_dict['fi_12'][2], const_dict['fi_12'][3])
    #print(fun(fi_2, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], const_dict['fi_12'][2], const_dict['fi_12'][3]))
    print(temp)
    if min[0] > abs(temp - (d + h)):
        min[0] = abs(temp - (d + h))
        min[1] = i
print(min)
exit()
'''
delta4_Y = {
    'fi_0':0,
    'fi_1':40000000,
    'fi_2':-264677552.28219286,
    'fi_3':0,
    'fi_4':5,
    'fi_5':0,
}
C = {
    'fi_01':C_fun(delta4_Y['fi_1'], delta4_Y['fi_0']),
    'fi_12':C_fun(delta4_Y['fi_2'], delta4_Y['fi_1']),
    'fi_23':None,
    'fi_34':C_fun(delta4_Y['fi_4'], delta4_Y['fi_3']),
    'fi_45':C_fun(delta4_Y['fi_5'], delta4_Y['fi_4']),
}
D = {
    'fi_01':D_fun(delta4_Y['fi_1'], delta4_Y['fi_0']),
    'fi_12':D_fun(delta4_Y['fi_2'], delta4_Y['fi_1']),
    'fi_23':None,
    'fi_34':D_fun(delta4_Y['fi_4'], delta4_Y['fi_3']),
    'fi_45':D_fun(delta4_Y['fi_5'], delta4_Y['fi_4']),
}

# Участок 1
const_dict = {
    'fi_01':[],
    'fi_12':[],
    'fi_23':None,
    'fi_34':[],
    'fi_45':[],
}
const_dict['fi_01'].append(const_1(fi_1, fi_0, C['fi_01'], D['fi_01'], delta3_Y['fi_0']))
const_dict['fi_01'].append(const_2(fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], delta2_Y['fi_0']))
const_dict['fi_01'].append(const_3(fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], delta1_Y['fi_0']))
const_dict['fi_01'].append(const_4(fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], const_dict['fi_01'][2], fun_Y['fi_0']))
print(delta4_fun(0, fi_1, fi_0, C['fi_01'], D['fi_01']))
delta3_Y['fi_1'] = delta3_fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0])
delta2_Y['fi_1'] = delta2_fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1])
delta1_Y['fi_1'] = delta1_fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], const_dict['fi_01'][2])
print(fun(fi_1, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], const_dict['fi_01'][2], const_dict['fi_01'][3]))

const_dict['fi_12'].append(const_1(fi_2, fi_1, C['fi_12'], D['fi_12'], delta3_Y['fi_1']))
const_dict['fi_12'].append(const_2(fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], delta2_Y['fi_1']))
const_dict['fi_12'].append(const_3(fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], delta1_Y['fi_1']))
const_dict['fi_12'].append(const_4(fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], const_dict['fi_12'][2], fun_Y['fi_1']))
print(delta4_fun(fi_1, fi_2, fi_1, C['fi_12'], D['fi_12']))
print(delta3_fun(fi_1, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0]))
print(delta2_fun(fi_1, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1]))
print(delta1_fun(fi_1, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], const_dict['fi_12'][2]))
print(fun(fi_2, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], const_dict['fi_12'][2], const_dict['fi_12'][3]))

X_1 = np.linspace(0, fi_1, 100)
Y_1 = []
for i in X_1:
    Y_1.append(fun(i, fi_1, fi_0, C['fi_01'], D['fi_01'], const_dict['fi_01'][0], const_dict['fi_01'][1], const_dict['fi_01'][2], const_dict['fi_01'][3]))

X_2 = np.linspace(fi_1, fi_2, 100)
Y_2 = []
for i in X_2:
    Y_2.append(fun(i, fi_2, fi_1, C['fi_12'], D['fi_12'], const_dict['fi_12'][0], const_dict['fi_12'][1], const_dict['fi_12'][2], const_dict['fi_12'][3]))

import matplotlib.pyplot as plt
plt.plot(X_1, Y_1)
plt.plot(X_2, Y_2)
plt.show()