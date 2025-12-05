from math import sin, cos, pi
import numpy as np

# -----------------------------
# Исходные данные
# -----------------------------

N_k = 1000                                # Количество оборотов кулачка в минут
T = 60 / N_k
d = 30.0 * 1e-3                         # Базовый диаметр кулака (мм)
r0 = d / 2
h = 12.0 * 1e-3                         # Максимальное перемещение толкателя (мм)
z = 0.25 * 1e-3                         # Тепловой зазор (мм)
f_pod = 80.0                     # Фаза подъёма (град)
f_v = 5.0                        # Фаза выдержки (град)
f_op = 75.0                      # Фаза опускания (град)
f_z = 25                         # Фаза теплового зазора (град)

# Диапазон значений от 1 до 100, принимают только целочисленные величины
k_1 = 15                         # коэффициент агрессивности первого участка (выбор зазора)
k_2 = 15                         # коэффициент агрессивности второго участка (Фаза подъёма)
k_3 = 15                         # коэффициент агрессивности пятого участка (Фаза опускания)
k_4 = 15                         # коэффициент агрессивности шестого участка (Фаза выбора зазора)

phi_0 = 0
phi_1 = phi_0 + f_z
phi_2 = phi_1 + f_pod
phi_3 = phi_2 + f_v
phi_4 = phi_3 + f_op
phi_5 = phi_4 + f_z
# -----------------------------
# Профилирование методом "полидайн"
# -----------------------------
p = 5
q = p + p - 2
r = q + p - 2
s = r + p - 2

# -----------------------------
# Коэффициенты C2, Cp, Cq, Cr, Cs
# -----------------------------

C2 = -p * q * r * s / ((p - 2) * (q - 2) * (r - 2) * (s - 2))
Cp = 2 * q * r * s / ((p - 2) * (q - p) * (r - p) * (s - p))
Cq = -2 * p * s * r / ((q - 2) * (r - q) * (q - p) * (s - q))
Cr = 2 * p * q * s / ((r - 2) * (s - r) * (r - p) * (r - q))
Cs = -2 * p * q * r / ((s - 2) * (s - p) * (s - q) * (s - r))

def C_fun(p):
    q = p + p - 2
    r = q + p - 2
    s = r + p - 2
    C2 = -p * q * r * s / ((p - 2) * (q - 2) * (r - 2) * (s - 2))
    Cp = 2 * q * r * s / ((p - 2) * (q - p) * (r - p) * (s - p))
    Cq = -2 * p * s * r / ((q - 2) * (r - q) * (q - p) * (s - q))
    Cr = 2 * p * q * s / ((r - 2) * (s - r) * (r - p) * (r - q))
    Cs = -2 * p * q * r / ((s - 2) * (s - p) * (s - q) * (s - r))
    return [C2, Cp, Cq, Cr, Cs], [2, p, q, r, s]

def h_phi(fi, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h):
    'Функция перемещения от угла поворота кулачка'
    temp = 1
    for i in range(0, len(C_list)):
        temp += C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**k_list[i])
    return temp * h_kn_max

def v_phi(fi, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h):
    'Функция скорости от угла поворота кулачка'
    temp = 1
    for i in range(0, len(C_list)):
        temp += k_list[i] * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**(k_list[i] - 1))
    return temp * h_kn_max

def a_phi(fi, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h):
    'Функция ускорения от угла поворота кулачка'
    temp = 1
    for i in range(0, len(C_list)):
        temp += k_list[i] * (k_list[i] - 1) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**(k_list[i] - 2))
    return temp * h_kn_max

def d_phi(fi, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h):
    'Функция рывка от угла поворота кулачка'
    temp = 1
    for i in range(0, len(C_list)):
        if C_list[i] <= 2:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**(k_list[i] - 3))
    return temp * h_kn_max

def k_phi(fi, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h):
    'Функция кракена от угла поворота кулачка'
    temp = 1
    for i in range(0, len(C_list)):
        if C_list[i] <= 3:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * (k_list[i] - 3) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**(k_list[i] - 4))
    return temp * h_kn_max

def h_t(t, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h, T = T):
    'Функция перемещения от времени'
    temp = 1
    fi = t / T * 360
    for i in range(0, len(C_list)):
        temp += C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**k_list[i])
    return temp * h_kn_max

def v_t(t, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h, T = T):
    'Функция скорости от времени'
    temp = 1
    fi = t / T * 360
    for i in range(0, len(C_list)):
        temp += k_list[i] * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**(k_list[i] - 1))
    return temp * h_kn_max

def a_t(t, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h, T = T):
    'Функция ускорения от времени'
    temp = 1
    fi = t / T * 360
    for i in range(0, len(C_list)):
        temp += k_list[i] * (k_list[i] - 1) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**(k_list[i] - 2))
    return temp * h_kn_max

def d_t(t, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h, T = T):
    'Функция рывка от времени'
    temp = 1
    fi = t / T * 360
    for i in range(0, len(C_list)):
        if C_list[i] <= 2:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**(k_list[i] - 3))
    return temp * h_kn_max

def k_t(t, C_list, k_list, fi_1, fi_0 = 0, h_kn_max = h, T = T):
    'Функция кракена от времени'
    temp = 1
    fi = t / T * 360
    for i in range(0, len(C_list)):
        if C_list[i] <= 3:
            continue
        temp += k_list[i] * (k_list[i] - 1) * (k_list[i] - 2) * (k_list[i] - 3) * C_list[i] * (((fi - fi_0) / (fi_1 - fi_0))**(k_list[i] - 4))
    return temp * h_kn_max

C_list_1, k_list_1  = C_fun(1 + round(k_1))
C_list_2, k_list_2  = C_fun(1 + round(k_2))
C_list_3, k_list_3  = C_fun(1 + round(k_3))
C_list_4, k_list_4  = C_fun(1 + round(k_4))

N = 1000
X_1 = np.linspace(0, phi_1, N)
Y_1 = []
V_1 = []
A_1 = []
D_1 = []
K_1 = []
for i in X_1:
    Y_1.append(-h_phi(i, C_list_1, k_list_1, phi_1, fi_0=0, h_kn_max = z) + r0)
    V_1.append(-v_phi(i, C_list_1, k_list_1, phi_1, fi_0=0, h_kn_max = z) + r0)
    A_1.append(-a_phi(i, C_list_1, k_list_1, phi_1, fi_0=0, h_kn_max = z) + r0)
    D_1.append(-d_phi(i, C_list_1, k_list_1, phi_1, fi_0=0, h_kn_max = z) + r0)
    K_1.append(-k_phi(i, C_list_1, k_list_1, phi_1, fi_0=0, h_kn_max = z) + r0)

X_2 = np.linspace(phi_1, phi_2, N)
Y_2 = []
V_2 = []
A_2 = []
D_2 = []
K_2 = []
for i in X_2:
    Y_2.append(-h_phi(i, C_list_2, k_list_2, phi_2, fi_0=phi_1, h_kn_max=h) + r0 + h)
    V_2.append(-v_phi(i, C_list_2, k_list_2, phi_2, fi_0=phi_1, h_kn_max=h) + r0 + h)
    A_2.append(-a_phi(i, C_list_2, k_list_2, phi_2, fi_0=phi_1, h_kn_max=h) + r0 + h)
    D_2.append(-d_phi(i, C_list_2, k_list_2, phi_2, fi_0=phi_1, h_kn_max=h) + r0 + h)
    K_2.append(-k_phi(i, C_list_2, k_list_2, phi_2, fi_0=phi_1, h_kn_max=h) + r0 + h)

X_3 = np.linspace(phi_3, phi_4, N)
Y_3 = []
V_3 = []
A_3 = []
D_3 = []
K_3 = []
for i in X_3:
    Y_3.append(-h_phi(phi_4 - i + phi_3, C_list_3, k_list_3, phi_4, fi_0=phi_3, h_kn_max = h) + r0 + h)
    V_3.append(v_phi(phi_4 - i + phi_3, C_list_3, k_list_3, phi_4, fi_0=phi_3, h_kn_max = h) + r0 + h)
    A_3.append(-a_phi(phi_4 - i + phi_3, C_list_3, k_list_3, phi_4, fi_0=phi_3, h_kn_max = h) + r0 + h)
    D_3.append(-d_phi(phi_4 - i + phi_3, C_list_3, k_list_3, phi_4, fi_0=phi_3, h_kn_max = h) + r0 + h)
    K_3.append(-k_phi(phi_4 - i + phi_3, C_list_3, k_list_3, phi_4, fi_0=phi_3, h_kn_max = h) + r0 + h)

X_4 = np.linspace(phi_4, phi_5, N)
Y_4 = []
V_4 = []
A_4 = []
D_4 = []
K_4 = []
for i in X_4:
    Y_4.append(-h_phi(phi_5 - i + phi_4, C_list_4, k_list_4, phi_5, fi_0=phi_4, h_kn_max = z) + r0)
    V_4.append(v_phi(phi_5 - i + phi_4, C_list_4, k_list_4, phi_5, fi_0=phi_4, h_kn_max = z) + r0)
    A_4.append(-a_phi(phi_5 - i + phi_4, C_list_4, k_list_4, phi_5, fi_0=phi_4, h_kn_max = z) + r0)
    D_4.append(-d_phi(phi_5 - i + phi_4, C_list_4, k_list_4, phi_5, fi_0=phi_4, h_kn_max = z) + r0)
    K_4.append(-k_phi(phi_5 - i + phi_4, C_list_4, k_list_4, phi_5, fi_0=phi_4, h_kn_max = z) + r0)

t_1 = X_1 * (T / 360)
t_2 = X_2 * (T / 360)
t_3 = X_3 * (T / 360)
t_4 = X_4 * (T / 360)

# Зависимости от грудуса поворота кулачка
'''V_1_grad = np.gradient(Y_1, X_1[1] - X_1[0])
V_2_grad = np.gradient(Y_2, X_2[1] - X_2[0])
V_3_grad = np.gradient(Y_3, X_3[1] - X_3[0])
V_4_grad = np.gradient(Y_4, X_4[1] - X_4[0])

A_1_grad = np.gradient(V_1_grad, X_1[1] - X_1[0])
A_2_grad = np.gradient(V_2_grad, X_2[1] - X_2[0])
A_3_grad = np.gradient(V_3_grad, X_3[1] - X_3[0])
A_4_grad = np.gradient(V_4_grad, X_4[1] - X_4[0])

D_1_grad = np.gradient(A_1_grad, X_1[1] - X_1[0])
D_2_grad = np.gradient(A_2_grad, X_2[1] - X_2[0])
D_3_grad = np.gradient(A_3_grad, X_3[1] - X_3[0])
D_4_grad = np.gradient(A_4_grad, X_4[1] - X_4[0])

K_1_grad = np.gradient(D_1_grad, X_1[1] - X_1[0])
K_2_grad = np.gradient(D_2_grad, X_2[1] - X_2[0])
K_3_grad = np.gradient(D_3_grad, X_3[1] - X_3[0])
K_4_grad = np.gradient(D_4_grad, X_4[1] - X_4[0])

# Зависимости от времени
V_1_sek = np.gradient(Y_1, t_1[1] - t_1[0])
V_2_sek = np.gradient(Y_2, t_2[1] - t_2[0])
V_3_sek = np.gradient(Y_3, t_3[1] - t_3[0])
V_4_sek = np.gradient(Y_4, t_4[1] - t_4[0])

A_1_sek = np.gradient(V_1_sek, t_1[1] - t_1[0])
A_2_sek = np.gradient(V_2_sek, t_2[1] - t_2[0])
A_3_sek = np.gradient(V_3_sek, t_3[1] - t_3[0])
A_4_sek = np.gradient(V_4_sek, t_4[1] - X_4[0])

D_1_sek = np.gradient(A_1_sek, t_1[1] - t_1[0])
D_2_sek = np.gradient(A_2_sek, t_2[1] - t_2[0])
D_3_sek = np.gradient(A_3_sek, t_3[1] - t_3[0])
D_4_sek = np.gradient(A_4_sek, t_4[1] - t_4[0])

K_1_sek = np.gradient(D_1_sek, t_1[1] - t_1[0])
K_2_sek = np.gradient(D_2_sek, t_2[1] - t_2[0])
K_3_sek = np.gradient(D_3_sek, t_3[1] - t_3[0])
K_4_sek = np.gradient(D_4_sek, t_4[1] - t_4[0])'''

import matplotlib.pyplot as plt

fig, axs = plt.subplots(5, 1, figsize=(8, 20))  # 5 строк, 1 столбец
'''
# --- 1. Координата ---
axs[0].plot(t_1, Y_1)
axs[0].plot(t_2, Y_2)
axs[0].plot([t_2[-1], t_3[0]], [r0 + h, r0 + h])
axs[0].plot(t_3, Y_3)
axs[0].plot(t_4, Y_4)
axs[0].set_xlabel('t, c')
axs[0].set_ylabel('координата толкателя, м')
axs[0].grid(True)

# --- 2. Скорость ---
axs[1].plot(t_1, V_1_sek)
axs[1].plot(t_2, V_2_sek)
axs[1].plot([t_2[-1], t_3[0]], [0, 0])
axs[1].plot(t_3, V_3_sek)
axs[1].plot(t_4, V_4_sek)
axs[1].set_xlabel('t, c')
axs[1].set_ylabel('Скорость')
axs[1].grid(True)

# --- 3. Ускорение ---
axs[2].plot(t_1, A_1_sek)
axs[2].plot(t_2, A_2_sek)
axs[2].plot([t_2[-1], t_3[0]], [0, 0])
axs[2].plot(t_3, A_3_sek)
axs[2].plot(t_4, A_4_sek)
axs[2].set_xlabel('t, c')
axs[2].set_ylabel('Ускорение')
axs[2].grid(True)

# --- 4. Рывок ---
axs[3].plot(t_1, D_1_sek)
axs[3].plot(t_2, D_2_sek)
axs[3].plot([t_2[-1], t_3[0]], [0, 0])
axs[3].plot(t_3, D_3_sek)
axs[3].plot(t_4, D_4_sek)
axs[3].set_xlabel('t, c')
axs[3].set_ylabel('Рывок')
axs[3].grid(True)

# --- 5. "Кракен" ---
axs[4].plot(t_1, K_1_sek)
axs[4].plot(t_2, K_2_sek)
axs[4].plot([t_2[-1], t_3[0]], [0, 0])
axs[4].plot(t_3, K_3_sek)
axs[4].plot(t_4, K_4_sek)
axs[4].set_xlabel('t, c')
axs[4].set_ylabel('Кракен')
axs[4].grid(True)

plt.tight_layout()
plt.show()
'''

# --- 1. Координата ---
axs[0].plot(t_1, Y_1)
axs[0].plot(t_2, Y_2)
axs[0].plot([t_2[-1], t_3[0]], [r0 + h, r0 + h])
axs[0].plot(t_3, Y_3)
axs[0].plot(t_4, Y_4)
axs[0].set_xlabel('t, c')
axs[0].set_ylabel('координата толкателя, м')
axs[0].grid(True)

# --- 2. Скорость ---
axs[1].plot(t_1, V_1)
axs[1].plot(t_2, V_2)
axs[1].plot([t_2[-1], t_3[0]], [0, 0])
axs[1].plot(t_3, V_3)
axs[1].plot(t_4, V_4)
axs[1].set_xlabel('t, c')
axs[1].set_ylabel('Скорость')
axs[1].grid(True)

# --- 3. Ускорение ---
axs[2].plot(t_1, A_1)
axs[2].plot(t_2, A_2)
axs[2].plot([t_2[-1], t_3[0]], [0, 0])
axs[2].plot(t_3, A_3)
axs[2].plot(t_4, A_4)
axs[2].set_xlabel('t, c')
axs[2].set_ylabel('Ускорение')
axs[2].grid(True)

# --- 4. Рывок ---
axs[3].plot(t_1, D_1)
axs[3].plot(t_2, D_2)
axs[3].plot([t_2[-1], t_3[0]], [0, 0])
axs[3].plot(t_3, D_3)
axs[3].plot(t_4, D_4)
axs[3].set_xlabel('t, c')
axs[3].set_ylabel('Рывок')
axs[3].grid(True)

# --- 5. "Кракен" ---
axs[4].plot(t_1, K_1)
axs[4].plot(t_2, K_2)
axs[4].plot([t_2[-1], t_3[0]], [0, 0])
axs[4].plot(t_3, K_3)
axs[4].plot(t_4, K_4)
axs[4].set_xlabel('t, c')
axs[4].set_ylabel('Кракен')
axs[4].grid(True)

plt.tight_layout()
plt.show()

X = []
Y = []
for i in range(0, N):
    X.append(Y_1[i] * sin(X_1[i] / 180 * pi))
    Y.append(Y_1[i] * cos(X_1[i] / 180 * pi))
for i in range(0, N):
    X.append(Y_2[i] * sin(X_2[i] / 180 * pi))
    Y.append(Y_2[i] * cos(X_2[i] / 180 * pi))
for i in np.linspace(phi_2, phi_3, N):
    X.append((r0 + h) * sin(i / 180 * pi))
    Y.append((r0 + h) * cos(i / 180 * pi))
for i in range(0, N):
    X.append(Y_3[i] * sin(X_3[i] / 180 * pi))
    Y.append(Y_3[i] * cos(X_3[i] / 180 * pi))
for i in range(0, N):
    X.append(Y_4[i] * sin(X_4[i] / 180 * pi))
    Y.append(Y_4[i] * cos(X_4[i] / 180 * pi))
for i in np.linspace(phi_5, 360, N):
    X.append((r0 - z) * sin(i / 180 * pi))
    Y.append((r0 - z) * cos(i / 180 * pi))

plt.plot(X, Y)
plt.scatter([0], [0])
plt.xlabel('X, м')
plt.ylabel('Y, м')
plt.grid(True)
plt.title("Профиль кулачка")
plt.show()