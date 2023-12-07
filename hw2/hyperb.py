# Created by Hannah at 07.12.2023 0:23
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# из аналитики
def u_analytic(x, t):
    return np.log(1 + (t-x) ** 2) + np.sin(x)


def show_graph(args, h, tau):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.plot_surface(args[0][0], args[0][1], args[0][2], label='Аналитическое решение')
    ax.plot_surface(args[1][0], args[1][1], args[1][2], label='Рассчетное решение')
    ax.set_xlabel('X')
    ax.set_ylabel('T')
    ax.set_zlabel('U')
    plt.legend()
    plt.title(f'U(x, t) для различного шагов h={h} и tau={tau}')
    plt.show()
    return


def show_graph_pl(x, t, u):
    #отрисовка u(x)
    plt.figure(figsize=(8, 6))
    u_new = u
    for i in range((t).shape[0]):
        plt.plot(x, u_new[i], label=f'u(x) при t = {int(i*tau*100)/100}')
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend()
    plt.title('Аналитическое решение u(x)')
    plt.show()

    # отрисовка u(t)
    #plt.figure(figsize=(8, 6))
    u_new = u.transpose()
    for i in range(x.shape[0]):
        plt.plot(t, u_new[i], label=f'u(t) при x = {int(i*h/100)*100}')
    plt.xlabel('t')
    plt.ylabel('u')
    plt.legend()
    plt.title('Аналитическое решение u(t)')
    plt.show()
    return


def show_last_layer(x_array, y1_array, y2_array, h):
    for i in range(len(x_array)):
        plt.plot(x_array[i], y1_array[i], label=f'U_analytics(x, t) при h={h[i]}')
        plt.plot(x_array[i], y2_array[i], label= f'U_progonka(x, t) при h={h[i]}')
    plt.xlabel('x')
    plt.ylabel('U(x, t=1)')
    plt.legend()
    plt.title('Слои аналитического решения и рассчетной сетки для U(x, t) при t = 1 для различного шага h')
    plt.show()
    return


def show_log(x, y):
    plt.plot(x, y, marker='o')
    plt.xlabel('log(h)')
    plt.ylabel('log(diff)')
    plt.title('Зависимость log(diff) от log(h)')
    plt.show()
    return


def analytics(x, t):

    X, T = np.meshgrid(x, t)
    u = u_analytic(X, T)

    return X, T, u


def progonka(x, t, h, tau):
    x_size = x.shape[0]
    t_size = t.shape[0]
    u = np.zeros((x_size, t_size))

    # начальные условия u_n_0 и u_n_1
    for n in range(t_size):
        u[0, n] = np.log(1 + (t[n]) ** 2)
        u[1, n] = u[0, n] + h * (np.cos(x[0]) - 2 * t[n] / (1 + t[n] ** 2)) + h ** 2 / 2 * (- np.sin(x[0]) + 2 * (1 - t[n]) / ((1 + t[n] ** 2) ** 2) - np.sin(x[0]))

    # начальное условие u_0_l
    for l in range(x_size):
        u[l, 0] = np.log(1 + x[l] ** 2) + np.sin(x[l])

    sigma = tau/h
    for l in range(2, x_size):
        for n in range(0, t_size - 1):
            # само уравнение переноса
            u[l, n+1] = u[l, n] + sigma/2 * (-u[l-2, n] + 4 * u[l-1, n] - 3 * u[l, n]) + (sigma ** 2)/2 * (u[l-2, n] - 2 * u[l-1, n] + u[l, n]) + tau * np.cos(x[l]) + (tau  ** 2)/ 2 * np.sin(x[l])
    return u

def count_t_h(h_kol):
    h = 1 / (h_kol - 1)
    tau = (2 * h) * 0.95
    t_kol = int(1 / tau) + 1
    tau = 1 / (t_kol - 1)
    return t_kol, h, tau


def max_diff_log(u_a, u_p):
    diff = 0
    for i in range(len(u_a)):
        diff = max(diff, abs(u_a[i] - u_p[i]))
    return np.log(diff)

h_kol = 3
t_kol, h, tau = count_t_h(h_kol)

u_p_ending = []
u_a_ending = []
x_ending = []
h_array = []
diff_log = []
h_log = []

for i in range(5):
    x = np.linspace(0, 1, h_kol)
    t = np.linspace(0, 1, t_kol)
    X_a, T_a, u_a = analytics(x, t)
    u_p = progonka(x, t, h, tau)
    u_p = u_p.transpose()

    u_a_ending.append(u_a[-1].tolist())
    u_p_ending.append(u_p[-1].tolist())
    x_ending.append(x.tolist())
    h_array.append(h)

    diff_log.append(max_diff_log(u_a_ending[-1], u_p_ending[-1]))
    h_log.append(np.log(h))

    show_graph([[X_a, T_a, u_a], [X_a, T_a, u_p]], h, tau)
    show_last_layer(x_ending, u_a_ending, u_p_ending, h_array)

    h_kol = (h_kol - 1) * 2 + 1
    t_kol, h, tau = count_t_h(h_kol)

show_log(h_log, diff_log)

