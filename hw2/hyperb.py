# Created by Hannah at 07.12.2023 0:23
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# из аналитики
def u_analytic(x, t):
    return np.log(1 + (t-x) ** 2) + np.sin(x)


# из устойчивости
h_kol = 21
t_kol = (h_kol - 1) * 2 - 3

h = 1/(h_kol - 1)
tau = 1/(t_kol - 1)


def show_graph(args):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    for i in args:
        ax.plot_surface(i[0], i[1], i[2])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    return


def analytics(x, t):
    x_size = x.shape[0]
    t_size = t.shape[0]

    X, T = np.meshgrid(x, t)
    u = u_analytic(X, T)

    # еще края из аналитики
    # # при t = 0
    # for l in x_size:
    #
    # # при x = 0
    # for n in t_size:

    return X, T, u


def progonka(x, t, h, tau):
    x_size = x.shape[0]
    t_size = t.shape[0]
    u = np.zeros((x_size, t_size))

    # начальные условия u_n_0 и u_n_1
    for n in range(t_size):
        u[0, n] = np.log(1 + (t[n]) ** 2)
        u[1, n] = (1 + t[n]) ** 2 + np.cos(1) + h * ((-2 * t[n])/(1 + t[n] ** 2) + np.cos(x[0])) - h ** 2 / 2 * (2 * (t[n] ** 2 - 1)/(1 + t[n] ** 2) ** 2) - np.sin(x[0])
    # начальное условие u_0_l
    for l in range(x_size):
        u[l, 0] = np.log(1 + x[l] ** 2) + np.sin(x[l])

    sigma = tau/h
    for l in range(2, x_size):
        for n in range(0, t_size - 1):
            # само уравнение переноса
            u[l, n+1] = u[l, n] + sigma/2 * (-u[l-2, n] + 4 * u[l-1, n] - 3 * u[l, n]) + (sigma ** 2)/2 * (u[l-2, n] - 2 * u[l-1, n] + u[l, n]) + tau * np.cos(x[l]) + (tau  ** 2)/ 2 * np.sin(x[l])
    return u


x = np.linspace(0, 1, h_kol)
t = np.linspace(0, 1, t_kol)
X_a, T_a, u_a = analytics(x, t)

u_p = progonka(x, t, h, tau)
u_p = u_p.transpose()

# print(X)
# print(T)
# print(u)

show_graph([[X_a, T_a, u_a]])
show_graph([[X_a, T_a, u_p]])
show_graph([[X_a, T_a, u_a], [X_a, T_a, u_p]])
