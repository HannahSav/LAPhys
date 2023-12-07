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


def show_graph_pl(x, t, u):
    #отрисовка u(x)
    plt.figure(figsize=(8, 6))
    u_new = u
    print(x.shape, u_new[0].shape)
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
    print('asd', t.shape, u_new[0].shape)
    for i in range(x.shape[0]):
        plt.plot(t, u_new[i], label=f'u(t) при x = {int(i*h/100)*100}')
    plt.xlabel('t')
    plt.ylabel('u')
    plt.legend()
    plt.title('Аналитическое решение u(t)')
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








h_kol1 = 161
t_kol1 = (h_kol1 - 1) * 2 - 3

h1 = 1/(h_kol1 - 1)
tau1 = 1/(t_kol1 - 1)
x1 = np.linspace(0, 1, h_kol1)
t1 = np.linspace(0, 1, t_kol1)
X_a1, T_a1, u_a1 = analytics(x1, t1)


u_p1 = progonka(x1, t1, h1, tau1)
u_p1 = u_p1.transpose()




print(u_p.shape)

# show_graph([[X_a, T_a, u_a]])
# show_graph([[X_a, T_a, u_p]])
show_graph([[X_a, T_a, u_p], [X_a1, T_a1, u_p1]])
# show_graph_pl(x, t, u_a)
