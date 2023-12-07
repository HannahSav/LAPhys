# Created by Hannah at 07.12.2023 3:02
# Created by Hannah at 07.12.2023 0:23
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# из аналитики
def u_analytic(x, t):
    print(x.shape, t.shape)
    x_size = x.shape[0]
    t_size = t.shape[0]
    # print(x)
    # print(t)


    u = np.zeros((x_size, t_size))
    for l in range(x_size):
      for n in range(t_size):
        if t[n] > np.exp(x[l]) - 1:
            u[l, n] = (np.exp(x[l]) - t[n]) * ((np.exp(x[l]) - t[n]) + np.log(np.exp(x[l]) - t[n]) - 1) + x[l] * t[n]
        if t[n] < np.exp(x[l]) - 1:
            u[l, n] = x[l] * t[n] + (-np.exp(x[l]) + t[n])**2 - 1
    return u


# из устойчивости
h_kol = 321
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

    u = u_analytic(x, t)
    X, T = np.meshgrid(x, t)


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
        u[0, n] = (t[n] * - 1) ** 2 - 1
        u[1, n] = t[n] * t[n] - 2 * t[n] - 2 * h - 2 * h * t[n] - h *h * t[n] + 2.5 *h*h
    # начальное условие u_0_l
    for l in range(x_size):
        u[l, 0] = np.exp(-x[l]) + (x[l] - 1) * np.exp(x[l])

    sigma = tau/h
    for l in range(2, x_size):
        for n in range(0, t_size - 1):
            # само уравнение переноса
            u[l, n+1] = u[l, n] + np.exp(-x[l])*(1+(tau/2)*np.exp(-x[l]))* (sigma/2) * (-u[l-2, n]+4*u[l-1, n]-3*u[l, n]) + np.exp(-2*x[l]) * ((sigma * sigma) /2 ) * (u[l-2, n] - 2*u[l-1, n] + u[l,n]) + tau * (x[l] - (tau*tau)/2 * np.exp(-x[l]))
                        # u[n + 1, l] = u[n, l] + np.exp(-x_l)*(1+(tau/2)*np.exp(-x_l))* (sigma/2) * (-u[n, l-2]+4*u[n, l-1]-3*u[n, l]) + np.exp(-2*x_l) * ((sigma * sigma) /2 ) * (u[n, l-2] - 2*u[n, l-1] + u[n,l]) + tau * (x_l - (tau*tau)/2 * np.exp(-x_l))
                                    # u[l, n+1] = u[l, n] + sigma/2 * (-u[l-2, n] + 4 * u[l-1, n] - 3 * u[l, n]) + (sigma ** 2)/2 * (u[l-2, n] - 2 * u[l-1, n] + u[l, n]) + tau * np.cos(x[l]) + (tau  ** 2)/ 2 * np.sin(x[l])
    return u


x = np.linspace(0, 1, h_kol)
t = np.linspace(0, 1, t_kol)
print(x.shape)
print(t.shape)

X_a, T_a, u_a = analytics(x, t)
u_a = u_a.transpose()

u_p = progonka(x, t, h, tau)
u_p = u_p.transpose()

# # print(X)
# # print(T)
# # print(u)

print("sdfgh", u_a.shape, u_p.shape)
show_graph([[X_a, T_a, u_a]])
show_graph([[X_a, T_a, u_p]])
show_graph([[X_a, T_a, u_a], [X_a, T_a, u_p]])