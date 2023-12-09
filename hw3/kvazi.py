# Created by Hannah at 07.12.2023 16:21
import numpy as np
import matplotlib.pyplot as plt

# из расчетной сетки
NU = 2
C0 = 5
C1 = 1

EPS = 10**(-4)


def u_analytics(t, x):
    return (C1 + x) ** (2/NU) * (C0 - 2 * (2 + NU)/NU * t) ** (-1/NU)


def show_graph(x_a, y_a, z_a, z_p, h, tau):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.plot_surface(x_a, y_a, z_a, label='Аналитическое решение')
    ax.plot_surface(x_a, y_a, z_p, label='Рассчетное решение')
    ax.set_xlabel('X')
    ax.set_ylabel('T')
    ax.set_zlabel('U')
    plt.legend()
    plt.title(f'U(x, t) для различного шагов h={int(h*10000)/10000} и tau={int(tau*10000)/10000}')
    plt.show()
    return


def show_s_res(x, ys, u_a):
    for i in range(len(ys)):
        plt.plot(x[i].tolist(), ys[i].tolist(), label=f'u_p, step={2*2**i+1}')
    plt.xlabel('x')
    plt.plot(x[-1], u_a, label=f'u_a')
    plt.ylabel('U(x, 1)')
    plt.legend()
    plt.title(f'Последний слои рассчетного решения ')
    plt.show()
    return


def show_log(x, y):
    plt.plot(x, y, marker='o')
    plt.xlabel('log(h)')
    plt.ylabel('log(diff)')
    plt.grid()
    plt.title('Зависимость log(diff) от log(h)')
    plt.show()
    return


def analytics(x, t):
    X, T = np.meshgrid(x, t)
    u = u_analytics(T, X)
    return X, T, u


def count_s(x_size, tau, h, prev_layer, s_prev, first_elem, last_elem):
    k1 = np.zeros((x_size))
    k2 = np.zeros((x_size))
    for l in range(1, x_size - 1):
        k1[l] = (s_prev[l+1] ** NU + s_prev[l] ** NU) / (2 * h ** 2)
        k2[l] = (s_prev[l] ** NU + s_prev[l-1] ** NU) / (2 * h ** 2)
    a = -k1
    b = 1/tau + k1 + k2
    c = -k2
    d = prev_layer/tau

    alpha = np.zeros((x_size))
    beta = np.zeros((x_size))
    alpha[0] = 0
    beta[0] = first_elem
    # alpha[-1] = 0
    # beta[-1] = last_elem
    for i in range(1, x_size):
        alpha[i] = -a[i] / (b[i] + c[i] * alpha[i - 1])
        beta[i] = (d[i] - c[i] * beta[i - 1]) / (b[i] + c[i] * alpha[i - 1])

    s = np.zeros((x_size))
    s[-1] = last_elem
    for i in range(x_size - 2, -1, -1):
        s[i] = alpha[i] * s[i + 1] + beta[i]
    return s


def progonka(x, t, h, tau, u_a):
    t_size = t.size
    x_size = x.size
    u = np.zeros((x_size, t_size))

    # граничные условия
    for l in range(x_size):
        u[l, 0] = (C1 + x[l]) ** (2/NU) * C0 ** (-1/NU)

    for n in range(t_size):
        u[0, n] = C1 ** (2 / NU) * (C0 - 2 * (2 + NU) * t[n] / NU) ** (-1/NU)
        u[-1, n] = (1 + C1) ** (2 / NU) * (C0 - 2 * (2 + NU) * t[n] / NU) ** (-1/NU)

    u = u.transpose()

    for now_layer in range(1, t_size):
        s_res = u[now_layer - 1] # в начале одно и то же
        s_res_arr = []
        s_res_arr.append(s_res)
        s_res_new = count_s(x_size, tau, h, u[now_layer - 1], s_res, u[now_layer, 0], u[now_layer, -1])
        s_res_arr.append(s_res_new)
        s_diff = abs(s_res_new - s_res)
        s_diff_otn = s_diff / abs(s_res)
        ITER = 0
        s_res = s_res_new

        while max(s_diff_otn) > EPS:
            ITER += 1
            s_res_new = count_s(x_size, tau, h, u[now_layer - 1], s_res, u[now_layer, 0], u[now_layer, -1])
            s_diff = abs(s_res_new - s_res)
            s_diff_otn = s_diff / s_res
            s_res = s_res_new
            s_res_arr.append(s_res)
        u[now_layer] = s_res
        # show_s_res(x, s_res_arr, now_layer, u_a[now_layer])
    return u


u_last_layer_p = []
x_s = []
log_diff = []
log_h = []

for i in range(5):
    h_kol = 2 * 2 ** i + 1
    t_kol = h_kol + 2
    h = 1/(h_kol - 1)
    tau = 1/(t_kol - 1)

    x = np.linspace(0, 1, h_kol)
    t = np.linspace(0, 1, t_kol)
    X_a, T_a, u_a = analytics(x, t)

    u_p = progonka(x, t, h, tau, u_a)

    show_graph(X_a, T_a, u_a, u_p, h, tau)

    x_s.append(x)
    u_last_layer_p.append(u_p[-1])
    log_diff.append(np.log(np.max(abs(u_a - u_p))))
    log_h.append(np.log(h))

show_s_res(x_s, u_last_layer_p, u_a[-1])
show_log(log_diff,log_h)
