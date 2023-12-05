# Created by Hannah at 04.12.2023 13:19
import math
import numpy as np
import matplotlib.pyplot as plt
from init_conditions_hannah import *
# from init_conditions_nast import *
# from init_conditions_alin import *
# from init_conditions_1 import *
# from init_conditions_2 import *
# from  init_conditions_3 import *
# from init_conditions_5 import *


def k_const(x):
    return k(0.5)


def q_const(x):
    return q(0.5)


def f_const(x):
    return f(0.5)


def draw(x, y_s, who):
    if who == "dyn":
        for i in range(len(y_s)):
            plt.plot(x[i], y_s[i], label='y'+str(10*2**i+1))
            plt.xlabel('x')
            plt.ylabel('y')
    elif who == 'const':
        plt.plot(x[0], y_s[0], label='analytic')
        plt.plot(x[0], y_s[1], label='const progonka')
        plt.xlabel('x')
        plt.ylabel('y')
    else:
        plt.plot(x, y_s, label = 'log(err)')
        plt.xlabel('log(h)')
        plt.ylabel('log(err)')
    plt.grid(True)
    plt.legend()
    plt.show()


# модельная
def analytics(x):
    c1 = ((k(x_num) * lambda2 + delta2) * (delta1 * f(x_num) - eps1 * q(x_num)) * math.exp(lambda2) + (
            k(x_num) * lambda2 - delta1) * (delta2 * f(x_num) - eps2 * q(x_num))) / \
         (q(x_num) * ((k(x_num) * lambda1 - delta1) * (k(x_num) * lambda2 + delta2) * math.exp(lambda2) - (
                 k(x_num) * lambda2 - delta1) * (k(x_num) * lambda1 + delta2) * math.exp(lambda1)))

    c2 = ((k(x_num) * lambda1 + delta2) * (delta1 * f(x_num) - eps1 * q(x_num)) * math.exp(lambda1) + (
            k(x_num) * lambda1 - delta1) * (delta2 * f(x_num) - eps2 * q(x_num))) / \
         (q(x_num) * ((k(x_num) * lambda2 - delta1) * (k(x_num) * lambda1 + delta2) * math.exp(lambda1) - (
                 k(x_num) * lambda2 + delta2) * (k(x_num) * lambda1 - delta1) * math.exp(lambda2)))

    # c1 = ((k(x_num)* lambda2 + delta2) * (delta1 * f(x_num) - eps1 * q(x_num)) * math.exp(lambda2) + (k(x_num) * lambda2 - delta1) * (
    #         delta2 * f(x_num) - eps2 * q(x_num))) / (
    #               q(x_num) * (k(x_num)* lambda1 - delta1) * (k(x_num)* lambda2 + delta2) * math.exp(lambda2) - (k(x_num)* lambda2 - delta1) * (
    #               k(x_num)* lambda1 + delta2) * math.exp(lambda1))
    # c2 = ((k(x_num)* lambda1 + delta2) * (delta1 * f(x_num) - eps1 * q(x_num)) * math.exp(lambda1) + (k(x_num)* lambda1 - delta1) * (
    #         delta2 * f(x_num) - eps2 * q(x_num))) / (
    #               q(x_num) * (k(x_num)* lambda2 - delta1) * (k(x_num)* lambda1 + delta2) * math.exp(lambda1) - (k(x_num)* lambda2 + delta2) * (
    #                k(x_num)* lambda1 - delta1) * math.exp(lambda2))

    y = []
    # y_a = []
    # print(c1, c2)
    for i in x:
        # y.append(c1 * math.exp(lambda1 * i) + c2 * math.exp(lambda2 * i) + f(i)/q(i))
        #print(i)
        y.append(c1 * np.exp(lambda1 * i) + c2 * np.exp(lambda2 * i) + f(x_num) / q(x_num))
        # y_a.append(c1_a * np.exp(lambda1 * i) + c2_a * np.exp(lambda2 * i) + f(x_num)/q(x_num))

    print('Значения для модельной задачи:\nС1 = {}, C2 = {}'.format(c1, c2))
    return y


# прогонки
def method_progonki(k, q, f, x):
    # массивчики коэффициентов
    a = [0] * num_of_points
    b = [0] * num_of_points
    alpha = [0] * num_of_points
    beta = [0] * num_of_points
    c = [0] * num_of_points
    d = [0] * num_of_points
    u = [0] * num_of_points

    # нулевые значения
    a[0] = k(x[0])
    b[0] = -(k(x[0]) + h * delta1)
    c[0] = 0
    d[0] = -eps1 * h
    alpha[0] = -(a[0] / b[0])
    beta[0] = d[0] / b[0]

    for i in range(1, num_of_points - 1):
        a[i] = k(x[i] + 1 / 2 * h)
        b[i] = -(k(x[i] - 1 / 2 * h) + k(x[i] + 1 / 2 * h) + q(x[i]) * (h ** 2))
        c[i] = k(x[i] - 1 / 2 * h)
        d[i] = -f(x[i]) * (h ** 2)
        alpha[i] = -a[i] / (b[i] + c[i] * alpha[i - 1])
        beta[i] = (d[i] - c[i] * beta[i - 1]) / (b[i] + c[i] * alpha[i - 1])

    a[-1] = 0
    b[-1] = -(k(x[-1]) + h * delta2)
    c[-1] = k(x[-1])
    d[-1] = -eps2 * h
    u[-1] = (d[-1] - c[-1] * beta[-2]) / (b[-1] + c[-1] * alpha[-2])

    # print("a:", a)
    # print("b:", b)
    # print("c:", c)
    # print("d:", d)
    # print("alpha:", alpha)
    # print("beta:", beta)

    for i in range(num_of_points - 2, -1, -1):
        u[i] = alpha[i] * u[i + 1] + beta[i]

    return u


x_start = 0
x_finish = 1
num_of_points = 10+1
h = (x_finish - x_start) / (num_of_points - 1)

x_num = 0.5
lambda1 = math.sqrt(q(x_num) / k(x_num))
lambda2 = -lambda1

x = []
x.append(list(np.linspace(x_start, x_finish, num_of_points)))

y = analytics(x[0])
u_const = method_progonki(k_const, q_const, f_const, x[0])


array = []
array.append(y)
array.append(u_const)

draw(x, array, 'const')


# print()
# print("res model:", y)
# print("res progonka const:", u_const)


u_dynamic = []
delta_e = 1000
need_e = 10**(-4)
h_prev = 1
max_diff = 100
log_e = []
log_h = []
while max_diff > need_e:#delta_e > need_e:
    new_u = method_progonki(k, q, f, x[len(u_dynamic)])
    u_dynamic.append(new_u)
    if len(u_dynamic) > 1:
        max_diff = 0
        for i in range(0, len(u_dynamic[-2])):#(len(u_dynamic[-2]) - 1)//(10)):
            max_diff = max(max_diff, abs(u_dynamic[-2][i] - u_dynamic[-1][2*i]))
        # delta_e = max_diff*h_prev
        # delta_e = math.log(max_diff)
        # print(max_diff, delta_e)
        log_h.append(math.log(h_prev))
        log_e.append(math.log(max_diff))
    num_of_points = (num_of_points - 1) * 2 + 1
    x_new = np.linspace(x_start, x_finish, num_of_points)
    x.append(x_new)
    h_prev = h
    h = (x_finish - x_start) / (num_of_points - 1)

draw(x, u_dynamic, 'dyn')
draw(log_h, log_e, 'err')
print('Order:', (log_e[0] - log_e[1])/(log_h[0] - log_h[1]))
# print(len(u_dynamic), len(u_dynamic[-1]))
# print("res progonka:", u_dynamic[-1][::(len(u_dynamic[-1]) - 1)//(len(x[0])-1)])

