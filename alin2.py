# Created by Hannah at 04.12.2023 23:48
import numpy as np
import matplotlib.pyplot as plt
import math as m


def method_of_analytics(x):
    k = m.exp(np.cos(0.5))
    q = m.exp(np.cos(0.5))
    f = 1
    l_1 = np.sqrt(q / k)
    l_2 = -np.sqrt(q / k)
    delta_1 = 1
    delta_2 = -1
    eps_1 = 0
    eps_2 = 0
    C_1 = ((k * l_2 + delta_2) * (delta_1 * f - eps_1 * q) * m.exp(l_2) + (k * l_2 - delta_1) * (
                delta_2 * f - eps_2 * q)) / (
                      q * (k * l_1 - delta_1) * (k * l_2 + delta_2) * m.exp(l_2) - (k * l_2 - delta_1) * (
                          k * l_1 + delta_2) * m.exp(l_1))
    C_2 = ((k * l_1 + delta_2) * (delta_1 * f - eps_1 * q) * m.exp(l_1) + (k * l_1 - delta_1) * (
                delta_2 * f - eps_2 * q)) / (
                      q * (k * l_2 - delta_1) * (k * l_1 + delta_2) * m.exp(l_1) - (k * l_2 + delta_2) * (
                          k * l_1 - delta_1) * m.exp(l_2))
    u = C_1 * np.exp(l_1 * x) + C_2 * np.exp(l_2 * x) + f / q
    print(C_1, C_2)
    return u


x = np.linspace(0, 1, 11)

y = method_of_analytics(x)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('method_of_analytics(x)')
plt.title('Analytics')
plt.grid()
plt.show()

# def method_progonki(x):
#   k = 1 + m.sin(x)**2
#   q = m.cos(x)
#   f = m.exp(x)
#   a_l = k
#   b_l = -2*k - q*(h**2)
#   c_l = k
#   d_l = -f*(h**2)

#   a_0 = k
#   a_L = 0


# прогонки

num_of_points = 11
delta1 = 1
delta2 = -1
eps1 = 0
eps2 = 0
h = 1 / num_of_points


def k(x):
    # return (1 + math.sin(x) ** 2)
    return m.exp(m.cos(x))


def q(x):
    # return math.cos(x)
    return m.exp(m.sin(x))


def f(x):
    # return math.exp(x)
    return 1


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

print("a:", a)
print("b:", b)
print("c:", c)
print("d:", d)
print("alpha:", alpha)
print("beta:", beta)

for i in range(num_of_points - 2, -1, -1):
    u[i] = alpha[i] * u[i + 1] + beta[i]


def draw_model():
    plt.plot(x, y)
    plt.plot(x, u)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()


draw_model()