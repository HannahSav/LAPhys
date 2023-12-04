# Created by Hannah at 04.12.2023 23:01
import numpy as np
import matplotlib.pyplot as plt
import math as m


def method_of_analytics(x):
    k = 1 + m.sin(0.5) ** 2
    q = m.cos(0.5)
    f = m.exp(0.5)
    l_1 = np.sqrt(q / k)
    l_2 = -np.sqrt(q / k)
    delta_1 = 1
    delta_2 = 1
    eps_1 = 1
    eps_2 = 1
    C_1 = ((k * l_2 + delta_2) * (delta_1 * f - eps_1 * q) * m.exp(l_2) + (k * l_2 - delta_1) * (
                delta_2 * f - eps_2 * q)) / (
                      q * (k * l_1 - delta_1) * (k * l_2 + delta_2) * m.exp(l_2) - (k * l_2 - delta_1) * (
                          k * l_1 + delta_2) * m.exp(l_1))
    C_2 = ((k * l_1 + delta_2) * (delta_1 * f - eps_1 * q) * m.exp(l_1) + (k * l_1 - delta_1) * (
                delta_2 * f - eps_2 * q)) / (
                      q * (k * l_2 - delta_1) * (k * l_1 + delta_2) * m.exp(l_1) - (k * l_2 - delta_2) * (
                          k * l_1 + delta_1) * m.exp(l_2))
    u = C_1 * np.exp(l_1 * x) + C_2 * np.exp(l_2 * x) + f / q
    return u


x = np.linspace(0, 1, 1000)

y = method_of_analytics(x)

plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('method_of_analytics(x)')
plt.title('Analytics')
plt.grid()
plt.show()
