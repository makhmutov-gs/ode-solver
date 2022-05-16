import numpy as np
from matplotlib import pyplot as plt


def f(x, y):
    res = np.zeros(y.size)
    res[0] = y[1]
    res[1] = y[2]
    res[2] = y[3]
    res[3] = y[4]
    res[4] = -(15 * y[4] + 90 * y[3] + 270 * y[2] + 405 * y[1] + 243 * y[0])

    return res


def ode_solver(f, xspan, y0):
    h = (xspan[1] - xspan[0]) / 1000

    x = np.linspace(xspan[0], xspan[1], 1001)
    y = np.zeros((1001, y0.size))

    y[0, :] = y0

    for i in range(1, 1001):
        k1 = f(x[i - 1], y[i - 1, :])
        k2 = f(x[i - 1] + h / 2, y[i - 1, :] + h / 2 * k1)
        k3 = f(x[i - 1] + h / 2, y[i - 1, :] + h / 2 * k2)
        k4 = f(x[i - 1] + h, y[i - 1, :] + h * k3)

        y[i, :] = y[i - 1, :] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    return x, y


if __name__ == '__main__':
    xspan = np.array([0, 5], float)
    y0 = np.array([0, 3, -9, -8, 0], float)

    x, y = ode_solver(f, xspan, y0)
    ref = -1 / 12 * x * np.exp(-3 * x) * (-36 - 54 * x + 16 * x**2 + 129 * x**3)
    cpp_results = np.loadtxt('cpp_results.txt')

    fig, ax1 = plt.subplots()
    ax1.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    ax1.plot(x, y[:, 0], label='Решение из Python', color='tab:orange')
    plt.legend()
    plt.savefig('python_solution.png', dpi=200)

    fig, ax2 = plt.subplots()
    ax2.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    ax2.plot(cpp_results[:, 0], cpp_results[:, 1], label='Решение из C++', color='tab:green')
    plt.legend()
    plt.savefig('cpp_solution.png', dpi=200)

    fig, ax3 = plt.subplots()
    ax3.grid(True)
    plt.xlabel('x')
    plt.ylabel('y')
    ax3.plot(x, ref, label='Аналитическое решение')
    ax3.plot(x, y[:, 0], marker='o', markevery=10, linewidth=0, label='Решение из Python')
    ax3.plot(cpp_results[:, 0], cpp_results[:, 1], marker='+', markevery=50, linewidth=0, label='Решение из C++')
    plt.legend()
    plt.savefig('total.png', dpi=200)

