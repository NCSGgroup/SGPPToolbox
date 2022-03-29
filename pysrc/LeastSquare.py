import numpy as np


def get_matrix_A(function, t):
    m = len(t)
    arg_count = function.__code__.co_argcount - 1
    identity_m = np.eye(m)
    identity_m_extension = [list(identity_m[i]) * arg_count for i in range(m)]
    identity_m_extension = np.array(identity_m_extension).reshape(m * arg_count, m)

    t_extension = np.dot(identity_m_extension, t)
    identity_arg = np.eye(arg_count)
    arg_extension = [np.array(list(identity_arg[i]) * m) for i in range(arg_count)]

    A = function(t_extension, *arg_extension).reshape(m, arg_count)
    return A


def dot_for_more(*mat):
    if len(mat) == 2:
        return np.dot(mat[0], mat[1])
    else:
        mat0 = mat[0]
        mat1 = mat[1]
        return dot_for_more(np.dot(mat0, mat1), mat[2:])


def curve_fit(function, t, *y):
    A = get_matrix_A(function, t)
    AT = A.T
    ATA_I = np.linalg.inv(np.dot(AT, A))
    A_ginv = np.dot(ATA_I, AT)
    A_ginv_A_ginv_T = np.dot(A_ginv, A_ginv.T)

    if len(y) == 1:
        y = y[0].reshape(-1, 1)
    else:
        y = np.vstack(y).T

    results = np.dot(A_ginv, y)

    epsilon = np.dot(A, results) - y

    if np.shape(epsilon)[1] == 1:
        var = np.var(epsilon) * A_ginv_A_ginv_T
    else:
        var = np.array([np.var(epsilon[:, i]) * A_ginv_A_ginv_T for i in range(np.shape(epsilon)[1])])

    return results.T, var


def fit_function(x, a, b, c, d):
    return a + b * x + c * np.sin(2 * np.pi * x) + d * np.cos(2 * np.pi * x)


def demo1_compare():
    import scipy.optimize as opt
    t = np.arange(2005, 2021, 1 / 12)
    y_true = fit_function(t, 200, 2, 5, 3)
    y_obs_list = []
    for i in range(10):
        y_obs_list.append(np.random.normal(y_true, 2))

    for i in range(len(y_obs_list)):
        z_opt = opt.curve_fit(fit_function, t, y_obs_list[i])
        z = curve_fit(fit_function, t, y_obs_list[i])

        print('experiment {}:'.format(i + 1))
        print('=' * 50)
        print('opt')
        print(z_opt[0], '\n', z[0], sep='')
        print('sigma')
        print(np.sqrt(np.diag(z_opt[1])), '\n', np.sqrt(np.diag(z[1])), sep='')
        print('=' * 50)


def demo0():
    import scipy.optimize as opt
    import time
    t = np.arange(2005, 2021, 1 / 12)
    y = []

    y_true = fit_function(t, 2000, 2.3, 5.8, 4.5)
    for i in range(10):
        noise = np.random.normal(0, 1, len(t))
        y.append(y_true - y_true[0] + noise)
        # plt.plot(t, y[i])

    time1 = time.time()
    # z = curve_fit(fit_function, t, y[0])
    z = curve_fit(fit_function, t, *y)
    time2 = time.time()
    print('cost {} s'.format(round(time2 - time1, 3)))
    # return 0

    for i in range(len(y)):
        print('=' * 50)
        print(np.sqrt(np.diag(z[1][i])))
        # print(z[0][i])

        z_opt = opt.curve_fit(fit_function, t, y[i])
        print(np.sqrt(np.diag(z_opt[1])))
        print(np.sqrt(np.diag(z_opt[1])) / np.sqrt(np.diag(z[1][i])))
        # print(z_opt[0])


if __name__ == '__main__':
    demo0()
