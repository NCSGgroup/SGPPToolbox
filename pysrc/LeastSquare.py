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


def fit_function(x, a, b, c, d, e, f):
    return a + b * x + c * np.sin(2 * np.pi * x) + d * np.cos(2 * np.pi * x) + e * np.sin(4 * np.pi * x) + f * np.cos(
        4 * np.pi * x)
