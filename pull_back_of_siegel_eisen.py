from sage.all import ZZ, floor, sqrt, matrix, block_matrix, QQ
from siegel_series.siegel_eisenstein import SiegelEisensteinSeries as sess
from siegel_series.utils import is_semi_positive_definite


def _r_iter(n, m):
    sq = int(floor(2 * sqrt(n * m)))
    for r in range(-sq, sq + 1):
        yield r


def r_2_2_iter(A1, A2):
    '''Return a generator of the tuple of 2 by 2 matrice R and
    mat = block_matrix([[A1, R/2], [R^t/2, A2]]) such that
    mat is half-integral, semi positive definite.
    '''
    a1_0, a1_1 = A1[0, 0], A1[1, 1]
    a2_0, a2_1 = A2[0, 0], A2[1, 1]
    for r00 in _r_iter(a1_0, a2_0):
        for r01 in _r_iter(a1_0, a2_1):
            for r10 in _r_iter(a1_1, a2_0):
                for r11 in _r_iter(a1_1, a2_1):
                    R = matrix([[ZZ(r00), ZZ(r01)], [ZZ(r10), ZZ(r11)]])
                    mat = block_matrix(
                        [[A1, R / ZZ(2)], [R.transpose() / ZZ(2), A2]])
                    if is_semi_positive_definite(mat):
                        yield (R, mat)


def r_2_1_iter(A1, A2):
    '''A1 is 2 by 2 half integral matrix, A2 is 1 by 1 matrix.
    Return a generator of the tuple of 2 by 1 matrice R and
    mat = block_matrix([[A1, R/2], [R^t/2, A2]]) such that
    mat is half-integral, semi positive definite.
    '''
    a1_0, a1_1 = A1[0, 0], A1[1, 1]
    a2 = A2[0, 0]
    for r0 in _r_iter(a1_0, a2):
        for r1 in _r_iter(a1_1, a2):
            R = matrix([[r0], [r1]])
            mat = block_matrix([[A1, R / ZZ(2)], [R.transpose() / ZZ(2), A2]])
            if is_semi_positive_definite(mat):
                yield (R, mat)


def eisenstein_pull_back_coeff_low_degree(n, k, A1, A2, func=None):
    '''Return Fourier coefficient of pull back of Siegel Eisenstein series of degree n,
    where n = 3, 4.
    return sum(func(a(block_matrix([[A1, R/2], [R^t/2, A2]]); E_{k, n}), R))
    Here a(T; E_{k, n}) is the Fourier coefficient of Siegel Eisenstein series.
    '''
    if func is None:
        func = lambda x, y: x
    es = sess(weight=k, degree=n)
    if n == 3:
        l = r_2_1_iter(A1, A2)
    elif n == 4:
        l = r_2_2_iter(A1, A2)
    else:
        raise NotImplementedError
    res = QQ(0)
    for R, mat in l:
        res += func(es.fourier_coefficient(mat), R)
    return res
