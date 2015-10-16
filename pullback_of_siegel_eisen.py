from sage.all import ZZ, floor, sqrt, matrix, block_matrix, QQ
from siegel_series.siegel_eisenstein import SiegelEisensteinSeries as sess
from siegel_series.utils import is_semi_positive_definite
import itertools


def _r_iter(n, m):
    sq = int(floor(2 * sqrt(n * m)))
    for r in range(-sq, sq + 1):
        yield r


def r_n_m_iter(A1, A2):
    '''
    A1: n by n half-integral, symmetric matrix.
    A2: m by m half-integral, symmetric matrix.
    Return a generator of the tuple of n by m matrice R and
    mat = block_matrix([[A1, R/2], [R^t/2, A2]]) such that
    mat is half-integral, semi positive definite.
    '''
    n = A1.ncols()
    m = A2.ncols()
    r_iters = [_r_iter(a1, a2) for a1, a2 in itertools.product(
        [A1[i, i] for i in range(n)], [A2[i, i] for i in range(m)])]
    for rs in itertools.product(*r_iters):
        R = matrix(n, [ZZ(r) for r in rs])
        mat = block_matrix([[A1, R / ZZ(2)],
                            [R.transpose() / ZZ(2), A2]])
        if is_semi_positive_definite(mat):
            yield (R, mat)


def eisenstein_pullback_coeff(k, A1, A2, func=None):
    '''Return Fourier coefficient of pull back of Siegel Eisenstein series of degree n.
    Here n = A1.ncols() + A2.ncols().
    return sum(func(a(mat; E_{k, n}), A1, A2, R, mat))
    Here a(T; E_{k, n}) is the Fourier coefficient of Siegel Eisenstein series and
    mat = block_matrix([[A1, R/2], [R^t/2, A2]]).
    '''
    n = A1.ncols() + A2.ncols()
    if func is None:
        func = lambda a, a1, a2, r, mat: a
    es = sess(weight=k, degree=n)
    res = QQ(0)
    for R, mat in r_n_m_iter(A1, A2):
        res += func(es.fourier_coefficient(mat), A1, A2, R, mat)
    return res
