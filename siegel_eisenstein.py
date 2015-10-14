# -*- coding: utf-8; mode: sage -*-
from sage.all import ZZ, prime_factors, zeta, quadratic_L_function__exact, QQ, MatrixSpace
from siegel_series.impl import siegel_series_polynomial, X
import operator


class SiegelEisensteinSeries(object):

    '''
    A class for Siegel Eisenstein series of general degrees.
    '''

    def __init__(self, weight=None, degree=None):
        self._weight = weight
        self._degree = degree
        self._fc_dct = {}

    @property
    def weight(self):
        return ZZ(self._weight)

    @property
    def degree(self):
        return self._degree

    def _fourier_coefficient_full_rank(self, mat):
        '''
        Assuming mat is a positive definite, half-integral matrix,
        it returns the Fourier coefficient of self at the matrix.
        We normalize the Siegel Eisenstein series so that the constant term is equal to 1.
        cf. S. Takemori, Siegel Eisenstein series of degree n and Lambda-adic Eisenstein series.
        '''
        n = self.degree
        k = self.weight
        if n % 2 == 0:
            _mat_det = ZZ(mat.det() * ZZ(2) ** n)
        else:
            _mat_det = ZZ(mat.det() * ZZ(2) ** (n - 1))
        unramfac = reduce(operator.mul,
                          (siegel_series_polynomial(mat, p).subs(
                              {X: p ** (k - n - 1)})
                           for p in prime_factors(_mat_det)), ZZ(1))
        return self._product_of_l_part(k, _mat_det, n) * unramfac

    def _product_of_l_part(self, k, mat_det, n):
        '''Factor of Fourier coefficient except unramified factors.
        '''
        if n % 2 == 0:
            quadl = quadratic_L_function__exact(
                1 + n // 2 - k, (-1) ** (n // 2) * mat_det)
            prod = reduce(operator.mul,
                          (zeta(1 + 2 * i - 2 * k)
                           for i in range(1, n // 2 + 1)), ZZ(1)) ** (-1)
            return ZZ(2) ** (n // 2) * zeta(1 - k) ** (-1) * quadl * prod
        else:
            prod = reduce(operator.mul,
                          (zeta(1 + 2 * i - 2 * k)
                           for i in range(1, (n - 1) // 2 + 1)), ZZ(1)) ** (-1)
            return ZZ(2) ** ((n + 1) // 2) * zeta(1 - k) ** (-1) * prod

    def fourier_coefficient(self, mat):
        '''Assuming mat is a half-integral, semi positive definite matrix,
        returns the Fourier coefficient at mat.
        '''
        if mat == 0:
            return QQ(1)
        n = mat.ncols()
        MZ = MatrixSpace(ZZ, n)
        _, u, _ = MZ(mat * 2, n).smith_form()
        r = mat.rank()
        if r == n:
            return self._fourier_coefficient_full_rank(mat)
        else:
            mat_rk_r = (u * mat * u.transpose()).submatrix(ncols=r, nrows=r)
            eisen_deg_r = SiegelEisensteinSeries(weight=self.weight, degree=r)
            return eisen_deg_r._fourier_coefficient_full_rank(mat_rk_r)
