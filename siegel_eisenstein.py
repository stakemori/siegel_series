# -*- coding: utf-8; mode: sage -*-
from sage.all import ZZ, prime_factors
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
        return self._weight

    @property
    def degree(self):
        return self._degree

    def _fourier_coefficient_full_rank(self, mat):
        '''
        Assuming mat is a positive definite, half-integral matrix,
        it returns the Fourier coefficient of self at the matrix.
        '''
        n = self.degree
        k = self.weight
        if n % 2 == 0:
            _mat_det = mat.det() * ZZ(2)**n
        else:
            _mat_det = mat.det() * ZZ(2)**(n - 1)
        unramfac = reduce(operator.mul,
                          (siegel_series_polynomial(mat, p).subs(
                              {X: p**(n + 1 - k)})
                           for p in prime_factors(_mat_det)))

