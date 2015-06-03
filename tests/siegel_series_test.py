# -*- coding: utf-8; mode: sage -*-
import unittest
from sage.all import (matrix, QQ, ZZ, prime_factors, mul,
                      fundamental_discriminant, valuation)
from degree2.all import eisenstein_series_degree2
from siegel_series.impl import siegel_series_polynomial, X

h_mat = matrix([[QQ(2), QQ(1)/QQ(2)],
                [QQ(1)/QQ(2), QQ(1)]])

y_mat = matrix([[QQ(1), QQ(1)/QQ(2)],
                [QQ(1)/QQ(2), QQ(1)]])

def mat_to_tuple(m):
    return tuple(int(a) for a in ((m[(0, 0)]), 2*m[(0, 1)], m[(1, 1)]))


class SiegelSeriesTest(unittest.TestCase):
    def assert_deg2_fc_eql(self, k, mat):
        n, r, m = mat_to_tuple(mat)
        content = max(n, r, m)
        es = eisenstein_series_degree2(k, prec=max(n, r, m))
        det_4 = ZZ(mat.det() * 4)
        fc_unram_fac = mul(
            siegel_series_polynomial(mat, p).subs({X: p**(k - 3)})
            for p in prime_factors(det_4))
        self.assertEqual(es._fc__unramfactor(content, det_4),
                         fc_unram_fac)

    def test_degree2_fc(self):
        self.assert_deg2_fc_eql(4, h_mat)
        self.assert_deg2_fc_eql(4, y_mat)
        self.assert_deg2_fc_eql(4, ZZ(2) * y_mat)
        self.assert_deg2_fc_eql(4, ZZ(2) * h_mat)

    def assert_degree2_func_eq(self, p, mat):
        p = ZZ(p)
        det_4 = mat.det() * 4
        fd = fundamental_discriminant(-det_4)
        f = (valuation(det_4, p) - valuation(fd, p)) / ZZ(2)
        pl = siegel_series_polynomial(mat, p)
        self.assertEqual(pl, p**(3*f) * X**(2*f) *
                         pl.subs({X: (p**3 * X)**(-1)}))

    def test_degree2_func_eq(self):
        self.assert_degree2_func_eq(2, ZZ(2) * y_mat)
        self.assert_degree2_func_eq(2, ZZ(2) * h_mat)
        self.assert_degree2_func_eq(2, ZZ(2**10) * y_mat)


suite = unittest.TestLoader().loadTestsFromTestCase(SiegelSeriesTest)
unittest.TextTestRunner(verbosity=2).run(suite)
