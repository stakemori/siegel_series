from sage.all import ZZ, matrix, diagonal_matrix, eisenstein_series_qexp
from ..pullback_of_siegel_eisen import eisenstein_pullback_coeff
from degree2.all import eisenstein_series_degree2
from degree2.basic_operation import PrecisionDeg2
import unittest


def tpl_to_half_int_mat(t):
    n, r, m = t
    return matrix([[ZZ(n), ZZ(r) / ZZ(2)], [ZZ(r) / ZZ(2), ZZ(m)]])


class PullBackTest(unittest.TestCase):

    def assert_pullback_m_1(self, k, A1, prec=10):
        f0 = eisenstein_pullback_coeff(k, A1, diagonal_matrix([0]))
        f = eisenstein_series_qexp(k, prec=prec, normalization='constant') * f0
        for a in range(prec):
            self.assertEqual(f[a],
                             eisenstein_pullback_coeff(k, A1, diagonal_matrix([a])))

    def assert_pullback_m_2(self, k, A1, prec=3):
        f0 = eisenstein_pullback_coeff(k, A1, tpl_to_half_int_mat((0, 0, 0)))
        f = eisenstein_series_degree2(k, prec=5) * f0
        for t in PrecisionDeg2(prec):
            self.assertEqual(f[t],
                             eisenstein_pullback_coeff(k, A1, tpl_to_half_int_mat(t)))

    def test_pullback_2_1(self):
        self.assert_pullback_m_1(6, tpl_to_half_int_mat((4, 0, 6)))
        self.assert_pullback_m_1(6, tpl_to_half_int_mat((1, 1, 1)))

    def test_pullback_3_1(self):
        self.assert_pullback_m_1(6, diagonal_matrix([1, 1, 1]))
        self.assert_pullback_m_1(6, matrix([[1, 1, 0],
                                            [1, 1, 0],
                                            [0, 0, 2]]))

    def test_pullback_4_1(self):
        self.assert_pullback_m_1(8, diagonal_matrix([1, 1, 1, 1]), prec=10)
        self.assert_pullback_m_1(8, diagonal_matrix([2, 1, 1, 1]), prec=8)

    def test_pullback_2_2(self):
        self.assert_pullback_m_2(6, tpl_to_half_int_mat((1, 0, 1)))
        self.assert_pullback_m_2(6, tpl_to_half_int_mat((1, 1, 1)))


suite = unittest.TestLoader().loadTestsFromTestCase(PullBackTest)
unittest.TextTestRunner(verbosity=2).run(suite)
