# -*- coding: utf-8; mode: sage -*-
import unittest
from sage.all import (matrix, QQ, ZZ, prime_factors, mul, QuadraticForm,
                      diagonal_matrix, block_diagonal_matrix,
                      least_quadratic_nonresidue)
from degree2.all import eisenstein_series_degree2
from siegel_series.impl import siegel_series_polynomial, X
from siegel_series.local_invariants import zeta_p, e_p
from siegel_series.tests.utils import random_even_symm_mat
from siegel_series.pull_back_of_siegel_eisen import r_n_m_iter

h_mat = matrix([[QQ(2), QQ(1) / QQ(2)],
                [QQ(1) / QQ(2), QQ(1)]])

y_mat = matrix([[QQ(1), QQ(1) / QQ(2)],
                [QQ(1) / QQ(2), QQ(1)]])


def mat_to_tuple(m):
    return tuple(int(a) for a in ((m[(0, 0)]), 2 * m[(0, 1)], m[(1, 1)]))


class SiegelSeriesTest(unittest.TestCase):

    def assert_deg2_fc_eql(self, k, mat):
        n, r, m = mat_to_tuple(mat)
        content = max(n, r, m)
        es = eisenstein_series_degree2(k, prec=max(n, r, m))
        det_4 = ZZ(mat.det() * 4)
        fc_unram_fac = mul(
            siegel_series_polynomial(mat, p).subs({X: p ** (k - 3)})
            for p in prime_factors(det_4))
        self.assertEqual(es._fc__unramfactor(content, det_4),
                         fc_unram_fac)

    def test_degree2_fc(self):
        self.assert_deg2_fc_eql(4, h_mat)
        self.assert_deg2_fc_eql(4, y_mat)
        self.assert_deg2_fc_eql(4, ZZ(2) * y_mat)
        self.assert_deg2_fc_eql(4, ZZ(2) * h_mat)

    def assert_func_eq(self, B, p):
        p = ZZ(p)
        n = B.ncols()
        q = QuadraticForm(ZZ, ZZ(2) * B)
        pl = siegel_series_polynomial(B, p)
        eb = e_p(q, p)
        self.assertEqual(pl.subs({X: p ** (- n - 1) * X ** (-1)}) *
                         p ** (((n + 1) * eb) // 2) * X ** (eb),
                         zeta_p(q, p) * pl)

    def test_degree1_func_eq(self):
        self.assert_func_eq(diagonal_matrix([ZZ(5) ** 10]), 5)

    def test_degree2_func_eq(self):
        self.assert_func_eq(ZZ(2) * y_mat, 2)
        self.assert_func_eq(ZZ(2) * h_mat, 2)
        self.assert_func_eq(ZZ(2 ** 10) * y_mat, 2)
        u = diagonal_matrix([1, 2])
        self.assert_func_eq(ZZ(3) * u, 3)
        self.assert_func_eq(ZZ(3 ** 2) * u, 3)

    def test_degree3_func_eq_2(self):
        self.assert_func_eq(block_diagonal_matrix(
            ZZ(2 ** 3) * diagonal_matrix([1]), h_mat), 2)
        self.assert_func_eq(block_diagonal_matrix(
            diagonal_matrix([1]), ZZ(2) ** 10 * h_mat), 2)
        self.assert_func_eq(block_diagonal_matrix(
            diagonal_matrix([1]), h_mat) * ZZ(2) ** 10, 2)
        self.assert_func_eq(block_diagonal_matrix(
            ZZ(2 ** 3) * diagonal_matrix([1]), h_mat) * ZZ(2) ** 10, 2)
        self.assert_func_eq(block_diagonal_matrix(
            diagonal_matrix([1]), ZZ(2 ** 3) * h_mat) * ZZ(2) ** 10, 2)

    def test_degree3_func_eq_odd(self):
        p = ZZ(23)
        u = least_quadratic_nonresidue(p)
        mat = diagonal_matrix([1, p ** 2 * u, p ** 3])
        mat1 = diagonal_matrix([1, u, u]) * p ** 10
        self.assert_func_eq(mat, p)
        self.assert_func_eq(mat1, p)

    def test_degree4_func_eq_2(self):
        self.assert_func_eq(
            block_diagonal_matrix(
                diagonal_matrix([1]) * ZZ(2) ** 2,
                diagonal_matrix([ZZ(3)]),
                h_mat), 2)
        self.assert_func_eq(
            block_diagonal_matrix(
                h_mat, h_mat) * ZZ(2) ** 3, 2)
        self.assert_func_eq(
            block_diagonal_matrix(
                ZZ(2) * y_mat,
                diagonal_matrix([ZZ(1), ZZ(3)])), 2)
        self.assert_func_eq(
            diagonal_matrix([ZZ(1), ZZ(3), ZZ(5), ZZ(7) * ZZ(2) ** 3])
            * ZZ(2) ** 10, 2)
        A1 = matrix(2, [ZZ(1), ZZ(0), ZZ(0), ZZ(2)])
        A2 = matrix(2, [ZZ(3), ZZ(0), ZZ(0), ZZ(3)])
        for _, mat in r_n_m_iter(A1, A2):
            if mat.rank() == 4:
                self.assert_func_eq(mat, 2)

    def test_degree4_func_eq_odd(self):
        p = ZZ(7)
        u = least_quadratic_nonresidue(p)
        self.assert_func_eq(
            diagonal_matrix([ZZ(1), p * u, p ** 2, p ** 3]) * p ** 10,
            p)
        self.assert_func_eq(
            diagonal_matrix([u, p * u, p ** 2, p ** 3]), p)

    def test_func_eq_general_degree_random(self):
        for _ in range(30):
            for n in [3, 4, 5, 6, 7]:
                m = random_even_symm_mat(n)
                for p in [ZZ(2), ZZ(3), ZZ(5), ZZ(7)]:
                    self.assert_func_eq(m, p)


suite = unittest.TestLoader().loadTestsFromTestCase(SiegelSeriesTest)
unittest.TextTestRunner(verbosity=2).run(suite)
