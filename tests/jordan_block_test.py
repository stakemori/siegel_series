# -*- coding: utf-8; mode: sage -*-

import unittest
from sage.all import (ZZ, QuadraticForm, diagonal_matrix, matrix, QQ,
                      block_diagonal_matrix, kronecker_symbol, random_matrix)

from siegel_series.impl import jordan_blocks_odd, jordan_blocks_2

def jordan_blcs_to_quad_form_odd(blcs, p):
    p = ZZ(p)
    return QuadraticForm(ZZ, 2 * diagonal_matrix([p**a * b for a, b in blcs]))

def random_even_symm_mat(n):
    while True:
        m = random_matrix(ZZ, n)
        m = m + m.transpose()
        for a in range(n):
            m[(a, a)] = m[(a, a)] * ZZ(2)
        if not m.is_singular():
            return m


def jordan_blcs_to_quad_form_2(blcs):
    def _to_mat(b):
        if isinstance(b, str):
            if b == "h":
                return matrix([[QQ(0), QQ(1)/QQ(2)],
                               [QQ(1)/QQ(2), QQ(0)]])
            else:
                return matrix([[QQ(1), QQ(1)/QQ(2)],
                               [QQ(1)/QQ(2), QQ(1)]])
        else:
            return matrix([[b]])
    blcs_mat = block_diagonal_matrix([ZZ(2)**a * _to_mat(b) for a, b in blcs])
    return QuadraticForm(ZZ, blcs_mat * 2)


def jordan_blocks_to_quad_form(blcs, p):
    if p == 2:
        return jordan_blcs_to_quad_form_2(blcs)
    else:
        return jordan_blcs_to_quad_form_odd(blcs, p)


class JordanBlockTest(unittest.TestCase):
    def assert_jordan_blcs(self, p, mat):
        p = ZZ(p)
        q = QuadraticForm(ZZ, 2 * mat)
        if p == 2:
            blcs = jordan_blocks_2(q)
        else:
            blcs = jordan_blocks_odd(q, p)
        q1 = jordan_blocks_to_quad_form(blcs, p)
        det_frac = q.det()/q1.det()
        if p == 2:
            should1 = det_frac % 8
        else:
            should1 = kronecker_symbol(det_frac, p)
        self.assertTrue(det_frac.valuation(p) == 0)
        self.assertEqual(q.content().valuation(p), q1.content().valuation(p))
        self.assertEqual(should1, 1)
        self.assertEqual(kronecker_symbol(q.det()/q1.det(), p), 1)
        self.assertEqual(q.hasse_invariant(p), q1.hasse_invariant(p))

    def test_jordan_blcs(self):
        for _ in range(30):
            m = random_even_symm_mat(5)
            for p in [2, 3, 5, 7]:
                self.assert_jordan_blcs(p, m)


suite = unittest.TestLoader().loadTestsFromTestCase(JordanBlockTest)
unittest.TextTestRunner(verbosity=2).run(suite)
