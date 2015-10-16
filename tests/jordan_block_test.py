# -*- coding: utf-8; mode: sage -*-

import unittest
from sage.all import ZZ, QuadraticForm, kronecker_symbol, valuation, gcd
from siegel_series.tests.utils import random_even_symm_mat
from siegel_series.impl import (jordan_blocks_odd, jordan_blocks_2,
                                _blocks_to_quad_form)


def _i_func(q):
    '''Return i(B) in Katsurada's paper.
    '''
    m = ZZ(2) * (q.matrix()) ** (-1)
    i = valuation(gcd(m.list()), ZZ(2))
    m = ZZ(2) ** (-i) * m
    if all(m[a, a] % 2 == 0 for a in range(m.ncols())):
        return - i - 1
    else:
        return - i


def q2_q2_blc(blcs):
    max_expt = blcs[0][0]
    non_diags = ("h", "y")
    unit_diags_first_blc = [qf for expt, qf in blcs
                            if expt == max_expt and qf not in non_diags]
    if len(unit_diags_first_blc) == 1:
        blcs_q2 = blcs[1:]
    else:
        if blcs[0][1] in non_diags:
            blcs_q2 = blcs[1:]
        else:
            blcs_q2 = blcs[2:]
    q2 = _blocks_to_quad_form(blcs_q2, ZZ(2))
    return (q2, blcs_q2)


class JordanBlockTest(unittest.TestCase):

    def assert_jordan_blcs(self, p, mat):
        p = ZZ(p)
        q = QuadraticForm(ZZ, 2 * mat)
        if p == 2:
            blcs = jordan_blocks_2(q)
        else:
            blcs = jordan_blocks_odd(q, p)
        q1 = _blocks_to_quad_form(blcs, p)
        det_frac = q.det() / q1.det()
        if p == 2:
            should1 = det_frac % 8
        else:
            should1 = kronecker_symbol(det_frac, p)
        self.assertTrue(det_frac.valuation(p) == 0)
        self.assertEqual(q.content().valuation(p), q1.content().valuation(p))
        self.assertEqual(should1, 1)
        self.assertEqual(kronecker_symbol(q.det() / q1.det(), p), 1)
        self.assertEqual(q.hasse_invariant(p), q1.hasse_invariant(p))

    def test_jordan_blcs(self):
        for _ in range(30):
            m = random_even_symm_mat(5)
            for p in [2, 3, 5, 7]:
                self.assert_jordan_blcs(p, m)

    def test_assumption(self):
        '''Test assumption of Theorem 4.1 and 4.2 in Katsurada's paper.
        '''
        for _ in range(100):
            for n in range(4, 11):
                m = random_even_symm_mat(n)
                q = QuadraticForm(ZZ, 2 * m)
                blcs = jordan_blocks_2(q)
                max_expt = blcs[0][0]
                q2, q2_blc = q2_q2_blc(blcs)
                while len(q2_blc) > 2:
                    self.assertTrue(max_expt >= _i_func(q2) + 1)
                    q2, q2_blc = q2_q2_blc(q2_blc)


suite = unittest.TestLoader().loadTestsFromTestCase(JordanBlockTest)
unittest.TextTestRunner(verbosity=2).run(suite)
