# -*- coding: utf-8; mode: sage -*-

import unittest
from sage.all import ZZ, QuadraticForm, kronecker_symbol
from siegel_series.tests.utils import random_even_symm_mat
from siegel_series.impl import (jordan_blocks_odd, jordan_blocks_2,
                                _blocks_to_quad_form)


class JordanBlockTest(unittest.TestCase):
    def assert_jordan_blcs(self, p, mat):
        p = ZZ(p)
        q = QuadraticForm(ZZ, 2 * mat)
        if p == 2:
            blcs = jordan_blocks_2(q)
        else:
            blcs = jordan_blocks_odd(q, p)
        q1 = _blocks_to_quad_form(blcs, p)
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
