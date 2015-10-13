import unittest
from sage.all import eisenstein_series_qexp, matrix, ZZ
from siegel_series.siegel_eisenstein import SiegelEisensteinSeries as sess
from degree2.all import eisenstein_series_degree2

class ForierCoeffsLowDegrees(unittest.TestCase):
    def assert_degree_1(self, k):
        es = eisenstein_series_qexp(k, prec=11, normalization="constant")
        es1 = sess(weight=k, degree=1)
        self.assertTrue(all(es[a] == es1.fourier_coefficient(matrix([[a]]))
                            for a in range(11)))

    def test_degree_1(self):
        '''
        Test Fourier coefficients of Siegel Eisenstein series of degree 1.
        '''
        for k in [4, 6, 8, 10]:
            self.assert_degree_1(k)

    def assert_degree_2(self, k):
        es = eisenstein_series_degree2(k, prec=10)
        es1 = sess(weight=k, degree=2)
        self.assertTrue(
            all(es[(n, r, m)] == es1.fourier_coefficient(matrix([[n, ZZ(r)/ZZ(2)],
                                                                 [ZZ(r)/ZZ(2), m]]))
                for n, r, m in es.prec))

    def test_degree_2(self):
        '''
        Test Fourier coefficients of Siegel Eisenstein series of degree 2.
        '''
        for k in [4, 6, 8, 10]:
            self.assert_degree_2(k)


suite = unittest.TestLoader().loadTestsFromTestCase(ForierCoeffsLowDegrees)
unittest.TextTestRunner(verbosity=2).run(suite)
