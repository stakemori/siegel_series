# -*- coding: utf-8; mode: sage -*-
'''
Compute Siegel series by using recursive equations given by
Katsurada, `An explicit formula for Siegel series,
American jornal of Mathematics 121 (199),415 - 452.'
'''

from sage.all import (QuadraticForm, ZZ, QQ,
                      matrix, fundamental_discriminant, valuation,
                      kronecker_symbol)
import operator

from siegel_series.rec_coeff import (cbb2_0, cbb2_1, _pol_ring, JordanBlock2,
                                     cbb2_dict)
from siegel_series.jordan_block import jordan_blocks_odd, jordan_blocks_2

X = _pol_ring().gens()[0]


def siegel_series_polynomial(B, p):
    '''
    B is a half integral matrix and p is a rational prime.
    It returns the polynomial part of the local Siegel series.
    '''
    p = ZZ(p)
    q = QuadraticForm(ZZ, B * 2)
    if p == 2:
        blcs = jordan_blocks_2(q)
        return _siegel_series_polynomial_2(tuple(blcs))
    else:
        blcs = jordan_blocks_odd(q, p)
        return _siegel_series_polynomial_odd(tuple(blcs), p)


def _blocks_to_quad_form(blcs, p):
    h = matrix([[QQ(0), QQ(1) / QQ(2)],
                [QQ(1) / QQ(2), QQ(0)]])
    y = matrix([[QQ(1), QQ(1) / QQ(2)],
                [QQ(1) / QQ(2), QQ(1)]])
    mat_dict = {"h": h, "y": y}
    mats_w_expt = [(expt, mat_dict[qf] if qf in ("h", "y") else matrix([[qf]]))
                   for expt, qf in blcs]
    qfs = [QuadraticForm(ZZ, m * ZZ(2) * p ** expt) for expt, m in mats_w_expt]
    return reduce(operator.add, qfs)


def swap_ls(l, i, j):
    l[i], l[j] = l[j], l[i]


# @cached_function
def _siegel_series_polynomial_2(blcs):
    non_diags = ("h", "y")
    max_expt = blcs[0][0]
    unit_diags_first_blc = [qf for expt, qf in blcs
                            if expt == max_expt and qf not in non_diags]
    q = _blocks_to_quad_form(blcs, 2)
    # Use known formulas if dim(q) in (1, 2).
    if q.dim() == 1:
        return siegel_series_dim1(q, ZZ(2))
    elif q.dim() == 2:
        return siegel_series_dim2(q, ZZ(2))

    if len(unit_diags_first_blc) == 1:
        # Use the recursive equation given in Theorem 4.1.

        blcs_q2 = blcs[1:]
        q2 = _blocks_to_quad_form(blcs_q2, 2)
        pol = _siegel_series_polynomial_2(tuple(blcs_q2))
        return (cbb2_1(q, q2, 2) * pol.subs({X: ZZ(2) * X}) +
                cbb2_0(q, q2, 2) * pol)

    # Else use the recursive equation given in Theorem 4.2.

    if blcs[0][1] in non_diags:
        m, units_or_hy = blcs[0]
        blcs_q2 = blcs[1:]
    else:
        (m, u1), (_, u2) = blcs[0], blcs[1]
        units_or_hy = [u1, u2]
        blcs_q2 = blcs[2:]

    q2 = _blocks_to_quad_form(blcs_q2, 2)
    b1 = JordanBlock2(units_or_hy, m)
    coeff_dict = cbb2_dict(b1, q2)
    c11 = coeff_dict['11']
    c10 = coeff_dict['10']
    c21 = coeff_dict['21']
    c20 = coeff_dict['20']
    times2 = {X: ZZ(2) * X}
    times4 = {X: ZZ(4) * X}
    coeffs = [c11 * c21.subs(times2),
              c11 * c20.subs(times2) + c10 * c21,
              c10 * c20]
    pol = _siegel_series_polynomial_2(tuple(blcs_q2))
    pols = [pol.subs(times4), pol.subs(times2), pol]
    return sum((a * b for a, b in zip(coeffs, pols)))


# @cached_function
def _siegel_series_polynomial_odd(blcs, p):
    q = _blocks_to_quad_form(blcs, p)
    # Use known formulas when dim(q) = 1 or = 2.
    if q.dim() == 1:
        return siegel_series_dim1(q, p)
    elif q.dim() == 2:
        return siegel_series_dim2(q, p)
    # Else use the recursive equation given in Theorem 4.1
    blcs_q2 = blcs[1:]
    pol = _siegel_series_polynomial_odd(tuple(blcs_q2), p)
    q2 = _blocks_to_quad_form(blcs_q2, p)
    return (cbb2_1(q, q2, p) * pol.subs({X: p * X}) +
            cbb2_0(q, q2, p) * pol)


def siegel_series_dim1(q, p):
    return sum(((p * X) ** a for a in range(q.content().valuation(p) + 1)))


def siegel_series_dim2(q, p):
    det_4 = q.Gram_det() * ZZ(4)
    c = q.content().valuation(p)
    fd = fundamental_discriminant(-det_4)
    f = (valuation(det_4, p) - valuation(fd, p)) / ZZ(2)
    return (_siegel_series_dim2(p, c, f + 1) -
            kronecker_symbol(fd, p) * p * X * _siegel_series_dim2(p, c, f))


def _siegel_series_dim2(p, a, b):
    if b == 0:
        return 0
    a = min(a, b - 1)
    r1 = (1 - (p ** 2 * X) ** (a + 1)) / (1 - p ** 2 * X)
    rn2 = (p ** 3 * X ** 2) ** b * p * X - p ** b * (p * X) ** (2 * b - a)
    rd2 = p * X - 1
    return (r1 - rn2 / rd2) / (1 - p ** 3 * X ** 2)
