# -*- coding: utf-8; mode: sage -*-
'''
Compute Siegel series by using recursive equations given by
Katsurada, `An explicit formula for Siegel series,
American jornal of Mathematics 121 (199),415 - 452.'
'''

from sage.all import (cached_function, QuadraticForm, ZZ, QQ,
                      least_quadratic_nonresidue, legendre_symbol, is_odd,
                      matrix, fundamental_discriminant, valuation,
                      kronecker_symbol)
from degree2.utils import list_group_by
import operator

from siegel_series.rec_coeff import (cbb2_0, cbb2_1, _pol_ring, JordanBlock2,
                                     cbb2_dict)

X = _pol_ring().gens()[0]

@cached_function
def does_2adically_rep_zero(a, b, c):
    '''
    Should a, b, c in [1, 3, 5, 7].
    Returns True if ax^2 + by^2 +cz^2 is 2-adically represents zero
    otherwise False.
    '''
    for x in [0, 1, 4]:
        for y in [0, 1, 4]:
            for z in [0, 1, 4]:
                if (any((w%2 for w in [x, y, z])) and
                    (a * x + b * y + c * z) % 8 == 0):
                    return True
    return False


def _trans_jordan_dec_2(us):
    '''
    us = (u1, ... , un)
    Returns a jordan decomp of diag(us) so that the number of diag elements
    is less that or equal to 2.
    For example,
    x1^2 + x2^2 - x3^2 ~ x1^2 + 2*x1*x2.
    x1^2 + x2^2 + x3^2 ~ 3*x1^2 + 2(x2^2 + x2*x3 + x3^2).
    '''
    hyp2_name = "h"
    y2_name = "y"
    if len(us) <= 2:
        return us
    u1, u2, u3 = us[:3]
    u = u1*u2*u3
    l = _trans_jordan_dec_2(us[3:])
    if does_2adically_rep_zero(u1, u2, u3):
        l.extend([(-u) % 8, hyp2_name])
    else:
        l.extend([(3*u)%8, y2_name])
    l_diag = [a for a in l if a not in [hyp2_name, y2_name]]
    if len(l_diag) <= 2:
        return l
    else:
        l_non_diag = [a for a in l if a in [hyp2_name, y2_name]]
        l1 = _trans_jordan_dec_2(l_diag)
        return l1 + l_non_diag

def jordan_blocks_odd(q, p):
    '''
    p: odd prime,
    q: quadratic form.
    Let q ~ sum p**a_i * b_i * x_i^2 be a jordan decomposition.
    Returns the list of (a_i, b_i) sorted so that a_0 >= a_1 >= ...
    '''
    p = ZZ(p)
    res = []
    u = least_quadratic_nonresidue(p)
    for a, q1 in q.jordan_blocks_by_scale_and_unimodular(p):
        for t in q1.Gram_matrix().diagonal():
            if legendre_symbol(t, p) == 1:
                res.append((a, 1))
            else:
                res.append((a, u))
    return list(reversed(res))


def jordan_blocks_2(q):
    '''
    q: quadratic form.
    Returns a list of tuples (a, b)
    `a' is a integer which represents an exponent.
    `b' corresponds a quadratic form with size <= 2,
    that is, `b' is in [1, 3, 5, 7] or is equal to
    strings 'h' or 'y'.
    If `b' is equal to 'h', then the gram matrix of
    the corresponding quadratic form is equal to
    matrix([[0, 1/2], [1/2, 0]]).
    If `b' is equal to 'y', then the gram matrix of
    the corresponding quadratic form is equal to
    matrix([[1, 1/2], [1/2, 1]]).
    If `b' is in [1, 3, 5, 7], then the corresponding
    gram matrix is equal to matrix([[b]]).
    These tuples are sorted so that assumptions of Theorem 4.1 or Theorem 4.2
    hold.
    '''
    ls = []
    for a, q1 in q.jordan_blocks_by_scale_and_unimodular(2):
        m = q1.Gram_matrix()
        while not m.is_zero():
            m00 = m[0, 0]
            if is_odd(m00):
                ls.append((a, m00 % 8))
                m = m.submatrix(1, 1)
            elif m00 == 0:
                ls.append((a + 1, "h"))
                m = m.submatrix(2, 2)
            else:
                ls.append((a + 1, "y"))
                m = m.submatrix(2, 2)
    # ls is a list of jordan blocks but it may contain diagonal entries of
    # units of length >= 3.
    # We take another jordan blocks by _trans_jordan_dec_2.
    res = [(a, b) for a, b in ls if isinstance(b, str)]
    diag_elts = [(a, b) for a, b in ls if not isinstance(b, str)]
    for a, us_w_idx in list_group_by(diag_elts, lambda x: x[0]):
        l = _trans_jordan_dec_2([_b for _, _b in us_w_idx])
        for b in l:
            if isinstance(b, str):
                res.append((a + 1, b))
            else:
                res.append((a, b))
    # Sort the result so that assumptions of Theorem 4.1 and 4.2 hold.
    non_diags = ("h", "y")
    res1 = []
    for m, blcs in sorted(list_group_by(res, lambda x: x[0]),
                          key=lambda x: -x[0]):
        unit_diags = [(a, qf) for a, qf in blcs if qf not in non_diags]
        non_unit_diags = [(a, qf) for a, qf in blcs if qf in non_diags]
        if len(unit_diags) == 1:
            blcs = unit_diags + non_unit_diags
        else:
            blcs = non_unit_diags + unit_diags
        res1.extend(blcs)
    return res1



def siegel_series_polynomial(B, p):
    '''
    B is a half integral matrix and p is a rational prime.
    It returns the polynomial part of the local Siegel series.
    '''
    p = ZZ(p)
    q = QuadraticForm(ZZ, B * 2)
    if p == 2:
        blcs = jordan_blocks_2(q)
        return _siegel_series_polynomial_2(blcs)
    else:
        blcs = jordan_blocks_odd(q, p)
        return _siegel_series_polynomial_odd(blcs, p)


def _blocks_to_quad_form(blcs, p):
    h = matrix([[QQ(0), QQ(1)/QQ(2)],
                [QQ(1)/QQ(2), QQ(0)]])
    y = matrix([[QQ(1), QQ(1)/QQ(2)],
                [QQ(1)/QQ(2), QQ(1)]])
    mat_dict = {"h": h, "y": y}
    mats_w_expt = [(expt, mat_dict[qf] if qf in ("h", "y") else matrix([[qf]]))
                  for expt, qf in blcs]
    qfs = [QuadraticForm(ZZ, m * ZZ(2) * p**expt) for expt, m in mats_w_expt]
    return reduce(operator.add, qfs)


def swap_ls(l, i, j):
    l[i], l[j] = l[j], l[i]

def _siegel_series_polynomial_2(blcs):
    non_diags = ("h", "y")
    max_expt = blcs[0][0]
    first_qfs_w_idx = [(idx, qf) for idx, (expt, qf) in enumerate(blcs)
                       if expt == max_expt]
    unit_diags_w_idx = [(idx, qf) for idx, qf in first_qfs_w_idx
                        if qf not in non_diags]
    q = _blocks_to_quad_form(blcs, 2)
    # Use known formulas if dim(q) in (1, 2).
    if q.dim() == 1:
        return _siegel_series_dim1(q, ZZ(2))
    elif q.dim() == 2:
        return _siegel_series_dim2(q, ZZ(2))

    if len(unit_diags_w_idx) == 1:
        # Use the recursive equation given in Theorem 4.1.

        blcs_q2 = blcs[1:]
        q2 = _blocks_to_quad_form(blcs_q2, 2)
        pol = _siegel_series_polynomial_2(blcs_q2)
        return (cbb2_1(q, q2, 2) * pol.subs({X: ZZ(2)*X}) +
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
    coeffs = [c11 * c21,
              c11*c20 + c10*c21,
              c10 * c20]
    pol = _siegel_series_polynomial_2(blcs_q2)
    pols = [pol.subs({X: ZZ(4) * X}),
            pol.subs({X: ZZ(2) * X}),
            pol]
    return sum((a*b for a, b in zip(coeffs, pols)))


def _siegel_series_polynomial_odd(blcs, p):
    q = _blocks_to_quad_form(blcs, p)
    # Use known formulas when dim(q) = 1 or = 2.
    if q.dim() == 1:
        return _siegel_series_dim1(q, p)
    elif q.dim() == 2:
        return _siegel_series_dim2(q, p)
    # Else use the recursive equation given in Theorem 4.1
    blcs_q2 = blcs[1:]
    pol = _siegel_series_polynomial_odd(blcs_q2, p)
    q2 = _blocks_to_quad_form(blcs_q2, p)
    return (cbb2_1(q, q2, p) * pol.subs({X: ZZ(2)*X}) +
            cbb2_0(q, q2, p) * pol)



def _siegel_series_dim1(q, p):
    return sum(((p*X)**a for a in range(q.content().valuation(p) + 1)))


def _siegel_series_dim2(q, p):
    det_4 = q.Gram_det() * ZZ(4)
    c = q.content().valuation(p)
    fd = fundamental_discriminant(-det_4)
    f = (valuation(det_4, p) - valuation(fd, p)) / ZZ(2)
    return (__siegel_series_dim2(p, c, f + 1) -
            kronecker_symbol(fd, p) * p * X * __siegel_series_dim2(p, c, f))


def __siegel_series_dim2(p, a, b):
    if b == 0:
        return 0
    a = min(a, b - 1)
    r1 = (1 - (p**2 * X)**(a + 1)) / (1 - p**2 * X)
    rn2 = (p**3 * X**2)**b * p*X - p**b * (p * X)**(2*b - a)
    rd2 = p * X - 1
    return (r1 - rn2/rd2) / (1 - p**3 * X**2)
