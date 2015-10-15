from sage.all import MatrixSpace, ZZ


def non_deg_submatrix(mat):
    '''mat is a degenrate n by n half-integral symmetric matrix.
    r is rank of mat.
    Return u in GL_n(Z) such that submatrix of first r rows and r columns of
    u * mat * u.transpose() is definite.
    '''
    n = mat.ncols()
    MZ = MatrixSpace(ZZ, n)
    _, u, _ = MZ(mat * 2, n).smith_form()
    return u


def is_semi_positive_definite(mat):
    '''Return true only when every eigenvalue of mat is non-negative.
    '''
    r = mat.rank()
    if r == mat.ncols():
        return mat.is_positive_definite()
    else:
        u = non_deg_submatrix(mat)
        mat1 = (u * mat * u.transpose()).submatrix(ncols=r, nrows=r)
        return mat1.is_positive_definite()
