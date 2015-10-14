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
