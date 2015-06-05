# -*- coding: utf-8; mode: sage -*-
from sage.all import random_matrix, ZZ

def random_even_symm_mat(n):
    while True:
        m = random_matrix(ZZ, n)
        m = m + m.transpose()
        for a in range(n):
            m[(a, a)] = m[(a, a)] * ZZ(2)
        if not m.is_singular():
            return m
