"""Negative controls: the audit checks must fail on wrong inputs."""
import numpy as np, flint
import a3_class_method as A3
import a2_dense_exact as A2
# (1) class lemma with the matrix of p = 3 but the class map of p = 5
orig = A3.phi_cls
A3.phi_cls = lambda c, m, p: orig(c, m, 5 if p == 3 else p)
try:
    A3.class_lemma_sparse(8, 3); print('control 1 FAILED TO FAIL')
except AssertionError:
    print('control 1 ok: class-lemma check rejects the wrong class map')
A3.phi_cls = orig
# (2) minpoly of a perturbed counting matrix (move one unit of mass in one row)
odd, rows = A2.build_C(7, 3)
r = rows[5]; j1 = next(j for j, x in enumerate(r) if x > 0); j0 = next(j for j, x in enumerate(r) if x == 0)
r[j1] -= 1; r[j0] += 1
mp = flint.fmpz_mat(rows).minpoly()
x = flint.fmpz_poly([0, 1])
print('control 2 ok: perturbed matrix has a different minimal polynomial' if mp != (x - 64) * x ** 6 else 'control 2 FAILED TO FAIL')
# (3) the Tao-chain 'rows depend only on x mod p^(J-n)' check on a chain with a non-geometric kernel
import b1_tao_chain as B1
orig_build = B1.build_M
def bad_build(p, J):
    rows, L = orig_build(p, J)
    rows[1][2], rows[1][4] = rows[1][4], rows[1][2]   # perturb one row
    return rows, L
B1.build_M = bad_build
try:
    B1.analyse(3, 3); print('control 3 FAILED TO FAIL')
except AssertionError:
    print('control 3 ok: Tao-chain checks reject a perturbed kernel')
