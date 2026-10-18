#!/usr/bin/env python3
"""Affine Kac-Moody vertex algebras V^k(g) for g = sl_n, so_{2n} (from matrices), as ModeAlgebras for lie_voa_ss.

Basis: Cartan H_i (diagonal) then root vectors; charges = eigenvalues under the Cartan basis (integer tuples).
Structure constants by solving [X_a, X_b] = sum_c f_ab^c X_c; invariant form (x|y) = tr(xy)/ind (ind = Dynkin
index of the defining representation: 1 for sl_n, 2 for so_N, N >= 5), so long roots have length 2.
Modes x_m: [x_m, y_n] = [x,y]_{m+n} + m k (x|y) delta_{m+n,0}.
positive_ops: zero modes of simple-root vectors and the mode-1 lowest-root vector (generate g^hat_+ with n_+).
"""
from fractions import Fraction
import itertools
import sympy as sp
from lie_voa_ss import ModeAlgebra

def mat(n, entries):
    M = sp.zeros(n, n)
    for (i, j), v in entries.items():
        M[i, j] += v
    return M

def sl_basis(n):
    B, names, cart = [], [], []
    for i in range(n - 1):
        B.append(mat(n, {(i, i): 1, (i + 1, i + 1): -1})); names.append(f"H{i+1}")
    cartan = list(B)
    for i in range(n):
        for j in range(n):
            if i != j:
                B.append(mat(n, {(i, j): 1})); names.append(f"E{i+1}{j+1}")
    simple_pos = [names.index(f"E{i+1}{i+2}") for i in range(n - 1)]
    simple_neg = [names.index(f"E{i+2}{i+1}") for i in range(n - 1)]
    lowest = names.index(f"E{n}1")  # -theta
    return B, names, cartan, Fraction(1), simple_pos, simple_neg, lowest

def so_even_basis(n):
    """so(2n) preserving the antidiagonal form; Cartan diag(t1..tn,-tn..-t1)."""
    N = 2 * n
    pr = lambda i: N - 1 - i
    B, names = [], []
    for i in range(n):
        B.append(mat(N, {(i, i): 1, (pr(i), pr(i)): -1})); names.append(f"H{i+1}")
    cartan = list(B)
    for i in range(n):
        for j in range(n):
            if i < j:
                B.append(mat(N, {(i, j): 1, (pr(j), pr(i)): -1})); names.append(f"e{i+1}-e{j+1}")
                B.append(mat(N, {(j, i): 1, (pr(i), pr(j)): -1})); names.append(f"e{j+1}-e{i+1}")
                B.append(mat(N, {(i, pr(j)): 1, (j, pr(i)): -1})); names.append(f"e{i+1}+e{j+1}")
                B.append(mat(N, {(pr(j), i): 1, (pr(i), j): -1})); names.append(f"-e{i+1}-e{j+1}")
    # simple roots of D_n: e1-e2, ..., e_{n-1}-e_n, e_{n-1}+e_n
    simple_pos = [names.index(f"e{i+1}-e{i+2}") for i in range(n - 1)] + [names.index(f"e{n-1}+e{n}")]
    simple_neg = [names.index(f"e{i+2}-e{i+1}") for i in range(n - 1)] + [names.index(f"-e{n-1}-e{n}")]
    lowest = names.index("-e1-e2")  # -theta, theta = e1+e2
    return B, names, cartan, Fraction(2), simple_pos, simple_neg, lowest

class AffineKM(ModeAlgebra):
    def __init__(self, basis_data, k, P):
        B, names, cartan, ind, sp_, sn_, low = basis_data
        self.names = names
        self.k = Fraction(k)
        dim = len(B)
        n = B[0].shape[0]
        # flatten basis for solving
        Mflat = sp.Matrix([[B[a][i, j] for a in range(dim)] for i in range(n) for j in range(n)])
        self.f = {}
        for a in range(dim):
            for b in range(dim):
                C = B[a] * B[b] - B[b] * B[a]
                if C == sp.zeros(n, n):
                    continue
                rhs = sp.Matrix([C[i, j] for i in range(n) for j in range(n)])
                sol, params = Mflat.gauss_jordan_solve(rhs)
                coeffs = {c: Fraction(int(sp.fraction(sol[c])[0]), int(sp.fraction(sol[c])[1]))
                          for c in range(dim) if sol[c] != 0}
                self.f[(a, b)] = coeffs
        self.form = {}
        for a in range(dim):
            for b in range(dim):
                t = (B[a] * B[b]).trace()
                if t != 0:
                    self.form[(a, b)] = Fraction(int(t)) / ind
        # charges: eigenvalue of ad(H_i) on X_a (root vectors are eigenvectors; Cartan elements 0)
        chg = []
        for a in range(dim):
            ch = []
            for Hh in cartan:
                C = Hh * B[a] - B[a] * Hh
                # find scalar lam with C = lam * B[a]
                lam = 0
                for i in range(n):
                    for j in range(n):
                        if B[a][i, j] != 0:
                            lam = C[i, j] / B[a][i, j]
                            break
                    else:
                        continue
                    break
                ch.append(int(lam))
            chg.append(tuple(ch))
        self.gens = [(names[a], 1, 0, chg[a]) for a in range(dim)]
        self.simple_pos, self.simple_neg, self.lowest = sp_, sn_, low
        super().__init__(P)
        self.fP = {key: [(c, self.F(v)) for c, v in d.items()] for key, d in self.f.items()}
        self.formP = {key: self.F(v) for key, v in self.form.items()}
        self.kP = self.F(self.k)

    def bracket(self, m1, m2):
        a_, ta = -m1[0], m1[1]
        b_, tb = -m2[0], m2[1]
        n = a_ + b_
        terms = [((-n, c), co) for c, co in self.fP.get((ta, tb), [])]
        cen = 0
        if n == 0 and (ta, tb) in self.formP:
            cen = self.kP * self.formP[(ta, tb)] % self.P * (a_ % self.P) % self.P
        return terms, cen

    def positive_ops(self):
        return [(0, t) for t in self.simple_pos] + [(-1, self.lowest)]
