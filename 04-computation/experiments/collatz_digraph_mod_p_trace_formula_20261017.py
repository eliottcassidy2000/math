#!/usr/bin/env python3
"""The Collatz affine digraph modulo a prime p (gcd(p, 6) = 1): vertices F_p, arcs x -> x/2 and x -> (3x+1)/2 (S39).

Trace formula (PROVED, elementary): tr(A_p^k) = sum over the 2^k words w of #Fix(g_w), where g_w(x) = (3^j x + c_w)/2^k
(j = number of odd steps) is the composite affine map; #Fix = 1 unless 3^j = 2^k mod p, in which case g_w is a
translation by c_w/2^k and #Fix = p [c_w = 0 mod p].  Hence

    tr(A_p^k) = 2^k + sum_{w : 3^j = 2^k (mod p)} (p [p | c_w] - 1),

the finite-place shadow of the parity-matched near-returns 2^k ~ 3^j of THM-4659: the closed walks of length k count
the 2-adic rational periodic points x_w = c_w/(2^k - 3^j) (one per word; Lagarias) reduced modulo p, and the words
whose denominator 2^k - 3^j vanishes modulo p contribute p or 0 according as the cycle numerator does.
Also computed: the arborescence counts tau_p = det of a reduced Laplacian (= number of Eulerian circuits, BEST, since
the digraph is 2-in 2-out), the characteristic polynomial on the nonconstant functions, and the block structure by
the cosets of <2,3> in F_p^x (the frequency basis; G_m = F_m x| <2,3> is the group of
collatz_affine_guarded_lifts_20261004.md).
Usage: python collatz_digraph_mod_p_trace_formula_20261017.py
Session opus-2026-10-17-S39.
"""
import itertools
import numpy as np
import sympy


def digraph(p):
    i2 = pow(2, -1, p)
    A = np.zeros((p, p), dtype=object)
    for x in range(p):
        A[x, (x * i2) % p] += 1
        A[x, ((3 * x + 1) * i2) % p] += 1
    return A


def word_data(k):
    """(j, c_w) over all words w in {0,1}^k, applying x -> x/2 (0) or (3x+1)/2 (1) in order: g_w(x) = (3^j x + c_w)/2^k."""
    out = []
    for w in itertools.product((0, 1), repeat=k):
        num_a, c, den = 1, 0, 1               # g(x) = (num_a x + c)/den
        for bit in w:
            if bit == 0:
                den *= 2
            else:
                num_a, c, den = 3 * num_a, 3 * c + den, 2 * den
        assert den == 2 ** k
        out.append((sum(w), c))
    return out


def trace_formula(p, k, words):
    t = 2 ** k
    for j, c in words:
        if (pow(3, j, p) - pow(2, k, p)) % p == 0:
            t += p * (1 if c % p == 0 else 0) - 1
    return t


if __name__ == '__main__':
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
        A = digraph(p)
        Ak = np.identity(p, dtype=object)
        for k in range(1, 15):
            Ak = Ak.dot(A)
            tr = int(np.trace(Ak))
            assert tr == trace_formula(p, k, word_data(k)), (p, k)
        L = sympy.Matrix(2 * np.identity(p, dtype=object) - A)
        tau = L[1:, 1:].det()
        z = sympy.symbols('z')
        cp = sympy.Matrix(A).charpoly(z).as_expr()
        q, r = sympy.div(cp, z - 2, z)
        assert r == 0
        H = len({pow(2, i, p) * pow(3, j, p) % p for i in range(p) for j in range(p)})
        fac = sympy.factor_list(q, z)[1]
        degs = sorted(sympy.degree(f, z) for f, m in fac for _ in range(m))
        print('p = %2d: trace formula holds for k <= 14; |<2,3>| = %2d; tau_p = %s = %s; nonconstant char. poly factor degrees %s'
              % (p, H, tau, sympy.factorint(tau), degs))
        assert q.subs(z, 2) == p * tau
    # equidistribution of the cycle numerators over the degenerate words (Perron-Frobenius: |tr A^k - 2^k| <= (p-1) rho^k)
    print('\ncycle numerators c_w over words with 3^j = 2^k (mod p): fraction divisible by p vs 1/p, and the bound')
    for p in [5, 7, 11, 13]:
        A = digraph(p).astype(float)
        ev = sorted(np.abs(np.linalg.eigvals(A)), reverse=True)
        assert abs(ev[0] - 2) < 1e-9
        rho = ev[1]
        line = []
        for k in range(8, 21, 4):
            words = word_data(k)
            deg = [(j, c) for j, c in words if (pow(3, j, p) - pow(2, k, p)) % p == 0]
            B = len(deg); Nk = sum(1 for j, c in deg if c % p == 0)
            if B:
                dev = abs(Nk - B / p)
                assert dev <= (p - 1) * rho ** k / p + 1e-9
                line.append('k=%d: %d/%d = %.4f (bound on |N - B/p| %.1f)' % (k, Nk, B, Nk / B, (p - 1) * rho ** k / p))
        print('p = %2d, rho_p = %.4f: %s' % (p, rho, '; '.join(line)))
    print('ALL CHECKS PASSED')
