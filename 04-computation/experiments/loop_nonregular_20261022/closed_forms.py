"""Closed forms of leaf classes in the five-matrix basis.  Session opus-2026-10-22-S54.

Every propagator f_d(Shat) of the loop engine has simple poles at the fourth roots of unity, so it is a combination of
  I,  G = -Shat^2 (I + Shat^2)^-1,  H = Shat^2 (I - Shat^2)^-1,  O = Shat (I - Shat^2)^-1,  Q = Shat (I + Shat^2)^-1.
For a core, expand each term W * prod_e f_(d_e)(sigma_e s_e) (sigma_e = -1 on pendant edges, whose leaf index is
contracted: (1^T f(Shat))_w = (f(-Shat) 1)_w) in this basis edge by edge, and collect the coefficient tensor
C[b_1, ..., b_E] / |G|. The class value is sum_b C[b] * contraction(basis matrices on the core graph).
Prints the nonzero coefficients for each core of the requested classes and checks the reassembled value against the
engine.  Usage: python closed_forms.py "e,k" ["e,k" ...]"""
import sys, os, itertools
from collections import defaultdict
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import leaf_engine as LE
from cluster_brute_classes import almost_regular

s = sp.symbols('s')
BASIS = ['I', 'G', 'H', 'O', 'Q']
BFUN = [sp.Integer(1), -s ** 2 / (1 + s ** 2), s ** 2 / (1 - s ** 2), s / (1 - s ** 2), s / (1 + s ** 2)]


def fsym(desc, x):
    R = lambda sg: sg * x / (1 - sg * x)
    if desc[0] == 'm1':
        return R(1) + R(-1) if desc[1] == '*' else R(desc[1])
    A = lambda sg: R(1) + R(-1) if sg == '*' else R(sg)
    return -A(desc[1]) * A(desc[2]) * (1 - x ** 2) / (1 + x ** 2)


def to_basis(expr):
    """coefficients c_b with expr = sum c_b BFUN[b] (rational function of s with poles in {+-1, +-i} and bounded at oo)."""
    cs = sp.symbols('c0:5')
    eq = sp.together(expr - sum(c * f for c, f in zip(cs, BFUN)))
    num = sp.Poly(sp.numer(eq), s)
    sol = sp.solve(num.coeffs(), cs, dict=True)
    assert sol, expr
    return [sp.nsimplify(sol[0].get(c, 0)) for c in cs]


_cache = {}


def basis_coeffs(desc, sigma):
    key = (desc, sigma)
    if key not in _cache:
        _cache[key] = to_basis(sp.together(fsym(desc, sigma * s)))
    return _cache[key]


def core_tensor(core, k):
    E = len(core.edges)
    C = defaultdict(lambda: sp.Integer(0))
    it = np.nditer(core.W, flags=['multi_index'])
    for w in it:
        wv = int(round(float(w)))
        if wv == 0:
            continue
        descs = [LE.LEF.DESCS[core.used[e][i]] for e, i in enumerate(it.multi_index)]
        coeffs = [basis_coeffs(d, -1 if core.edges[e][0] < k else 1) for e, d in enumerate(descs)]
        for b in itertools.product(range(5), repeat=E):
            c = wv
            for e in range(E):
                c = c * coeffs[e][b[e]]
                if c == 0:
                    break
            if c != 0:
                C[b] += c
    return {b: sp.nsimplify(v / core.G) for b, v in C.items() if v != 0}


def matrices(Sh):
    n = Sh.shape[0]; I = np.eye(n); S2 = Sh @ Sh
    return [I, -S2 @ np.linalg.inv(I + S2), S2 @ np.linalg.inv(I - S2), Sh @ np.linalg.inv(I - S2), Sh @ np.linalg.inv(I + S2)]


def contract(core, b, M, k, n):
    letters = 'abcdefgh'
    terms, ops = [], []
    for e, (a, c) in enumerate(core.edges):
        X = M[b[e]]
        if a == c:
            terms.append(letters[a]); ops.append(np.diag(X))
        elif a < k:                                      # pendant edge: the basis was taken at -s (transpose), so the
            terms.append(letters[c] + letters[a]); ops.append(X)   # vector is (B 1)_w = sum_leaf B[w, leaf]
        else:
            terms.append(letters[a] + letters[c]); ops.append(X)
    spec = ','.join(terms) + '->'
    return float(np.einsum(spec, *ops, optimize='greedy'))


if __name__ == '__main__':
    rng = np.random.default_rng(9)
    A = almost_regular(12, rng)
    Sh = (A - A.T).astype(float) / 12
    M = matrices(Sh)
    for arg in sys.argv[1:]:
        e, k = map(int, arg.split(','))
        total = 0.0
        for core in LE.get(e, k):
            C = core_tensor(core, k)
            val = sum(float(c) * contract(core, b, M, k, 12) for b, c in C.items())
            total += val
            print(f'class ({e},{k}) core edges {core.edges} (leaves 0..{k - 1}), |G| = {core.G}:')
            for b, c in sorted(C.items(), key=lambda t: str(t[0])):
                print(f'    {c} * [' + ', '.join(f'{BASIS[i]}{"(leaf)" if core.edges[j][0] < k else ""}' for j, i in enumerate(b)) + ']')
        eng = LE.X(Sh, e, k).real
        print(f'class ({e},{k}): reassembled {total:+.15e}, engine {eng:+.15e}, diff {total - eng:+.1e}', flush=True)
