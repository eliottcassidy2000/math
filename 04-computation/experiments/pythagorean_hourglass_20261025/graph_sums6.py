# The nine dense six-vertex Paley graph sums (minimum degree >= 3, even edge count; every other six-vertex graph sum
# reduces to five or fewer free vertices or vanishes at p = 3 mod 4).  S_G(p) = sum over distinct x3..x6 not in {0,1}
# of prod_{(i<j) in G} chi(x_j - x_i), with x1 = 0, x2 = 1 (so T_G = p(p-1) S_G).  Computed in O(p^4) by matrix
# products, then each S_G is tested for exact membership in the span of
#   1, p, p^2, p^3, a_p(f8), p a_p(f8), |b(p)|^2, c6(p), CM(-2) weights 3 and 5, and (2/p) times the polynomials.
# mac-mini-2026-10-25-pyth.
import sys
from fractions import Fraction
import numpy as np
from sympy import primerange, jacobi_symbol

K6 = [(i, j) for i in range(1, 7) for j in range(i + 1, 7)]
def minus(removed):
    rem = {tuple(sorted(e)) for e in removed}
    return [e for e in K6 if e not in rem]
GRAPHS = {
    'K6-e':        minus([(1, 2)]),
    'K6-K3':       minus([(1, 2), (1, 3), (2, 3)]),
    'K6-P4':       minus([(1, 2), (2, 3), (3, 4)]),
    'K6-(P3+K2)':  minus([(1, 2), (2, 3), (4, 5)]),
    'K6-3K2':      minus([(1, 2), (3, 4), (5, 6)]),
    'K6-(C5+K1)':  minus([(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)]),
    'K6-(C4+K2)':  minus([(1, 2), (2, 3), (3, 4), (1, 4), (5, 6)]),
    'K6-(C3+P3)':  minus([(1, 2), (1, 3), (2, 3), (4, 5), (5, 6)]),
    'K6-P6':       minus([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]),
}
for name, E in GRAPHS.items():
    deg = [sum(1 for e in E if v in e) for v in range(1, 7)]
    assert min(deg) >= 3 and len(E) % 2 == 0, (name, deg, len(E))

def S(p, E):
    Es = set(E)
    chi = np.zeros(p); chi[[(x * x) % p for x in range(1, p)]] = 1; chi[chi == 0] = -1; chi[0] = 0
    idx = np.arange(p)
    D = (idx[None, :] - idx[:, None]) % p                  # D[r, c] = c - r
    ChiD = chi[D]                                          # chi(col - row)
    C = ChiD if (5, 6) in Es else (1.0 - np.eye(p))
    tot = 0.0
    for a in range(2, p):                                  # x3 = a
        s = 1.0
        if (1, 3) in Es: s *= chi[a]
        if (2, 3) in Es: s *= chi[(a - 1) % p]
        if s == 0: continue
        w = np.ones(p); w[[0, 1, a]] = 0                   # x4 distinct from 0, 1, a
        if (1, 4) in Es: w *= chi[idx]
        if (2, 4) in Es: w *= chi[(idx - 1) % p]
        if (3, 4) in Es: w *= chi[(idx - a) % p]
        def leaf(v):
            U = np.ones((p, p))                            # rows x4, cols x_v
            m = np.ones(p); m[[0, 1, a]] = 0
            U *= m[None, :]
            U *= (1.0 - np.eye(p))
            if (1, v) in Es: U *= chi[idx][None, :]
            if (2, v) in Es: U *= chi[(idx - 1) % p][None, :]
            if (3, v) in Es: U *= chi[(idx - a) % p][None, :]
            if (4, v) in Es: U *= ChiD                     # chi(x_v - x4)
            return U
        U5, U6 = leaf(5), leaf(6)
        inner = ((U5 @ C) * U6).sum(axis=1)
        tot += s * float(w @ inner)
    return int(round(tot))

if __name__ == '__main__':
    Pmax = int(sys.argv[1])
    ps = [p for p in primerange(7, Pmax) if p % 4 == 3]
    data = {name: {} for name in GRAPHS}
    for p in ps:
        for name, E in GRAPHS.items():
            data[name][p] = S(p, E)
        print(f"p={p} " + " ".join(f"{nm}={data[nm][p]}" for nm in GRAPHS), flush=True)
    # modular data
    sys.argv = [sys.argv[0]]
    from tt6_adjoint_search import f8, im_b, c6, eta_product
    cm3 = eta_product([(1, 2), (2, 1), (4, 1), (8, 2)])    # weight 3 level 8 chi_-8, CM by Q(sqrt-2)
    def cm5(p):
        for x in range(0, int(p ** 0.5) + 1):
            y2 = p - x * x
            if y2 % 2 == 0:
                y = int(round((y2 // 2) ** 0.5))
                if y * y == y2 // 2: return 2 * (x ** 4 - 12 * x * x * y * y + 4 * y ** 4)
        return 0
    cols = {'1': lambda p: 1, 'p': lambda p: p, 'p^2': lambda p: p * p, 'p^3': lambda p: p ** 3,
            'a(f8)': lambda p: f8[p], 'p a(f8)': lambda p: p * f8[p], '|b|^2': lambda p: im_b[p] ** 2,
            'c6': lambda p: c6[p], 'cm3(-2)': lambda p: cm3[p], 'cm5(-2)': cm5,
            '(2/p)': lambda p: int(jacobi_symbol(2, p)), '(2/p) p': lambda p: int(jacobi_symbol(2, p)) * p,
            '(2/p) p^2': lambda p: int(jacobi_symbol(2, p)) * p * p}
    names = list(cols)
    for gname in GRAPHS:
        A = [[Fraction(cols[c](p)) for c in names] + [Fraction(data[gname][p])] for p in ps]
        m, n = len(A), len(names); r = 0; piv = []
        for c in range(n + 1):
            pr = next((i for i in range(r, m) if A[i][c] != 0), None)
            if pr is None: continue
            A[r], A[pr] = A[pr], A[r]; inv = 1 / A[r][c]; A[r] = [x * inv for x in A[r]]
            for i in range(m):
                if i != r and A[i][c] != 0:
                    f = A[i][c]; A[i] = [x - f * y for x, y in zip(A[i], A[r])]
            piv.append(c); r += 1
        if n in piv: print(f"{gname:12s}: NOT in the span ({len(ps)} primes)")
        else:
            sol = {names[c]: A[i][n] for i, c in enumerate(piv)}
            print(f"{gname:12s}: EXACT on {len(ps)} primes ({n} columns): S = " +
                  " + ".join(f"({v})*{k}" for k, v in sol.items() if v != 0))
