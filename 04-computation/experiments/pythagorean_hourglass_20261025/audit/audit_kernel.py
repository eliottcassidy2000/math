# Audit check F: kernel-partition conversion between all-maps sums U and injective sums T, by brute force.
#   U(W4)  = T_W4 - q(q-1)(q-4)
#   U(oct) = T_oct + 3(q-4) T_C4 + 3(q)_4 + (q)_3,  T_C4 = -q(q-1)(q-3)
#   HYP-9275: U(oct) = q(q-1)((5-q) a_q + 3q);  the all-maps wheel U(W4) = -q(q-1)(q + a_q)
#   Ahlgren-Ono moment: sum_b chi(b) t(b)^2 = -a_p - p, t(b) = sum_c chi(c(c-1)(c-b)), and U(W4) = q(q-1) * moment
#   THM-4780 (6) Legendre form: S_W4 = sum_(x != 0,1) chi(x) a_x^2 + (q-3), a_x = -t(x)
# a_p from PARI (modular_coeffs_gp.out).  Labellings as in the documents (W4: hub 0, non-edges {1,2},{3,4};
# octahedron: antipodal {1,2},{3,4},{5,6}; order orientation).  Audit 2026-10-25.
import itertools, re, sys
import numpy as np
sys.path.insert(0, '.')
from audit_graph_sums import S_norm, W4
from sympy import primerange

gp = {}
for l in open('modular_coeffs_gp.out'):
    m = re.match(r'p=(\d+) a8=(-?\d+) bre=(-?\d+) bim=(-?\d+)', l)
    if m: gp[int(m.group(1))] = int(m.group(2))

def chi_np(p):
    c = np.full(p, -1, dtype=np.int64); c[0] = 0; c[(np.arange(1, p) ** 2) % p] = 1; return c

def U_all(p, nv, E):
    """all maps x: [nv] -> F_p, translation-normalised x_0 = 0 (U = p * that); last two coordinates vectorised"""
    c = chi_np(p); tot = 0
    X, Y = np.meshgrid(np.arange(p), np.arange(p), indexing='ij')
    for rest in itertools.product(range(p), repeat=nv - 3):
        x = [0] + list(rest)
        prod = np.ones((p, p), dtype=np.int64)
        s = 1
        for i, j in E:
            xi = x[i] if i < nv - 2 else (X if i == nv - 2 else Y)
            xj = x[j] if j < nv - 2 else (X if j == nv - 2 else Y)
            if isinstance(xi, int) and isinstance(xj, int): s *= int(c[(xj - xi) % p])
            else: prod = prod * c[(np.asarray(xj) - np.asarray(xi)) % p]
        if s: tot += s * int(prod.sum())
    return p * tot

OCT = [(i, j) for i in range(6) for j in range(i + 1, 6) if (i, j) not in [(0, 1), (2, 3), (4, 5)]]
out = []
ff = lambda q, k: int(np.prod([q - i for i in range(k)]))
for p in (5, 7, 11, 13, 17):
    a = gp[p]
    Tw = p * (p - 1) * S_norm(p, 5, W4); Uw = U_all(p, 5, W4)
    To = p * (p - 1) * S_norm(p, 6, OCT); Uo = U_all(p, 6, OCT)
    TC4 = -p * (p - 1) * (p - 3)
    out.append(f'q = {p} (mod 4 = {p % 4}), a_q = {a}: U(W4) = {Uw}; T_W4 - q(q-1)(q-4) = {Tw - p * (p - 1) * (p - 4)} '
               f'[{"OK" if Uw == Tw - p * (p - 1) * (p - 4) else "FAIL"}]; -q(q-1)(q+a_q) = {-p * (p - 1) * (p + a)} '
               f'[{"OK" if Uw == -p * (p - 1) * (p + a) else "FAIL"}]')
    rhs = To + 3 * (p - 4) * TC4 + 3 * ff(p, 4) + ff(p, 3)
    out.append(f'      U(oct) = {Uo}; T_oct + 3(q-4)T_C4 + 3(q)_4 + (q)_3 = {rhs} [{"OK" if Uo == rhs else "FAIL"}]; '
               f'q(q-1)((5-q)a_q + 3q) = {p * (p - 1) * ((5 - p) * a + 3 * p)} [{"OK" if Uo == p * (p - 1) * ((5 - p) * a + 3 * p) else "FAIL"}]; '
               f'T_oct = (q-5) T_W4: {To == (p - 5) * Tw}')
# Ahlgren-Ono moment and the Legendre form at all odd primes 3..199
badAO, badL = [], []
for p in primerange(3, 200):
    c = chi_np(p); y = np.arange(p)
    t = np.array([int(c[(y * (y - 1) * (y - b)) % p].sum()) for b in range(p)])
    mom = int((c * t * t).sum())
    if mom != -gp[p] - p: badAO.append(p)
    if p >= 5:
        leg = int(sum(c[x] * t[x] ** 2 for x in range(2, p))) + (p - 3)
        if leg != S_norm(p, 5, W4): badL.append(p)
out.append(f'Ahlgren-Ono moment sum_b chi(b) t(b)^2 = -a_p - p at all odd primes 3..199: failures {badAO}')
out.append(f'Legendre form S_W4 = sum_(x != 0,1) chi(x) a_x^2 + (q-3) at all primes 5..199: failures {badL}')
# U(W4) = q(q-1) * moment, exactly (no collision correction between U and the moment)
for p in (5, 7, 11):
    c = chi_np(p); y = np.arange(p)
    t = np.array([int(c[(y * (y - 1) * (y - b)) % p].sum()) for b in range(p)])
    out.append(f'q = {p}: U(W4) = {U_all(p, 5, W4)}, q(q-1) * moment = {p * (p - 1) * int((c * t * t).sum())}')
print('\n'.join(out))
with open('audit_kernel.out', 'w') as fh: fh.write('\n'.join(out) + '\n')
