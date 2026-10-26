# Audit check B: injective graph sums by brute force with EXACT injectivity (every pair of vertices distinct), written
# independently of graph_sums6.py / paley_graph_sums.py.
#   T_G(q) = sum over injective x : V -> F_q of prod_{(i -> j)} chi(x_j - x_i);   S_G = sum with x_a = 0, x_b = 1.
# (i) unnormalised T_G by full enumeration at small q, to confirm T_G = q(q-1) S_G (even |E|);
# (ii) S_W4 (hub 0, non-edges {1,2},{3,4}, order orientation; x0 = 0, x1 = 1) at all odd primes 5..199, compared with
#      the K5-2K2 column of paley_graph_sums.out;
# (iii) the nine dense even six-vertex sums (labels 1..6, order orientation, x1 = 0, x2 = 1) at primes 7..59
#      (q = 3 mod 4) and 5..29 (q = 1 mod 4), compared with graph_sums6.out.   Audit 2026-10-25.
import itertools, re, sys
import numpy as np
from sympy import primerange

def chi_table(p):
    c = np.full(p, -1, dtype=np.int64); c[0] = 0
    c[(np.arange(1, p) ** 2) % p] = 1
    return c

def T_full(p, nv, E):
    """unnormalised injective sum, full enumeration (small p only)"""
    c = chi_table(p); tot = 0
    for x in itertools.permutations(range(p), nv):
        v = 1
        for i, j in E:
            v *= c[(x[j] - x[i]) % p]
            if v == 0: break
        tot += v
    return int(tot)

def S_norm(p, nv, E, a=0, b=1):
    """sum over injective x with x_a = 0, x_b = 1; loops over the first nv-4 free vertices in Python, vectorises the
    last two over a p x p grid; every pair constrained distinct"""
    c = chi_table(p)
    free = [v for v in range(nv) if v not in (a, b)]
    head, (u, w) = free[:-2], free[-2:]
    X, Y = np.meshgrid(np.arange(p), np.arange(p), indexing='ij')   # x_u = X, x_w = Y
    tot = 0
    for vals in itertools.permutations([z for z in range(p) if z not in (0, 1)], len(head)):
        x = {a: 0, b: 1}; x.update(dict(zip(head, vals)))
        fixed = set(x.values())
        mask = (X != Y)
        for z in fixed:
            mask &= (X != z) & (Y != z)
        prod = np.ones((p, p), dtype=np.int64)
        s = 1
        for i, j in E:
            xi = X if i == u else (Y if i == w else x[i])
            xj = X if j == u else (Y if j == w else x[j])
            if isinstance(xi, int) and isinstance(xj, int):
                s *= int(c[(xj - xi) % p])
            else:
                prod = prod * c[(np.asarray(xj) - np.asarray(xi)) % p]
        if s == 0: continue
        tot += s * int(prod[mask].sum())
    return tot

W4 = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)]       # hub 0; non-edges {1,2},{3,4}
K6 = [(i, j) for i in range(1, 7) for j in range(i + 1, 7)]
def minus(rem):
    R = {tuple(sorted(e)) for e in rem}
    return [(i - 1, j - 1) for (i, j) in K6 if (i, j) not in R]           # 0-based internally, order orientation
SIX = {
    'K6-e': minus([(1, 2)]), 'K6-K3': minus([(1, 2), (1, 3), (2, 3)]), 'K6-P4': minus([(1, 2), (2, 3), (3, 4)]),
    'K6-(P3+K2)': minus([(1, 2), (2, 3), (4, 5)]), 'K6-3K2': minus([(1, 2), (3, 4), (5, 6)]),
    'K6-(C5+K1)': minus([(1, 2), (2, 3), (3, 4), (4, 5), (1, 5)]), 'K6-(C4+K2)': minus([(1, 2), (2, 3), (3, 4), (1, 4), (5, 6)]),
    'K6-(C3+P3)': minus([(1, 2), (1, 3), (2, 3), (4, 5), (5, 6)]), 'K6-P6': minus([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]),
}

if __name__ == '__main__':
    out = []
    # (i) normalisation
    for p in (7, 11):
        for name, (nv, E) in {'W4': (5, W4), 'K6-P6': (6, SIX['K6-P6']), 'K6-3K2': (6, SIX['K6-3K2'])}.items():
            T = T_full(p, nv, E); S = S_norm(p, nv, E)
            out.append(f'(i) q={p} {name}: T (full enumeration) = {T}, q(q-1) S = {p * (p - 1) * S}: {"OK" if T == p * (p - 1) * S else "FAIL"}')
    T19 = T_full(19, 5, W4); S19 = S_norm(19, 5, W4)
    out.append(f'(i) q=19 W4: T = {T19}, q(q-1)S = {19 * 18 * S19}: {"OK" if T19 == 342 * S19 else "FAIL"}')

    # (ii) S_W4 vs the lane
    lane = {}
    for l in open('../paley_graph_sums.out'):
        m = re.match(r'p=(\d+) mod4=\d K4=-?\d+ K5=-?\d+ K5-e=-?\d+ K5-2K2=(-?\d+)', l)
        if m: lane[int(m.group(1))] = int(m.group(2))
    mine_w4 = {p: S_norm(p, 5, W4) for p in primerange(5, 200)}
    bad = [p for p in mine_w4 if p in lane and lane[p] != mine_w4[p]]
    out.append(f'(ii) S_W4 at {len(mine_w4)} odd primes 5..199: disagreements with paley_graph_sums.out: {bad}')
    gp = {}
    for l in open('modular_coeffs_gp.out'):
        m = re.match(r'p=(\d+) a8=(-?\d+) bre=(-?\d+) bim=(-?\d+)', l)
        if m: gp[int(m.group(1))] = (int(m.group(2)), int(m.group(3)), int(m.group(4)))
    badw = [p for p in mine_w4 if mine_w4[p] != -4 - gp[p][0]]
    out.append(f'(ii) S_W4 = -4 - a_p(8.4.a.a, PARI) at all {len(mine_w4)} primes 5..199 (both classes): failures {badw}')

    # (iii) six-vertex sums
    lane6 = {}
    for l in open('../graph_sums6.out'):
        if l.startswith('p='):
            p = int(l.split()[0][2:])
            lane6[p] = {kv.split('=')[0]: int(kv.split('=')[1]) for kv in l.split()[1:]}
    res6 = {}
    for p in [7, 11, 19, 23, 31, 43, 47, 59] + [5, 13, 17, 29]:
        res6[p] = {name: S_norm(p, 6, E) for name, E in SIX.items()}
        agree = (p not in lane6) or all(lane6[p][nm] == res6[p][nm] for nm in SIX)
        out.append(f'(iii) q={p} (mod 4 = {p % 4}): ' + ' '.join(f'{nm}={v}' for nm, v in res6[p].items())
                   + ('' if p not in lane6 else f'   lane: {"AGREE" if agree else "DISAGREE"}'))
    # predicted values: S_(K6-P6) = |b|^2 - q^2 - 7 at q = 3 mod 4 (HYP-9276), S_oct = (q-5) S_W4
    for p in res6:
        a8, bre, bim = gp[p]
        b2 = bre * bre + bim * bim
        if p % 4 == 3:
            out.append(f'      q={p}: S_(K6-P6) = |b|^2 - q^2 - 7 ? {res6[p]["K6-P6"] == b2 - p * p - 7};  S_oct = (q-5) S_W4 ? {res6[p]["K6-3K2"] == (p - 5) * mine_w4[p]}')
        else:
            out.append(f'      q={p}: S_(K6-P6) = |b|^2 - q^2 + 21 (HYP-9276, q = 1 mod 4) ? {res6[p]["K6-P6"] == b2 - p * p + 21};  S_oct = (q-5) S_W4 ? {res6[p]["K6-3K2"] == (p - 5) * mine_w4[p]}')
    print('\n'.join(out))
