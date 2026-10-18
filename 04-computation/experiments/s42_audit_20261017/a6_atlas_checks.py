#!/usr/bin/env python3
"""S42 audit a6: independent checks of the elementary claims L1-L4, L6 and of the 'tests run' section of
preprint_connections_solfeggio_digits_20261017.md (own code; no session code imported).
"""
import os, sys, json, math, itertools, time
import numpy as np
from numba import njit
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from a2_independent_incremental import parity_table, fmap

out = {}


# ---------------- L1 ----------------
@njit(cache=False)
def l1_scan(Nmax, par):
    H4 = 0; H1 = 0; Rprev = 0
    lip = 0; lip3 = 0; ident_fail = 0
    maxratio = 0.0
    K = 40
    D = np.zeros(K); mx = np.zeros(K)
    for N in range(1, Nmax):
        for k in range(4 * N - 3, 4 * N + 1):
            H4 += (2 if k % 3 == 2 else -1) * (1 - 2 * par[k])
        H1 += (2 if N % 3 == 2 else -1) * (1 - 2 * par[N])
        R = H4 - H1
        if N > 1:
            d = R - Rprev
            if abs(d) > lip:
                lip = abs(d)
            # identity: R(N) - R(N-1) = sum_{k=4N-3}^{4N-1} w(k) f(k)   (the 4N term cancels the N term)
            s3 = 0
            for k in range(4 * N - 3, 4 * N):
                s3 += (2 if k % 3 == 2 else -1) * (1 - 2 * par[k])
            if s3 != d:
                ident_fail += 1
        Rprev = R
        r = abs(R) / N
        if r > maxratio:
            maxratio = r
        kk = 0
        x = N
        while x > 1:
            x >>= 1
            kk += 1
        D[kk] += R * R
        if abs(R) > mx[kk]:
            mx[kk] = abs(R)
    return lip, ident_fail, maxratio, D, mx


t0 = time.time()
KMAX = 25
par = parity_table(4 * 2 ** (KMAX + 1))
lip, ident_fail, maxratio, D, mx = l1_scan(2 ** (KMAX + 1), par)
rows = []
for k in range(6, KMAX + 1):
    b = 4 * D[k] ** (1 / 3)
    rows.append(dict(k=k, maxR=float(mx[k]), bound=float(b), ratio=float(b / mx[k]), lemma_holds=bool(mx[k] <= b)))
out['L1'] = dict(max_abs_dR=int(lip), identity_dR_equals_three_term_sum_failures=int(ident_fail),
                 max_abs_R_over_N=float(maxratio), octaves=rows,
                 side_condition_arithmetic='A <= 9 N0 and N0 < 2X give only A/16 < 18X/16 = 1.125 X; '
                                           'the trivial bound |R(N)| <= 4N + 2 gives A/16 < (8X + 2)/16 < X for X >= 1')
print('L1: max |R(N+1)-R(N)| over N < 2^26 = %d; identity Delta R = w f summed over 4N+1..4N+3 fails at %d N; '
      'max |R(N)|/N = %.4f; ratio bound/max at 2^6 = %.2f, 2^25 = %.2f; lemma holds in all octaves: %s [%.0f s]'
      % (lip, ident_fail, maxratio, rows[0]['ratio'], rows[-1]['ratio'], all(r['lemma_holds'] for r in rows), time.time() - t0))


# ---------------- L2 ----------------
@njit(cache=False)
def l2_check(M):
    bad_sign = 0; bad_prod = 0
    for m in range(1, M, 2):
        v = 3 * m + 1
        a = 0
        while v % 2 == 0:
            v //= 2
            a += 1
        leg = 1 if v % 3 == 1 else (-1 if v % 3 == 2 else 0)
        if leg != (1 - 2 * (a & 1)):
            bad_sign += 1
    for n in range(3, M, 2):
        x = n; steps = 0; prod = 1
        while x != 1:
            # one odd-to-odd move with shortcut steps
            v = 3 * x + 1
            a = 0
            while v % 2 == 0:
                v //= 2
                a += 1
            steps += a          # T: (3x+1)/2 counts one step, then a-1 halvings
            r = v % 3
            prod *= (1 if r == 1 else (-1 if r == 2 else 0))
            x = v
        if prod != 1 - 2 * (steps & 1):
            bad_prod += 1
    return bad_sign, bad_prod


bs, bp = l2_check(10 ** 6)
out['L2'] = dict(odd_m_lt_1e6_sign_identity_failures=int(bs), odd_n_lt_1e6_product_identity_failures=int(bp))
print('L2: (-1)^a = (u/3) fails for %d odd m < 10^6; (-1)^sigma_T(n) = prod (u/3) fails for %d odd n < 10^6' % (bs, bp))


# ---------------- L3 ----------------
z = sp.symbols('z')
l3_gf = sp.simplify(z / (1 - z) - (z + z ** 2) / (1 - z ** 2)) == 0
seeds_a = {r: (lambda v: (v & -v).bit_length() - 1)(3 * r + 1) for r in (1, 3, 7)}


def cword(n, L):
    w = []
    for _ in range(L):
        w.append(n & 1)
        n = n // 2 if n % 2 == 0 else 3 * n + 1
    return tuple(w)


fib = [0, 1]
for i in range(80):
    fib.append(fib[-1] + fib[-2])
cnt_ok = {}
for L in range(1, 17):
    words = set(cword(n, L) for n in range(1, 2 ** (L + 2)))
    golden = sum(1 for w in itertools.product((0, 1), repeat=L) if all(not (w[i] and w[i + 1]) for i in range(L - 1)))
    cnt_ok[L] = (len(words), golden, fib[L + 2])
Q = sp.Matrix([[1, 1], [1, 0]])
Q12 = (Q ** 12).applyfunc(lambda x: x % 9)
pis = next(k for k in range(1, 200) if fib[k] % 9 == 0 and fib[k + 1] % 9 == 1)
# Q-orbits on primitive pairs mod 9 and norm classes x^2 - xy - y^2
prim = [(x, y) for x in range(9) for y in range(9) if (x % 3, y % 3) != (0, 0)]
seen = set(); orbits = []
for p0 in prim:
    if p0 in seen:
        continue
    orb = []; p = p0
    while p not in orb:
        orb.append(p); p = ((p[0] + p[1]) % 9, p[0])
    seen |= set(orb)
    norms = sorted(set((x * x - x * y - y * y) % 9 for x, y in orb))
    orbits.append((len(orb), norms))
out['L3'] = dict(gf_identity=bool(l3_gf), seed_a_values={str(k): v for k, v in seeds_a.items()},
                 parity_word_counts={str(k): v for k, v in cnt_ok.items()},
                 counts_are_fibonacci=all(a == b == c for a, b, c in cnt_ok.values()),
                 Q12_mod9=[[int(x) for x in row] for row in Q12.tolist()], pisano_9=pis, Q_orbits_mod9=orbits)
print('L3: z/(1-z) = (z+z^2)/(1-z^2): %s; a(s) for s = 1,3,7 mod 8: %s; #distinct C-parity words of length L (L<=16) = '
      '#golden words = F_(L+2): %s; Q^12 mod 9 = %s; Pisano(9) = %d; Q-orbits on primitive pairs mod 9 (size, norms): %s'
      % (l3_gf, seeds_a, out['L3']['counts_are_fibonacci'], Q12.tolist(), pis, orbits))


# ---------------- L4 ----------------
x = sp.symbols('x')
B = lambda t: 4 * t + 1
B3 = sp.expand(B(B(B(x))))
sol = [396, 417, 528, 639, 741, 852, 963, 174, 285]
digs = {s: [int(c) % 9 for c in str(s)] for s in sol}
fibre_ok = all(sorted(d) == sorted({d[0] % 3 + 3 * i for i in range(3)}) and
               all((d[i] - d[i + 1]) % 9 == 3 for i in range(2)) for d in digs.values())
div3 = [s for s in sol if s % 3 == 0]; div9 = [s for s in sol if s % 9 == 0]
tab = [(k, pow(2, k, 9), 'fibre147' if pow(2, k, 9) % 3 == 1 else 'fibre258',
        {1: '+-1', 8: '+-1', 2: '+-2', 7: '+-2', 4: '+-4', 5: '+-4'}[pow(2, k, 9)]) for k in range(6)]
units = [1, 2, 4, 5, 7, 8]
edges = list(itertools.combinations(units, 2))
tri = [e for e in edges if e[0] % 3 == e[1] % 3]
bip = [e for e in edges if e[0] % 3 != e[1] % 3]
fdr = [(fib[n] - 1) % 9 + 1 if fib[n] else 9 for n in range(1, 25)]
out['L4'] = dict(B3=str(B3), B3_mod9='x + %d' % (21 % 9) if 64 % 9 == 1 else 'n/a', solfeggio_fibres_descending_by_3=fibre_ok,
                 divisible_by_3=div3, divisible_by_9=div9, two_power_table=tab, K6_edges=len(edges), triangle_edges=len(tri),
                 bipartite_edges=len(bip), fib_digital_roots=fdr,
                 fib_half_complement=all((fib[n] + fib[n + 12]) % 9 == 0 for n in range(1, 60)),
                 pos_369=[n for n in range(1, 25) if fib[n] % 3 == 0])
print('L4: B^3(x) = %s = x + %d mod 9; Solfeggio digit triples are descending-by-3 fibres: %s; div by 3: %d of 9, by 9: %s; '
      '2^k table %s; K6 edges %d = triangles %d + K33 %d; Fibonacci digital roots %s; F_n + F_(n+12) = 0 mod 9: %s; '
      'positions of 3,6,9: %s' % (B3, 21 % 9, fibre_ok, len(div3), div9, tab, len(edges), len(tri), len(bip), fdr,
                                  out['L4']['fib_half_complement'], out['L4']['pos_369']))


# ---------------- L6 ----------------
crt = sorted((3 * a + 2 * b) % 6 for a in range(2) for b in range(3)) == list(range(6))
out['L6'] = dict(crt_bijection=crt)
print('L6: (c2, c3) -> 3c2 + 2c3 mod 6 is a bijection: %s' % crt)


# ---------------- THM-4667 statement 8 / Solfeggio branching types in the SEED TRANSPORT ----------------
def P(s):
    v = 3 * s + 1
    while v % 2 == 0:
        v //= 2
    while v % 8 == 5:
        v = (v - 1) // 4
    return v


LIM = 3 * 10 ** 6
pre_min = {}
for s in range(1, LIM, 2):
    if s % 8 == 5 or s == 1:
        continue
    u = P(s)
    if u not in pre_min or s < pre_min[u]:
        pre_min[u] = s
summary = {0: [0, 0, 0, None], 1: [0, 0, 0, None], 2: [0, 0, 0, None]}
violations = 0
for u in range(3, 20001, 2):
    if u % 8 == 5:
        continue
    r = u % 3
    mn = pre_min.get(u)
    summary[r][0] += 1
    if mn is not None:
        summary[r][1] += 1
        if mn < u:
            summary[r][2] += 1
        ratio = mn / u
        summary[r][3] = ratio if summary[r][3] is None else min(summary[r][3], ratio)
    if r == 2 and mn != (2 * u - 1) // 3:
        violations += 1
out['statement8'] = dict(seeds_checked_3_to_20000=sum(v[0] for v in summary.values()),
                         by_residue={str(r): dict(seeds=v[0], with_preimage_below_3e6=v[1], with_preimage_below_u=v[2],
                                                  min_ratio_minpre_over_u=v[3]) for r, v in summary.items()},
                         u2mod3_min_preimage_is_2u_minus_1_over_3_failures=violations,
                         example_preimages_of_3=[s for s in range(1, 200, 2) if s % 8 != 5 and s != 1 and P(s) == 3])
print('statement 8 (seeds 3..20000): by u mod 3 -> [#seeds, #with a preimage < 3e6, #with a preimage < u, min(min preimage/u)]: %s; '
      'u = 2 mod 3 whose least preimage is not (2u-1)/3: %d; preimages of 3 below 200: %s'
      % (summary, violations, out['statement8']['example_preimages_of_3']))


# ---------------- 'tests run' item 4: shares and Legendre correlation ----------------
from a1_frontier_bruteforce import brute, fpar
shares = {}
for N in (7910, 66463, 1048576, 6904411):
    for Y in (100, 1000):
        c, E, R, rdir, bad = brute(N, Y)
        us = np.array([u for u in range(1, Y + 1) if u % 2 == 1 and u % 8 != 5])
        cu = c[us].astype(float)
        Ef = float(cu @ cu)
        chi = np.array([0 if u % 3 == 0 else (1 if u % 3 == 1 else -1) for u in us], dtype=float)
        rho_session = float(cu @ chi) / math.sqrt(len(us) * Ef)
        rho_cos = float(cu @ chi) / math.sqrt(Ef * (chi @ chi))
        sh = {name: float((cu[sel] ** 2).sum() / Ef) for name, sel in
              (('0', us % 3 == 0), ('1', us % 3 == 1), ('2', us % 3 == 2),
               ('g1', np.isin(us % 9, (1, 8))), ('g2', np.isin(us % 9, (2, 7))), ('g4', np.isin(us % 9, (4, 5))))}
        shares['%d_%d' % (N, Y)] = dict(rho_session_normalisation=rho_session, cosine_with_chi=rho_cos, shares=sh)
        print('N = %d Y = %d: rho (session normalisation sqrt(F_Y E)) = %+.4f, true cosine with (u/3) = %+.4f; shares %s'
              % (N, Y, rho_session, rho_cos, {k: round(v, 3) for k, v in sh.items()}))
out['shares'] = shares
json.dump(out, open(os.path.join(HERE, 'a6_atlas_checks.json'), 'w'), indent=1, default=str)
