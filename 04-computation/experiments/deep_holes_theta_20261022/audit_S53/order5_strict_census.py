"""Audit check (S53): strict Fricke points of order 5 and their unit multiples (HYP-9368 (A)).

HYP-9368 defines a strict Fricke point of order n by: Lambda + v has a unique vector of norm 2 - 2/n (and none
smaller), and every proper multiple jv (j | n, 1 < j < n) is deep or in Lambda. For prime n there are no proper
multiples, so nothing constrains the unit multiples 2v, 3v, ... . Rationality of the twining (needed for
Gamma_0(n)-invariance, THM-4760 (0')) forces theta(Lambda + 2v) = theta(Lambda + v).

Translate so that the unique shallow vector is v itself: v = y/5 with y in Lambda, y^2 = 40, and v strictly inside
the Voronoi cell of 0 (no lambda != 0 with |v - lambda|^2 <= 8/5). We sample y of norm 40 near random directions,
keep the strict ones, and record the norm profile of Lambda + 2v below and at norm 2. A strict point whose 2v-profile
is not {8/5: 1, 2: 24} has irrational twining (THM-4760 (2), necessary condition), so it is NOT equivalent to the
rational (h+1)-point c''(A3^8); it is then a counterexample to Conjecture (A) and to the sentence 'THM-4759 (10)
proves that then Z_v - 24 is the Fricke Hauptmodul'.
Usage: python order5_strict_census.py [N_DIRECTIONS=3000] [seed=1]"""
import sys, time, math
import numpy as np
from audit_common import GLAM_F, GLAM_I, shell_vecs, coset_profile

ND = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 1
P = 5
NORM = 2 * P * (P - 1)          # |y|^2 = p^2 (2 - 2/p)
rng = np.random.default_rng(SEED)
L = np.linalg.cholesky(GLAM_F)   # G = L L^T
Linv = np.linalg.inv(L)
DEN = 10 ** 6
t0 = time.time()
seen = set()
n_y = n_strict = 0
tally = {}
examples = {}
for it in range(ND):
    # random direction on the sphere of radius sqrt(NORM): x = L^T c has |x| = sqrt(NORM)
    g = rng.standard_normal(24)
    g *= math.sqrt(NORM) / np.linalg.norm(g)
    c = Linv.T @ g                                    # coefficient vector with c^T G c = NORM
    Yd = np.round(-c * DEN).astype(np.int64)
    vals, vecs = shell_vecs(GLAM_F, GLAM_I, Yd, DEN, 3.0)
    for a in vecs:
        a = np.asarray(a, dtype=np.int64)
        if int(a @ GLAM_I @ a) != NORM:
            continue
        key = tuple(a)
        if key in seen or tuple(-a) in seen:
            continue
        seen.add(key)
        n_y += 1
        # strictness: all lambda with |lambda - y/5|^2 <= 8/5; only lambda = 0 allowed
        pr1 = coset_profile(-a, P, 2.0)               # profile of Lambda - y/5 (same as Lambda + y/5)
        fr = P * P * 2 - 2 * P                         # 25 * 8/5 = 40
        if min(pr1) != fr or pr1[fr] != 1:
            continue
        n_strict += 1
        pr2 = coset_profile(-2 * a, P, 2.0)
        key2 = tuple(sorted(pr2.items()))
        key1 = pr1.get(2 * P * P, 0)
        k = (key1, key2)
        tally[k] = tally.get(k, 0) + 1
        examples.setdefault(k, a.copy())
print(f'order {P}: sampled {n_y} distinct y of norm {NORM} (up to sign) from {ND} directions; strict Fricke points: {n_strict} [{time.time() - t0:.0f}s]')
for (n2, k2), cnt in sorted(tally.items(), key=lambda kv: -kv[1]):
    prof2 = {round(s / 25, 4): c for s, c in k2}
    rational_ok = (prof2 == {1.6: 1, 2.0: 24}) and n2 == 24
    print(f'  x{cnt:5d}: #norm-2 in Lambda+v = {n2:3d}; profile of Lambda+2v = {prof2}   '
          f'{"(consistent with rational, c\'\'-like)" if rational_ok else "(2v NOT shallow-unique => twining irrational, not c\'\'(A3^8))"}')
for k, a in examples.items():
    print('  example y (LAM coords):', list(map(int, a)), ' class', k[0], dict((round(s / 25, 4), c) for s, c in k[1]))
