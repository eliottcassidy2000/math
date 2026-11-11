#!/usr/bin/env python3
"""Verifier checks for HYP-9661 (D): the class law as stated there, and the F1 facts quoted in the ANALOGY.

(1) Exact (rational) test of the first-order class law of a nondegenerate bilinear gluing:
    P(B(u, v) = c) = 1/p + (pi_1 + pi_2 - p^-r)([c = 0] - 1/p) + Rhat(c), with Rhat bilinear in the
    deviations D1 - U, D2 - U; consequences: classes c != 0 equal to first order, excess of class 0 over
    their mean = pi_1 + pi_2 - p^-r + (p/(p-1)) Rhat(0); uniform states: 1/p + p^-r([c = 0] - 1/p) > 0.
(2) The one-apex states: psi recursion, the two integer matrices on (psi, h, floor(h/2), 1).
(3) F1 values N = x1 x2 + 2^s y1 y2 - [s = 0]: self-joins are never 7 (mod 8); among the primes 4q - 1 < 3000
    with q prime, r_F1 = 0 exactly at the eight listed in HYP-9657 (B) (7, 43, 67, 163, ...).
(4) F1 as a source family: s = 0 is D - E + U with E = 1, s >= 1 is D + 2^s U (D = x1 x2, U = y1 y2).
"""
import itertools
import random
import sys
from fractions import Fraction

FAILS = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (" : " + detail if detail else ""))
    if not cond:
        FAILS.append(name)


# ---------------------------------------------------------------- (1) class law
random.seed(9661)


def rand_dist(vecs, zero_mass=None):
    w = [Fraction(random.randint(1, 50)) for _ in vecs]
    tot = sum(w)
    return {v: x / tot for v, x in zip(vecs, w)}


def rand_form(p, r):
    while True:
        M = [[random.randrange(p) for _ in range(r)] for _ in range(r)]
        # determinant mod p (r <= 3)
        if r == 1:
            det = M[0][0]
        elif r == 2:
            det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
        else:
            det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                   - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                   + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        if det % p:
            return M


def law(p, r, M, D1, D2):
    P = [Fraction(0)] * p
    for u, a in D1.items():
        for v, b in D2.items():
            c = sum(u[i] * M[i][j] * v[j] for i in range(r) for j in range(r)) % p
            P[c] += a * b
    return P


n_inst = 0
ok_first = ok_bilinear = ok_excess = ok_uniform = ok_halfuniform = True
for p in (2, 3, 5, 7):
    for r in (1, 2, 3):
        if p ** r > 130:
            continue
        vecs = list(itertools.product(range(p), repeat=r))
        zero = tuple([0] * r)
        U = {v: Fraction(1, p ** r) for v in vecs}
        for _ in range(6):
            M = rand_form(p, r)
            D1, D2 = rand_dist(vecs), rand_dist(vecs)
            P = law(p, r, M, D1, D2)
            pi1, pi2 = D1[zero], D2[zero]
            first = [Fraction(1, p) + (pi1 + pi2 - Fraction(1, p ** r)) * ((1 if c == 0 else 0) - Fraction(1, p))
                     for c in range(p)]
            Rhat = [P[c] - first[c] for c in range(p)]
            # Rhat is the law of the deviations: sum_{u,v} a(u) b(v) [B(u,v) = c], a = D1 - U, b = D2 - U
            a = {v: D1[v] - U[v] for v in vecs}
            b = {v: D2[v] - U[v] for v in vecs}
            dev = law(p, r, M, a, b)
            ok_bilinear &= (dev == Rhat) and sum(Rhat) == 0
            # excess of class 0 over the mean of the others
            if p > 1:
                mean_others = sum(P[1:]) / (p - 1)
                ok_excess &= (P[0] - mean_others == pi1 + pi2 - Fraction(1, p ** r) + Fraction(p, p - 1) * Rhat[0])
            # second-order smallness: scale both deviations by t; remainder scales by t^2
            t = Fraction(1, 7)
            D1t = {v: U[v] + t * a[v] for v in vecs}
            D2t = {v: U[v] + t * b[v] for v in vecs}
            Pt = law(p, r, M, D1t, D2t)
            firstt = [Fraction(1, p) + (D1t[zero] + D2t[zero] - Fraction(1, p ** r))
                      * ((1 if c == 0 else 0) - Fraction(1, p)) for c in range(p)]
            ok_first &= all(Pt[c] - firstt[c] == t * t * Rhat[c] for c in range(p))
            # one state uniform: remainder vanishes
            Ph = law(p, r, M, U, D2)
            ok_halfuniform &= all(Ph[c] == Fraction(1, p) + D2[zero] * ((1 if c == 0 else 0) - Fraction(1, p))
                                  for c in range(p))
            # both uniform
            Pu = law(p, r, M, U, U)
            ok_uniform &= all(Pu[c] == Fraction(1, p) + Fraction(1, p ** r) * ((1 if c == 0 else 0) - Fraction(1, p))
                              and Pu[c] > 0 for c in range(p))
            n_inst += 1
check("class law: remainder = law of the two deviations (bilinear), sums to 0 [%d instances]" % n_inst, ok_bilinear)
check("class law: remainder scales with the product of the two deviations (second order)", ok_first)
check("excess of class 0 over the mean of the others = pi1 + pi2 - p^-r + (p/(p-1)) Rhat(0)", ok_excess)
check("one uniform state: remainder vanishes identically", ok_halfuniform)
check("both uniform: 1/p + p^-r([c = 0] - 1/p), positive in every class", ok_uniform)

# degenerate form: the statement needs nondegeneracy (rank-1 form on (Z/3)^2, uniform states)
p, r = 3, 2
vecs = list(itertools.product(range(p), repeat=r))
Uu = {v: Fraction(1, 9) for v in vecs}
Pd = law(p, r, [[1, 0], [0, 0]], Uu, Uu)
check("degenerate form (rank 1) does not satisfy the uniform law (nondegeneracy is needed)",
      Pd[0] != Fraction(1, 3) + Fraction(1, 9) * (1 - Fraction(1, 3)), str(Pd))


# ---------------------------------------------------------------- (2) one-apex states
def psi(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return psi(h // 2) + h // 2
    k = h // 2
    return psi(k) + k // 2 + 1


M0 = [[1, 1, 0, 0], [0, 2, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]   # digit 0: h -> 2h
M1 = [[1, 0, 1, 1], [0, 2, 0, 1], [0, 1, 0, 0], [0, 0, 0, 1]]   # digit 1: h -> 2h + 1
ok = True
for h in range(0, 4096):
    st = [psi(h), h, h // 2, 1]
    for d, Mx in ((0, M0), (1, M1)):
        new = [sum(Mx[i][j] * st[j] for j in range(4)) for i in range(4)]
        h2 = 2 * h + d
        ok &= new == [psi(h2), h2, h2 // 2, 1]
check("(psi, h, floor(h/2), 1) is mapped by two integer matrices, one per binary digit (h < 4096)", ok)


def xy(h):
    return 2 * psi(h) + 1, 2 * h + 1


check("x <= y, x odd, y odd, y/2 < x (h < 4096)", all(xy(h)[0] % 2 == 1 and xy(h)[1] / 2 < xy(h)[0] <= xy(h)[1]
                                                    for h in range(4096)))

# ---------------------------------------------------------------- (3) F1 values
LIM = 3000
r = {}
st = [xy(h) for h in range(0, LIM)]
selfjoin_mod8 = set()
for s in range(0, 14):
    for (x1, y1) in st:
        if x1 + (2 ** s) * y1 - (s == 0) > LIM:
            break
        for (x2, y2) in st:
            N = x1 * x2 + (2 ** s) * y1 * y2 - (1 if s == 0 else 0)
            if N > LIM:
                break
            r[N] = r.get(N, 0) + 1
for s in range(0, 40):
    for (x1, y1) in st[:2000]:
        selfjoin_mod8.add((x1 * x1 + (2 ** s) * y1 * y1 - (1 if s == 0 else 0)) % 8)
check("a self-join x^2 + 2^s y^2 - [s = 0] is never 7 (mod 8)", 7 not in selfjoin_mod8, str(sorted(selfjoin_mod8)))


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


# note: y is increasing in h and x > y/2, so the breaks above are safe only if N is monotone in h2; recount plainly
r = {}
for s in range(0, 13):
    for (x1, y1) in st:
        if (2 ** s) * y1 > LIM:
            break
        for (x2, y2) in st:
            if (2 ** s) * y1 * y2 > LIM:
                break
            N = x1 * x2 + (2 ** s) * y1 * y2 - (1 if s == 0 else 0)
            if N <= LIM:
                r[N] = r.get(N, 0) + 1
cand = [N for N in range(3, LIM, 4) if is_prime(N) and is_prime((N + 1) // 4)]
missed = [N for N in cand if r.get(N, 0) == 0]
check("30 primes 4q - 1 < 3000 with q prime; F1 misses 7, 43, 67, 163, 211, 547, 1171, 1723",
      len(cand) == 30 and missed == [7, 43, 67, 163, 211, 547, 1171, 1723], "%d candidates, missed %s" % (len(cand), missed))
check("r_F1(11) = 2, r_F1(19) = 4, r_F1(7) = r_F1(21) = 0", (r.get(11, 0), r.get(19, 0), r.get(7, 0), r.get(21, 0)) == (2, 4, 0, 0),
      str((r.get(11, 0), r.get(19, 0), r.get(7, 0), r.get(21, 0))))
odd_missed_small = [N for N in range(3, 200, 2) if r.get(N, 0) == 0]
print("     odd N < 200 missed by F1:", odd_missed_small)

# ---------------------------------------------------------------- (4) sheets as head and tails
ok = True
for h1 in range(0, 40):
    for h2 in range(0, 40):
        x1, y1 = xy(h1)
        x2, y2 = xy(h2)
        D, U, E = x1 * x2, y1 * y2, 1
        ok &= (x1 * x2 + y1 * y2 - 1 == D - E + U) and all(x1 * x2 + 2 ** s * y1 * y2 == D + 2 ** s * U for s in (1, 2, 3))
        ok &= U >= D
check("F1: sheet s = 0 is D - E + U with E = 1, sheets s >= 1 are D + 2^s U, and U >= D", ok)

print()
if FAILS:
    print("FAILED:", FAILS)
    sys.exit(1)
print("ALL VERIFIER CHECKS PASSED")
