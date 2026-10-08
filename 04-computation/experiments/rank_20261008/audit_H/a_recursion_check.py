#!/usr/bin/env python3
"""audit H, task A: THM-4611 statement 1 by brute force over all digit words with exact Fractions.

For random exact states (M in Gamma, e an integer or an element of Z_(d) with denominator prime to d) and k = 2, 3, 4:
  S1 = sum over all d^k words of eta eta^T               (eta = debt change over the k steps, from exact M's)
  S2 = sum_{s<k} d^(k-1-s) sum_{words of length s} D(pi_s)  (definition of d^k * sum_s E[C_pi(t+s)])
  S3 = A_k(h), h = (M mod d^k, e mod d^k), from the residue recursion
must agree; sum of eta over words must be 0.  Also: the coupling sequence along every word of length k, and
(M', e') mod d^(s-1), depend only on (M, e) mod d^s (checked by perturbing M by g in Gamma with g = 1 mod d^k and
e by multiples of d^k / random Z_(d) multiples of d^k).
Usage: python3 a_recursion_check.py [nstates]"""
import itertools, random, sys, time
from fractions import Fraction as Fr
from hcore import Map, expvec, residue

def run(mp, M0, e0, k):
    d, n = mp.d, mp.n
    x0 = expvec(M0, mp.primes)
    S1 = [[0] * n for _ in range(n)]
    S2 = [[0] * n for _ in range(n)]
    mean = [0] * n
    coupl = {}
    # walk the word tree level by level, keeping exact states
    level = [((), M0, e0)]
    for s in range(k):
        nxt = []
        for w, M, e in level:
            a, b = residue(M, d), residue(e, d)
            coupl[w] = (a, b)
            D = mp.Dmat(a, b)
            for p in range(n):
                for q in range(n):
                    S2[p][q] += d ** (k - 1 - s) * D[p][q]
            for j in range(d):
                M2, e2, i = mp.step(M, e, j)
                # debt change equals the root v_i - v_j (independent check via factorisation)
                dx = tuple(a1 - b1 for a1, b1 in zip(expvec(M2, mp.primes), expvec(M, mp.primes)))
                assert dx == tuple(p1 - q1 for p1, q1 in zip(mp.v[i], mp.v[j])), (dx, i, j)
                nxt.append((w + (j,), M2, e2))
        level = nxt
    for w, M, e in level:
        eta = [a1 - b1 for a1, b1 in zip(expvec(M, mp.primes), x0)]
        for p in range(n):
            mean[p] += eta[p]
            for q in range(n):
                S1[p][q] += eta[p] * eta[q]
    S3 = mp.A_rec(k, residue(M0, d ** k), residue(e0, d ** k), {})
    return S1, S2, S3, mean, coupl

def rand_state(mp, rnd, kind):
    # M: random element of Gamma (product of random powers of the ratios m_i/m_0)
    M = Fr(1)
    for i in range(1, mp.d):
        M *= Fr(mp.m[i], mp.m[0]) ** rnd.randint(-3, 3)
    if kind == 0:
        e = Fr(rnd.randint(-10 ** 6, 10 ** 6))
    else:
        q = rnd.choice([3, 7, 9, 11, 13, 21, 77, 2, 4, 6]) if mp.d % 2 else rnd.choice([3, 5, 7, 9, 11])
        while q % mp.d == 0 or any(q % p == 0 for p in range(2, mp.d + 1) if mp.d % p == 0):
            q += 1
        e = Fr(rnd.randint(-10 ** 6, 10 ** 6), q)
    return M, e

if __name__ == '__main__':
    nst = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    rnd = random.Random(4611)
    maps = [(5, [1, 2, 3, 7, 1], [0, 3, 4, 4, 1]), (5, [1, 1, 2, 3, 7], None), (7, [1, 1, 1, 2, 3, 5, 11], None),
            (4, [1, 3, 5, 7], None)]
    t0 = time.time()
    for d, m, r in maps:
        mp = Map(d, m, r)
        print(f"Z_{d} m = {m} r = {mp.r} primes {mp.primes} rank {mp.rank}", flush=True)
        ks = (2, 3, 4) if d ** 4 <= 2401 and d <= 5 else (2, 3)
        nbad = 0; ntest = 0
        for k in ks:
            for t in range(nst):
                M0, e0 = rand_state(mp, rnd, t % 2)
                S1, S2, S3, mean, coupl = run(mp, M0, e0, k)
                ok = (S1 == S2 == S3) and not any(mean)
                # perturbation: same residues mod d^k, different exact state
                g = Fr(1)
                # find g in Gamma with g = 1 mod d^k: a power of the first nontrivial ratio
                base = next(Fr(mp.m[i], mp.m[0]) for i in range(1, d) if mp.m[i] != mp.m[0])
                pw = 1
                while residue(base ** pw, d ** k) != 1:
                    pw += 1
                g = base ** pw
                M1 = M0 * g
                e1 = e0 + Fr(d ** k * rnd.randint(-50, 50), rnd.choice([1, 3]) if d % 3 else 1)
                S1b, S2b, S3b, meanb, couplb = run(mp, M1, e1, k)
                ok2 = (S1b == S1) and (couplb == coupl)
                ntest += 1
                if not (ok and ok2):
                    nbad += 1
                    print("   MISMATCH", k, M0, e0, S1, S2, S3, mean, ok2, flush=True)
            print(f"   k = {k}: {nst} random exact states (+{nst} perturbed twins): S1 = S2 = A_k(h) and mean 0: "
                  f"{'all OK' if nbad == 0 else str(nbad) + ' failures'}  [{time.time() - t0:.1f}s]", flush=True)

# ---- special states: identity runs, M = 1 mod d^k with M != 1, e = 0; and k = 5 on Z_5 ----
def special():
    rnd = random.Random(7)
    out = []
    for d, m, r in [(5, [1, 2, 3, 7, 1], [0, 3, 4, 4, 1]), (5, [1, 2, 3, 7, 11], None), (4, [1, 3, 5, 7], None)]:
        mp = Map(d, m, r)
        for k in ((2, 3, 4, 5) if d == 5 else (2, 3, 4, 5, 6)):
            base = next(Fr(mp.m[i], mp.m[0]) for i in range(1, d) if mp.m[i] != mp.m[0])
            pw = 1
            while residue(base ** pw, d ** k) != 1:
                pw += 1
            g = base ** pw          # g in Gamma, g = 1 mod d^k, g != 1
            states = [(Fr(1), Fr(b * d ** (k - 1))) for b in range(1, d)] + [(g, Fr(0)), (g, Fr(d ** k)), (g, Fr(d ** (k - 1))),
                      (Fr(1), Fr(d ** k)), (1 / g, Fr(3 * d ** (k - 1), 1 if d % 3 == 0 else 1))]
            bad = 0
            for M0, e0 in states:
                S1, S2, S3, mean, coupl = run(mp, M0, e0, k)
                if not (S1 == S2 == S3 and not any(mean)):
                    bad += 1
                    print("   MISMATCH special", d, m, k, M0, e0, flush=True)
            print(f"special states Z_{d} m = {m} k = {k}: {len(states)} states, mismatches {bad}", flush=True)

special()
