#!/usr/bin/env python3
"""H(P_31) exactly: the compute-node lead "run H(T_31) exactly - falsify/confirm HYP-2371" (kind-pasteur 2026-06-10).

Method (opus-2026-10-14-S34):
  * H(P_q)/|Aff| = N_1 = number of Hamiltonian paths of the Paley tournament P_q starting with the arc 0 -> 1
    (the stabiliser {x -> a x : a a nonzero square} of 0 in Aff acts freely on paths from 0 and transitively on
    their second vertices; S33 note section 5c).
  * N_1 by a LAYERED pull DP over subsets of W = F_q \\ {0, 1} in colex order, two layers in memory (the layout of the
    S33 audit referee's mod-256 DP), here in uint32 residues modulo three primes p_j close to 2^31; the exact N_1 is
    recovered by the Chinese remainder theorem. A successor-permutation injection and Bregman's bound
    certify N_1 <= (d!)^((q-1)/d), d=(q+1)/2, and the exact integer comparison with the CRT modulus.
  * Validation chain (the design doc's acceptance list, paley_H31_compute_design_kpc1.md):
      - the same code must reproduce the exact orbit counts at q = 7, 11, 19, 23 (H(19) = 1172695746915,
        H(23) = 15760206976379349);
      - at q = 27 (GF(27)) and q = 31 the residues mod 256 must match the S33 / audit values (1 and 41);
      - H odd (Redei), 465 | H, H/465 odd (HYP-1713), and the verdict on HYP-2371:
        R(31) = H(T_31) 2^30 / 31! in [2.58949, 2.60249].
Reproduce: python3 04-computation/experiments/paley_H31_exact_20261014.py [qlist=7,11,19,23,27,31]  (q = 31: ~9 GB)
"""
import sys, time, math
from fractions import Fraction
import numpy as np
import numba as nb

QLIST = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [7, 11, 19, 23, 27, 31]
PRIMES = [2147483647, 2147483629, 2147483587]          # primes just below 2^31 (checked below)
FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


@nb.njit(cache=True)
def unrank_colex(r, size, binom):
    mask = np.int64(0)
    for j in range(size, 0, -1):
        s = j - 1
        while binom[s + 1, j] <= r:
            s += 1
        r -= binom[s, j]
        mask |= np.int64(1) << s
    return mask


@nb.njit(parallel=True, cache=True)
def layer_step(f_old, i, m, inmask, binom, nchunks, mod):
    """from layer i (subsets of size i of W, rows in colex order, column = position of the end vertex) to layer i+1"""
    size = i + 1
    total = binom[m, size]
    f_new = np.zeros((total, size), dtype=np.uint32)
    chunk = (total + nchunks - 1) // nchunks
    for c in nb.prange(nchunks):
        r0 = c * chunk
        r1 = min(total, r0 + chunk)
        if r0 >= r1:
            continue
        mask = unrank_colex(r0, size, binom)
        el = np.zeros(size, dtype=np.int64)
        pre = np.zeros(size + 1, dtype=np.int64)
        suf = np.zeros(size + 1, dtype=np.int64)
        for r in range(r0, r1):
            t = 0
            x = mask
            b = 0
            while x:
                if x & 1:
                    el[t] = b
                    t += 1
                x >>= 1
                b += 1
            pre[0] = 0
            for k in range(size):
                pre[k + 1] = pre[k] + binom[el[k], k + 1]
            suf[size - 1] = 0
            for k in range(size - 2, -1, -1):
                suf[k] = suf[k + 1] + binom[el[k + 1], k + 1]
            for k in range(size):
                v = el[k]
                rk = pre[k] + suf[k]                     # colex rank of S minus its k-th element
                acc = np.int64(0)
                im = inmask[v]
                for j in range(size):
                    if j == k:
                        continue
                    u = el[j]
                    if (im >> u) & 1:
                        pos = j if j < k else j - 1
                        acc += np.int64(f_old[rk, pos])
                f_new[r, k] = np.uint32(acc % mod)
            if r + 1 < r1:                                # Gosper: next subset of the same size in colex order
                cbit = mask & (-mask)
                rr = mask + cbit
                mask = (((rr ^ mask) >> 2) // cbit) | rr
    return f_new


def n1_mod(elements, arc, mod):
    W = elements[2:]
    one = elements[1]
    m = len(W)
    inmask = np.zeros(m, dtype=np.int64)
    for a, v in enumerate(W):
        for b, u in enumerate(W):
            if a != b and arc(u, v):
                inmask[a] |= 1 << b
    binom = np.zeros((m + 2, m + 2), dtype=np.int64)
    for n in range(m + 2):
        binom[n, 0] = 1
        for k in range(1, n + 1):
            binom[n, k] = binom[n - 1, k - 1] + (binom[n - 1, k] if k <= n - 1 else 0)
    f = np.zeros((m, 1), dtype=np.uint32)
    for a, v in enumerate(W):
        f[a, 0] = 1 if arc(one, v) else 0
    for i in range(1, m):
        f = layer_step(f, i, m, inmask, binom, 384, np.int64(mod))
    return int(f[0].astype(np.int64).sum() % mod)


def prime_paley(q):
    sq = {(x * x) % q for x in range(1, q)}
    return list(range(q)), (lambda x, y: (y - x) % q in sq)


def gf27_paley():
    """GF(27) = F_3[y]/(y^3 + 2y + 1); element a = a0 + 3 a1 + 9 a2 stands for a0 + a1 y + a2 y^2; 0 = zero, 1 = one"""
    def mul(a, b):
        A = [(a // 3 ** i) % 3 for i in range(3)]
        B = [(b // 3 ** i) % 3 for i in range(3)]
        pr = [0] * 5
        for i in range(3):
            for j in range(3):
                pr[i + j] += A[i] * B[j]
        for k in (4, 3):                               # y^3 = y + 2
            c = pr[k]
            pr[k] = 0
            pr[k - 3] += 2 * c
            pr[k - 2] += c
        return sum((pr[i] % 3) * 3 ** i for i in range(3))

    def sub(a, b):
        return sum(((((a // 3 ** i) % 3) - ((b // 3 ** i) % 3)) % 3) * 3 ** i for i in range(3))
    sq = {mul(x, x) for x in range(1, 27)}
    assert len(sq) == 13
    return list(range(27)), (lambda x, y: sub(y, x) in sq)


def crt(residues, moduli):
    x, M = 0, 1
    for r, m in zip(residues, moduli):
        t = ((r - x) * pow(M, -1, m)) % m
        x += M * t
        M *= m
    return x, M


EXACT_N1 = {7: 9, 11: 1729, 19: 6857869865, 23: 62293308207033}
MOD256 = {3: 1, 7: 9, 11: 193, 19: 41, 23: 185, 27: 1, 31: 41}       # S33 note section 5c (and the S33 audit)


def main():
    for p in PRIMES:
        check(f'prime {p}', is_prime(p))
    print('moduli', PRIMES, 'all prime:', all(is_prime(p) for p in PRIMES))
    rows = {}
    for q in QLIST:
        els, arc = gf27_paley() if q == 27 else prime_paley(q)
        aff = q * (q - 1) // 2
        res = []
        t0 = time.time()
        for p in PRIMES:
            res.append(n1_mod(els, arc, p))
        n1, M = crt(res, PRIMES)
        H = aff * n1
        dt = time.time() - t0
        R = Fraction(H * 2 ** (q - 1), math.factorial(q))
        line = f'q = {q:2d}: N_1 = H/|Aff| = {n1}  (residues {res}; {dt:.1f} s)\n         H(P_{q}) = {H}\n         R(q) = H 2^(q-1)/q! = {float(R):.9f}'
        if q in EXACT_N1:
            line += f'\n         reproduces the known exact value: {check(f"exact q{q}", n1 == EXACT_N1[q])}'
        if q in MOD256:
            line += f'\n         N_1 mod 256 = {n1 % 256} (S33 value {MOD256[q]}): {check(f"mod256 q{q}", n1 % 256 == MOD256[q])}'
        # Paths with fixed 0->1 inject into successor permutations on q-1 rows,
        # each with at most d allowed columns (out-neighbours plus a terminal).
        # Bregman: actual N_1 <= (d!)^((q-1)/d). A small reconstructed residue
        # alone is not a bound for the unreduced count. See wheel_run_interface_20261014h.md.
        d = (q + 1) // 2
        if not check(f'crt uniqueness q{q}', math.factorial(d) ** (q - 1) < M ** d):
            raise ArithmeticError('CRT modulus does not exceed the certified count bound')
        print(line, flush=True)
        rows[q] = (n1, H, R)
    if 31 in rows:
        n1, H, R = rows[31]
        lo, hi = Fraction(258949, 100000), Fraction(260249, 100000)
        Hlo, Hhi = 19830629617139608462365775, 19930130881568868002912737
        print('\n== Acceptance tests for H(T_31) (paley_H31_compute_design_kpc1.md section 0)')
        print('   1. H odd:', check('H odd', H % 2 == 1))
        print('   2. 465 | H:', check('465|H', H % 465 == 0))
        print('   3. H/465 odd (HYP-1713):', check('H/465 odd', (H // 465) % 2 == 1))
        print('   H = 465 mod 930:', H % 930 == 465)
        inside = lo <= R <= hi
        print(f'   4. HYP-2371: R(31) = {float(R):.7f} against the predicted interval [2.58949, 2.60249] (central 2.59599): '
              f'{"INSIDE - CONFIRMED" if inside else "OUTSIDE - REFUTED"};  H in [{Hlo}, {Hhi}]: {Hlo <= H <= Hhi}')
        print(f'   distance from the central value 2.59599: {float(R) - 2.59599:+.6f} (half-width 0.00650)')
        print('   5. the same code reproduced H(19), H(23) exactly:', all(f'exact q{q}' not in FAIL for q in (19, 23) if q in rows))
    print('\nR(q) = H(P_q) 2^(q-1)/q! for all computed q (e = 2.718281828):')
    for q, (n1, H, R) in sorted(rows.items()):
        print(f'   q = {q:2d}: R = {float(R):.9f},  e - R = {math.e - float(R):.6f},  q (1 - R/e) = {q * (1 - float(R) / math.e):.6f}')
    print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
    print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')


if __name__ == '__main__':
    main()
