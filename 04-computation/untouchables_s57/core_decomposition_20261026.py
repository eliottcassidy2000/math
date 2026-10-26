#!/usr/bin/env python3
"""Exact core decomposition of the preimage counts of even targets, by residue class (opus-2026-10-26-S57).

Every preimage n of an even k >= 8 is (i) an even n <= 2k - 2, or (ii) an odd square.  Write an even n as n = a*P with
P = the largest prime factor of n.  If P appears to the first power, the core is a = n/P and k = s(a)P + sigma(a)
(prime-extension identity); otherwise n is "squareful at the top".  For k in the late window [N/10, N] we tally, for
each residue r of k mod l (l = 3, 5, 7, 13), the number of preimages per core a (a <= 64 individually, larger cores
lumped), the odd prime squares (k - 1 prime, "core 1"), the odd composite squares and the squareful tops.

Kill/serve rule (PROVED, elementary): if l does not divide s(a), core a can serve k = sigma(a) (mod l) only through
the prime P = l, i.e. never once k > l s(a) + sigma(a); if l divides s(a), core a serves only k = sigma(a) (mod l).
The table shows which cores carry the nonzero-residue spread of the untouchable share.
Usage: python core_decomposition_20261026.py N
"""
import sys, os, json, math, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from untouchable_sieve_20261026 import primes_upto

LS = np.array([3, 5, 7, 13], np.int64)
AMAX = 64
# class index: 0..AMAX (core a, even a only used; index a), AMAX+1 = big core, AMAX+2 = squareful top,
# AMAX+3 = odd prime square (core 1), AMAX+4 = odd composite square
NC = AMAX + 5


@nb.njit(cache=False)
def seg_cores(lo, cnt, primes, N, kmin, LS, tab):
    rem = np.empty(cnt, np.int64)
    sig = np.ones(cnt, np.int64)
    pmax = np.ones(cnt, np.int64)
    emax = np.zeros(cnt, np.int64)
    for i in range(cnt):
        rem[i] = lo + 2 * i
    for i in range(cnt):
        pk = 1
        s = 1
        e = 0
        while rem[i] % 2 == 0:
            rem[i] //= 2
            pk *= 2
            s += pk
            e += 1
        sig[i] = s
        pmax[i] = 2
        emax[i] = e
    hi = lo + 2 * cnt
    for t in range(1, primes.shape[0]):
        p = primes[t]
        if p * p >= hi:
            break
        step = 2 * p
        first = ((lo + step - 1) // step) * step
        for m in range(first, hi, step):
            i = (m - lo) // 2
            pk = 1
            s = 1
            e = 0
            while rem[i] % p == 0:
                rem[i] //= p
                pk *= p
                s += pk
                e += 1
            sig[i] *= s
            pmax[i] = p
            emax[i] = e
    for i in range(cnt):
        n = lo + 2 * i
        r = rem[i]
        sg = sig[i]
        if r > 1:
            sg *= (r + 1)
        k = sg - n
        if k < kmin or k > N or (k & 1) == 1:
            continue
        if r > 1:
            a = n // r
            c = a if a <= 64 else 65
        else:
            if emax[i] == 1:
                a = n // pmax[i]
                c = a if a <= 64 else 65
            else:
                c = 66
        for t in range(LS.shape[0]):
            tab[t, k % LS[t], c] += 1


@nb.njit(cache=False)
def odd_squares(N, kmin, LS, tab):
    Q = int(N ** 0.5) + 2
    # odd prime q: k = q + 1 (need primes up to N); done by a simple segmented sieve
    L = int(N ** 0.5) + 2
    base = np.ones(L + 1, np.uint8)
    base[:2] = 0
    for i in range(2, int(L ** 0.5) + 1):
        if base[i]:
            base[i * i::i] = 0
    bp = np.nonzero(base)[0]
    seg = 1 << 22
    lo = 3
    while lo < N:
        hi = min(lo + seg, N)
        fl = np.ones(hi - lo, np.uint8)
        for t in range(bp.shape[0]):
            p = bp[t]
            if p * p >= hi:
                break
            st = max(p * p, ((lo + p - 1) // p) * p)
            for m in range(st, hi, p):
                fl[m - lo] = 0
        for i in range(hi - lo):
            q = lo + i
            if fl[i] and (q & 1) == 1:
                k = q + 1
                if kmin <= k <= N:
                    for t in range(LS.shape[0]):
                        tab[t, k % LS[t], 68] += 1
        lo = hi
    Qc = int(N ** (2.0 / 3.0)) + 2
    for q in range(9, Qc + 1, 2):
        x = q
        s = 1
        p = 3
        comp = False
        while p * p <= x:
            if x % p == 0:
                comp = True
                e = 0
                while x % p == 0:
                    x //= p
                    e += 1
                tt = 1
                pk = 1
                for _ in range(2 * e):
                    pk *= p
                    tt += pk
                s *= tt
            p += 2
        if not comp:
            continue
        if x > 1:
            s *= (1 + x + x * x)
        k = s - q * q
        if kmin <= k <= N and (k & 1) == 0:
            for t in range(LS.shape[0]):
                tab[t, k % LS[t], 69] += 1


def sigma(n):
    s, d = 0, 1
    while d * d <= n:
        if n % d == 0:
            s += d
            if d * d != n:
                s += n // d
        d += 1
    return s


if __name__ == '__main__':
    N = int(float(sys.argv[1]))
    kmin = N // 10
    t0 = time.time()
    primes = primes_upto(int(math.isqrt(2 * N)) + 2)
    tab = np.zeros((LS.shape[0], 13, 70), np.int64)
    M = 2 * N - 2
    lo = 2
    while lo <= M:
        cnt = min(1 << 22, (M - lo) // 2 + 1)
        seg_cores(lo, cnt, primes, N, kmin, LS, tab)
        lo += 2 * cnt
    odd_squares(N, kmin, LS, tab)
    print(f'tabulated ({time.time() - t0:.0f}s)', flush=True)
    # number of even k per residue in the window, and untouchable shares from the exact array if available
    out = {'N': N, 'window': [kmin, N]}
    names = {68: 'odd prime square (core 1)', 69: 'odd composite square', 65: 'core > 64', 66: 'squareful top'}
    for t, l in enumerate(LS):
        nk = np.zeros(l, np.int64)
        for r in range(l):
            # even k in [kmin, N] with k = r mod l
            cntr = 0
            for k0 in range(kmin + (kmin & 1), kmin + (kmin & 1) + 2 * l, 2):
                if k0 % l == r:
                    cntr = (N - k0) // (2 * l) + 1 if k0 <= N else 0
            nk[r] = cntr
        rows = {}
        for c in list(range(2, 65, 2)) + [65, 66, 68, 69]:
            per = [tab[t, r, c] / nk[r] for r in range(l)]
            if sum(tab[t, :, c]) == 0:
                continue
            nm = names.get(c, f'core {c}')
            if c <= 64:
                sa = sigma(c) - c
                rule = (f'serves only r = {sigma(c) % l}' if sa % l == 0 else f'kills r = {sigma(c) % l}')
                nm += f' (s={sa}, sigma={sigma(c)}; {rule})'
            rows[nm] = [round(x, 6) for x in per]
        tot = [tab[t, r, :].sum() / nk[r] for r in range(l)]
        out[int(l)] = {'mean_multiplicity_by_residue': [round(x, 5) for x in tot], 'per_core': rows}
        print(f'\nmod {l}: mean multiplicity by residue r = 0..{l - 1}: ' + ' '.join(f'{x:.4f}' for x in tot))
        if l in (3, 5, 7):
            for nm, per in rows.items():
                if max(per) >= 0.004:
                    print(f'   {nm:<60} ' + ' '.join(f'{x:.4f}' for x in per))
    json.dump(out, open(os.path.join(HERE, f'core_decomposition_{N}.json'), 'w'), indent=1)
    print(f'done ({time.time() - t0:.0f}s)')
