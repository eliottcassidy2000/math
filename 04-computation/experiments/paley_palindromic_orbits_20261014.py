#!/usr/bin/env python3
"""Paley Hamiltonian paths split by the reversal involution (THM-582): H(P_q)/|Aff| = 2F + P.

For a prime q = 3 mod 4 the Paley tournament P_q (i -> j iff j - i is a nonzero square) has the involutory
anti-automorphisms phi_b(x) = -x + b, and rho_b(P) = phi_b(reverse(P)) is an involution on Hamiltonian paths whose fixed
points are the phi_b-palindromic paths (THM-582).  A path fixed by two of them would be fixed by the translation
phi_b phi_b', impossible since automorphisms act freely on Hamiltonian paths; so the extended group G = Aff x| <rho_0>
(order q(q-1)) has orbits of size q(q-1) (F of them) or q(q-1)/2 (the P orbits of palindromic-type paths), and

    H = q(q-1) F + q(q-1)/2 P,   P = 2 pal/(q-1),   pal = #{phi_0-palindromic Hamiltonian paths},

so H/|Aff| = 2F + P and HYP-9273 (B) reads 2F + P = 1 (mod 8).  A phi_0-palindromic path has middle vertex 0 and is
determined by its first half, which takes one element from each pair {x, -x} and ends at a vertex v with v -> 0;
so pal is a Hamiltonian-path count on (q-1)/2 pairs: a cheap dynamic programme (2^((q-1)/2) subsets).

Uses: (1) an independent consistency check of H(P_19), H(P_23) (F must be a nonnegative integer); (2) P mod 8 at
q = 31, 43, 47, where H is out of reach, so HYP-9273 (B) there becomes a statement about F mod 4 alone.

Session opus-2026-10-14-S33.  Note: 05-knowledge/results/number_tournament_bridges_20261014.md
Reproduce: python3 04-computation/experiments/paley_palindromic_orbits_20261014.py   (about 1 minute, 3.2 GB at q = 47)
"""
import time
import numpy as np
import numba as nb
import sympy

FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


@nb.njit(cache=True)
def pal_count(q, isqr, M):
    """number of phi_0-palindromic Hamiltonian paths of P_q, modulo M (M = 0: exact int64)"""
    m = (q - 1) // 2
    # vertices 1..q-1; pair index of x: index of min(x, q-x) among 1..m
    pair = np.zeros(q, dtype=np.int64)
    for x in range(1, q):
        y = x if x <= m else q - x
        pair[x] = y - 1
    full = (1 << m) - 1
    dp = np.zeros(((1 << m), q), dtype=np.int64)
    for v in range(1, q):
        dp[1 << pair[v], v] = 1
    for mask in range(1, full + 1):
        for v in range(1, q):
            c = dp[mask, v]
            if c == 0:
                continue
            for w in range(1, q):
                if (mask >> pair[w]) & 1:
                    continue
                if isqr[(w - v) % q]:
                    nm = mask | (1 << pair[w])
                    if M > 0:
                        dp[nm, w] = (dp[nm, w] + c) % M
                    else:
                        dp[nm, w] += c
    tot = 0
    for v in range(1, q):
        if isqr[(q - v) % q]:            # v -> 0
            tot += dp[full, v]
            if M > 0:
                tot %= M
    return tot


H_KNOWN = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349}


def main():
    print('Paley tournaments P_q, q prime = 3 mod 4:  H/|Aff| = 2F + P,  P = 2 pal/(q-1)')
    rows = []
    for q in (3, 7, 11, 19, 23, 31, 43, 47):
        t0 = time.time()
        isqr = np.zeros(q, dtype=np.bool_)
        for x in range(1, q):
            isqr[x * x % q] = True
        exact = q <= 31
        M = 0 if exact else 8 * (q - 1) // 2 * 1024       # pal mod 4(q-1) suffices for P mod 8; keep extra bits
        pal = int(pal_count(q, isqr, M))
        aff = q * (q - 1) // 2
        if exact:
            assert (2 * pal) % (q - 1) == 0
            P = 2 * pal // (q - 1)
            Pm8 = P % 8
        else:
            # pal = P (q-1)/2 exactly; modulo M = 4(q-1) * 1024 the residue determines P mod 8 * 1024
            assert (2 * pal) % (q - 1) == 0
            P = None
            Pm8 = (2 * pal // (q - 1)) % 8
        line = f'   q = {q:2d}: pal = {pal if exact else str(pal) + " (mod " + str(M) + ")"}, P = {P if exact else "?"} (P mod 8 = {Pm8})'
        if q in H_KNOWN:
            H = H_KNOWN[q]
            num = H - q * pal
            ok = num >= 0 and num % (q * (q - 1)) == 0
            F = num // (q * (q - 1)) if ok else None
            orb = H // aff
            line += f', H/|Aff| = {orb} = 2F + P with F = {F}: {check(f"F int {q}", ok and orb == 2 * F + P)};  F mod 4 = {F % 4}, orbit mod 8 = {orb % 8}'
        else:
            line += f';  HYP-9273 (B) at q = {q} <=> 2F + P = 1 (mod 8) <=> F = {((1 - Pm8) % 8) // 2} (mod 4)'
        print(line + f'   [{time.time() - t0:.1f} s]')
        rows.append((q, pal, P, Pm8))
    check('THM-582 P7', rows[1][1] == 9)
    check('THM-582 P11', rows[2][1] == 185)
    print('   (THM-582 lists 1, 9, 185 palindromic paths for q = 3, 7, 11)')
    print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
    print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')


if __name__ == '__main__':
    main()
