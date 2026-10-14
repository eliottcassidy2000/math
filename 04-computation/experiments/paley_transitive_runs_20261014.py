#!/usr/bin/env python3
"""Paley tournaments, consecutive residues and the least quadratic non-residue (a 'linked list' inside P_p).

For a prime p = 3 mod 4 and the Paley tournament P_p (i -> j iff j - i is a nonzero square), the run {0, 1, ..., n(p)-1}
is a transitive subtournament, where n(p) is the least quadratic non-residue: all differences 1, ..., n(p)-1 are squares.
So the largest transitive subtournament trans(P_p) is at least n(p) (elementary; the Paley-graph analogue
omega(G_p) >= n(p) is classical).  Consequences, all CITED:
  * Erdos-Moser: every tournament on p vertices has trans >= floor(log2 p) + 1; random tournaments have trans about
    2 log2 p.
  * Graham-Ringrose (1990): n(p) >> log p log log log p for infinitely many p, so trans(P_p)/log p is unbounded:
    Paley tournaments are not Erdos-Moser-extremal infinitely often.
  * GRH (Ankeny 1952; Bach 1990): n(p) < 2 log^2 p, so the consecutive run is short.  The discrepancy bound
    |e(A,B) - |A||B|/2| <= (1/2) sqrt(|A||B| p) only gives trans(P_p) <= 2 sqrt(p) + O(1).
This script computes n(p) and trans(P_p) exactly (branch and bound; translation lets the source be 0) for p = 3 mod 4
up to PMAX, and checks trans >= n(p) and trans >= floor(log2 p) + 1.  (Only the INITIAL run of residues 1, ..., n(p)-1
gives a transitive set: transitivity needs all pairwise differences to be residues.)
Session opus-2026-10-14-S34.  Reproduce: python3 04-computation/experiments/paley_transitive_runs_20261014.py [PMAX=400]
"""
import sys, math, time

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 400
FAIL = []


def check(label, ok):
    if not ok:
        FAIL.append(label)
    return ok


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def trans_paley(p):
    """largest transitive subtournament.  A transitive set X has a unique source s(X), and X - s(X) lies in the
    out-neighbourhood of s(X); so every transitive set is reached exactly once by 'choose the source, recurse into its
    out-neighbourhood' (no vertex may be discarded after its branch: unlike cliques, chains are ordered).  By
    translation the global source may be 0.  Bound: trans(C) <= 1 + max_v |C cap out(v)|."""
    qr = {x * x % p for x in range(1, p)}
    out = [sum(1 << ((v + r) % p) for r in qr) for v in range(p)]
    best = [1]

    def extend(size, cand):
        if size > best[0]:
            best[0] = size
        if not cand:
            return
        degs = []
        c = cand
        while c:
            low = c & -c
            v = low.bit_length() - 1
            c ^= low
            degs.append((bin(cand & out[v]).count('1'), v))
        degs.sort(reverse=True)
        if size + 1 + degs[0][0] <= best[0]:
            return
        for d, v in degs:
            if size + 1 + d <= best[0]:
                break
            extend(size + 1, cand & out[v])
    extend(1, out[0])
    return best[0]


def trans_bruteforce(p):
    """exhaustive: largest acyclic vertex set, by a DP over all subsets (p <= 19)"""
    qr = {x * x % p for x in range(1, p)}
    out = [sum(1 << ((v + r) % p) for r in qr) for v in range(p)]
    N = 1 << p
    ok = bytearray(N)
    ok[0] = 1
    best = 0
    for X in range(1, N):
        x = X
        while x:
            low = x & -x
            s = low.bit_length() - 1
            x ^= low
            rest = X ^ low
            if (rest & out[s]) == rest and ok[rest]:
                ok[X] = 1
                break
        if ok[X]:
            best = max(best, bin(X).count('1'))
    return best


def least_nonresidue(p):
    for a in range(2, p):
        if pow(a, (p - 1) // 2, p) == p - 1:
            return a


def main():
    ps = [p for p in primes_upto(PMAX) if p % 4 == 3]
    print(f'Paley tournaments P_p, p = 3 mod 4, p <= {PMAX}')
    print('    p   n(p)  trans(P_p)  floor(log2 p)+1  2 log2 p   2 log^2 p (Bach, GRH)')
    t0 = time.time()
    for p in ps:
        n = least_nonresidue(p)
        tr = trans_paley(p)
        lo = int(math.floor(math.log2(p))) + 1
        check(f'trans>=n {p}', tr >= n)
        check(f'EM {p}', tr >= lo)
        if p <= 19:
            check(f'brute {p}', tr == trans_bruteforce(p))
        print(f'  {p:4d}  {n:4d}  {tr:6d}       {lo:6d}         {2 * math.log2(p):6.2f}   {2 * math.log(p) ** 2:8.1f}')
    print(f'  ({time.time() - t0:.1f} s)')
    print('\nFAILED CHECKS:', FAIL if FAIL else 'none')
    print('ALL CHECKS PASSED' if not FAIL else 'SOME CHECKS FAILED')


if __name__ == '__main__':
    main()
