#!/usr/bin/env python3
"""The derivative rules D_c(S) = {x : [x in S] != [cx in S]} on circulant tournaments of prime order.
K_max(p) = odd part of Z_p^* (the 2^e-th powers, 2^e || p-1). Facts checked: (i) if c = -1 mod K_max then
D_c(S) = Z_p^* (the odd complete graph) for every K_max-invariant tournament S; if c in K_max then D_c(S) = empty.
(ii) For every c != 1 the primes with c = -1 mod K_max(p) seem infinite (non-squares: p = 3 mod 4 with (c/p) = -1;
squares d^2: p = 5 mod 8 with (d/p) = -1; etc.). Table for 2 <= c <= 12 and primes p < 400 with K_max != 1."""
from sympy import isprime, primitive_root
def kmax(p):
    e = 0; m = p - 1
    while m % 2 == 0: m //= 2; e += 1
    return {pow(x, 2 ** e, p) for x in range(1, p)}, e
rows = {}
for p in [q for q in range(5, 400) if isprime(q)]:
    K, e = kmax(p)
    if len(K) == 1: continue
    minusK = {(p - k) % p for k in K}
    for c in range(2, 13):
        cc = c % p
        tag = 'K' if cc in K else ('-K' if cc in minusK else '.')
        rows.setdefault(c, []).append((p, e, tag))
for c in range(2, 13):
    bad = [p for p, e, t in rows[c] if t == '-K']
    good = [p for p, e, t in rows[c] if t == 'K']
    print(f"c = {c:2d}: primes p < 400 (K_max != 1) with c = -1 mod K_max (D_c sends the K_max-symmetric tournaments to K_p): "
          f"{len(bad)} of {len(rows[c])}, e.g. {bad[:12]}; with c in K_max (sent to empty): {len(good)}, e.g. {good[:8]}")
# direct check of (i) for small p by brute force on all K_max-invariant tournament sets
import itertools
for p in [7, 11, 13, 19, 23, 29, 31, 37, 41, 43]:
    K, e = kmax(p); units = list(range(1, p))
    cosets = []; seen = set()
    for x in units:
        if x in seen: continue
        C = frozenset(x * k % p for k in K); seen |= C; cosets.append(C)
    # pair cosets with their negatives
    pairs = []; used = set()
    for C in cosets:
        if C in used: continue
        negC = frozenset((p - x) % p for x in C); used |= {C, negC}; pairs.append((C, negC))
    ok = True
    for choice in itertools.product([0, 1], repeat=len(pairs)):
        S = frozenset().union(*[pr[ch] for pr, ch in zip(pairs, choice)])
        for c in range(2, 13):
            D = frozenset(x for x in units if (x in S) != ((c * x) % p in S))
            cc = c % p
            if cc in K: ok &= (D == frozenset())
            elif cc in {(p - k) % p for k in K}: ok &= (D == frozenset(units))
    print(f"p = {p} (e = {e}, |K_max| = {len(K)}): fact (i) verified on all {2 ** len(pairs)} K_max-invariant tournament sets: {ok}")
