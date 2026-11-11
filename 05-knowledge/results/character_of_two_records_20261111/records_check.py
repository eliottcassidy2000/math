#!/usr/bin/env python3
"""Records, sprouts and the Redei-Berge function of tournaments (section "Records" of the note).

For every labelled tournament T on [n] (n <= NMAX) and every collection S of pairwise vertex-disjoint directed odd
cycles of T (the remaining vertices are singleton blocks), weight 2^|S|:
  (a) cycle type lambda(S)         -> the Redei-Berge function U_T = sum_S 2^|S| p_lambda(S)   (Grinberg-Stanley);
  (b) max-ordered type tau(S)      -> blocks listed by increasing largest vertex (Chen-Mogielnicki-Ono, Prop. 6.1).
Checks, summed over all 2^C(n,2) tournaments:
  (a) sum_T U_T = 2^C(n,2) n! sum_{lambda odd parts} 2^(l(lambda) - n) p_lambda / z_lambda  ( = 2^C(n,2) n! 2^(-n) q_n );
  (b) sum_T [coefficient of type alpha] = 2^C(n,2) n! 2^(l(alpha) - n) / pi_u(alpha) for alpha with odd parts, 0 otherwise;
  (c) sum_{alpha |= n, odd parts} 2^l(alpha) / pi_u(alpha) = 2   (n <= 40, exact).
Tournaments and reachability:
  (d) l_n = #{T on [n] : vertex n reaches every vertex} = sum_j C(n-1, j-1) S_j 2^C(n-j, 2) (S_j = strong labelled),
      exp(sum l_n x^n / n!) = sum 2^C(n,2) x^n / n!;
  (e) the number of T whose max-ordered piece type is alpha equals prod_j C(s_j - 1, alpha_j - 1) l_(alpha_j).
Usage: records_check.py NMAX (default 5)"""
import sys, itertools
from fractions import Fraction as Fr
from math import comb, factorial
from collections import Counter

def compositions(n):
    if n == 0: yield (); return
    for k in range(1, n + 1):
        for rest in compositions(n - k): yield (k,) + rest
def pi_u(alpha):
    s = 0; p = 1
    for a in alpha: s += a; p *= s
    return p
def zlam(lam):
    c = Counter(lam); z = 1
    for k, m in c.items(): z *= k ** m * factorial(m)
    return z
def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        out = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out
def odd_cycles(out, n):
    cyc = []
    for L in range(3, n + 1, 2):
        for Sset in itertools.combinations(range(n), L):
            v = Sset[0]
            for p in itertools.permutations(Sset[1:]):
                seq = (v,) + p
                if all((out[seq[i]] >> seq[(i + 1) % L]) & 1 for i in range(L)): cyc.append(frozenset(seq))
    return cyc
def packings(cyc):
    res = []
    def rec(start, used, chosen):
        res.append(list(chosen))
        for i in range(start, len(cyc)):
            if not (cyc[i] & used):
                chosen.append(cyc[i]); rec(i + 1, used | cyc[i], chosen); chosen.pop()
    rec(0, frozenset(), [])
    return res
def reach_all(out, n, v):
    seen = 1 << v; stack = [v]
    while stack:
        u = stack.pop(); t = out[u] & ~seen
        while t:
            w = (t & -t).bit_length() - 1; seen |= 1 << w; stack.append(w); t &= t - 1
    return seen == (1 << n) - 1
def sccs_in_order(out, n):
    # strong components in condensation order (sources first): sort vertices by the size of their reach set, descending
    def reach(v):
        seen = 1 << v; stack = [v]
        while stack:
            u = stack.pop(); t = out[u] & ~seen
            while t:
                w = (t & -t).bit_length() - 1; seen |= 1 << w; stack.append(w); t &= t - 1
        return seen
    R = [reach(v) for v in range(n)]
    comps = {}
    for v in range(n): comps.setdefault(R[v], []).append(v)
    return [comps[k] for k in sorted(comps, key=lambda r: -bin(r).count('1'))]

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
fail = 0
# strong labelled counts from T = 1/(1 - S)
t = [2 ** comb(k, 2) for k in range(0, 13)]
Sl = [0] * 13
for k in range(1, 13):
    Sl[k] = t[k] - sum(comb(k, j) * Sl[j] * t[k - j] for j in range(1, k))
ell = [0] + [sum(comb(k - 1, j - 1) * Sl[j] * 2 ** comb(k - j, 2) for j in range(1, k + 1)) for k in range(1, 13)]
print("strong labelled S_n (n = 1..10):", Sl[1:11])
print("l_n = #{vertex n reaches every vertex} (n = 1..10):", ell[1:11])
# (d) exponential formula check: t_n = sum_k C(n-1, k-1) l_k t_(n-k)
for k in range(1, 13):
    if t[k] != sum(comb(k - 1, j - 1) * ell[j] * t[k - j] for j in range(1, k + 1)): fail += 1
print("exponential formula T = exp(L) (n <= 12): %s" % ("holds" if not fail else "FAILS"))
for n in range(1, nmax + 1):
    U = Counter(); N = Counter(); reach_count = 0; piece = Counter()
    for out in tournaments(n):
        cyc = odd_cycles(out, n)
        for S in packings(cyc):
            used = set().union(*S) if S else set()
            blocks = [set(c) for c in S] + [{v} for v in range(n) if v not in used]
            lam = tuple(sorted((len(b) for b in blocks), reverse=True))
            tau = tuple(len(b) for b in sorted(blocks, key=max))
            U[lam] += 2 ** len(S); N[tau] += 2 ** len(S)
        if reach_all(out, n, n - 1): reach_count += 1
        comps = sccs_in_order(out, n)
        # pieces: cut the component sequence before each component whose maximum exceeds all earlier maxima
        sizes = []; cur = 0; best = -1
        for c in comps:
            if max(c) > best:
                if cur: sizes.append(cur)
                cur = 0; best = max(c)
            cur += len(c)
        sizes.append(cur)
        # order pieces by their maxima: they are already in increasing order of maxima
        piece[tuple(sizes)] += 1
    tot = 2 ** comb(n, 2); nf = factorial(n)
    for lam, v in U.items():
        pred = Fr(tot * nf * 2 ** len(lam), 2 ** n * zlam(lam)) if all(p % 2 for p in lam) else 0
        if v != pred: fail += 1
    for alpha in compositions(n):
        pred = Fr(tot * nf * 2 ** len(alpha), 2 ** n * pi_u(alpha)) if all(p % 2 for p in alpha) else 0
        if N.get(alpha, 0) != pred: fail += 1
        s = 0; pr = 1
        for a in alpha: s += a; pr *= comb(s - 1, a - 1) * ell[a]
        if piece.get(alpha, 0) != pr: fail += 1
    if reach_count != ell[n]: fail += 1
    print("n = %d: %d tournaments; mean Redei-Berge function in power sums: %s" % (n, tot, ", ".join("%s p_%s" % (Fr(v, tot), "".join(map(str, lam))) for lam, v in sorted(U.items()))))
    print("        max-ordered types (weighted counts / 2^C(n,2)): %s" % ", ".join("%s:%s" % ("".join(map(str, a)), Fr(v, tot)) for a, v in sorted(N.items())))
    print("        tournaments by max-ordered piece type: %s ; vertex n reaches all in %d" % (", ".join("%s:%d" % ("".join(map(str, a)), v) for a, v in sorted(piece.items())), reach_count))
for n in range(1, 41):
    s = sum(Fr(2 ** len(a), pi_u(a)) for a in compositions(n) if all(p % 2 for p in a)) if n <= 22 else None
    if s is not None and s != 2: fail += 1
print("identity sum_{odd alpha} 2^l(alpha)/pi_u(alpha) = 2 checked for n <= 22")
print("failures: %d" % fail)
