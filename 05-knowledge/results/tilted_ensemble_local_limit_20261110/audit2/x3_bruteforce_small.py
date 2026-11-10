#!/usr/bin/env python3
"""AUDIT 2 (auditor's own code).  Third, DP-free layer: every Hamiltonian-path count below is obtained by testing all
permutations; cycles and paths are enumerated as vertex sequences.  All labelled tournaments, exact rationals.

For every labelled R on m = 2..5 vertices:
  H(R + v_X) for all X (permutations of m+1 vertices); a_i, mu_ij from first / second differences; (4.2) on every X;
  Lemma 4.5: sum mu = sum a = sum_S H(R[S]) H(R[V-S]) = 2 mod 4 (the split sum from permutation counts of the induced
  subtournaments); the same total = (m-1) H(R) + F_1(R) (THM-4111 (6)); the number of directed paths on an even number of vertices is odd; (4.4) beta = P_ij + P_ji mod 2;
  b_R, Corollary 4.6 (a); E_R[b_R^2], E_R|b_R|, D_(m+1).
For every labelled T on n = 4..6 vertices (v = 0, v' = 1): both right-hand sides of (4.3),
  (-1)^((H(T-v) - H(T-v'))/2)  and  (-1)^(N_1(T)), N_1 by enumeration of the odd cycles through exactly one of v, v'."""
from fractions import Fraction
from itertools import combinations, permutations

def arcs_index(n):
    prs = list(combinations(range(n), 2)); return prs, {p: k for k, p in enumerate(prs)}

def perm_tests(verts, idx):
    """for every ordering of verts: (care mask, wanted bits) so that the ordering is a path iff code & care == want;
    bit of pair (i<j) set means i -> j"""
    out = []
    for p in permutations(verts):
        care = want = 0
        for a, b in zip(p, p[1:]):
            k = idx[(min(a, b), max(a, b))]; care |= 1 << k
            if a < b: want |= 1 << k
        out.append((care, want))
    return out

def beats(code, idx, a, b):
    k = idx[(min(a, b), max(a, b))]; return bool(code >> k & 1) == (a < b)

def part1(m):
    n = m + 1; prs, idx = arcs_index(n); v = m
    tests_full = perm_tests(range(n), idx)
    sub_tests = {S: perm_tests([x for x in range(m) if S >> x & 1], idx) for S in range(1, 1 << m)}
    prsR = [p for p in prs if p[1] < m]
    fails = dict(q42=0, l45=0, thm4111=0, thm4184=0, evenpaths=0, f44=0, cor46a=0, parity=0, redei=0)
    sum_b2 = Fraction(0); sum_abs = Fraction(0); sum_D = Fraction(0); cnt = 0
    for rc in range(1 << len(prsR)):
        base = 0
        for t, p in enumerate(prsR):
            if rc >> t & 1: base |= 1 << idx[p]
        Hx = []
        for X in range(1 << m):
            code = base
            for i in range(m):
                if not X >> i & 1: code |= 1 << idx[(i, v)]          # i -> v iff i notin X ; v -> i iff i in X
            Hx.append(sum(1 for care, want in tests_full if code & care == want))
        hS = {0: 1}
        for S in range(1, 1 << m): hS[S] = sum(1 for care, want in sub_tests[S] if base & care == want)
        full = (1 << m) - 1; h0 = Hx[0]
        if h0 != hS[full] or Hx[full] != h0 or any(h % 2 == 0 for h in Hx): fails['redei'] += 1
        a = [Hx[1 << i] - h0 for i in range(m)]
        mu = {(i, j): Hx[1 << i] + Hx[1 << j] - Hx[(1 << i) | (1 << j)] - h0 for i, j in combinations(range(m), 2)}
        if any(x % 2 for x in a) or any(x % 2 or x < 2 for x in mu.values()): fails['parity'] += 1
        for X in range(1 << m):
            pred = h0 + sum(a[i] for i in range(m) if X >> i & 1) - sum(u for (i, j), u in mu.items() if X >> i & 1 and X >> j & 1)
            if pred != Hx[X]: fails['q42'] += 1; break
        split = sum(hS[S] * hS[full ^ S] for S in range(1, full))
        if not (sum(mu.values()) == sum(a) == split and split % 4 == 2): fails['l45'] += 1
        # prior canon: THM-4111 (6) / THM-4115 (7) / Session 210 (8): the same total is (m-1) H(R) + F_1(R),
        # F_1 = number of orderings of V(R) with exactly one backward adjacency
        F1 = sum(1 for p in permutations(range(m)) if sum(1 for x, y in zip(p, p[1:]) if not beats(base, idx, x, y)) == 1)
        if split != (m - 1) * h0 + F1: fails['thm4111'] += 1
        # prior canon: THM-4184 (5) (P^0 = W/2) and Session 210 Lemma 3: half the total is the sum over the non-empty
        # EVEN subsets S of H(R[S]) H(R[V-S]) (H(empty) = 1); with Redei this gives the congruence at once
        evsum = sum(hS[S] * hS[full ^ S] for S in range(1, full + 1) if bin(S).count('1') % 2 == 0)
        if 2 * evsum != split: fails['thm4184'] += 1
        # directed paths on an even number of vertices, by enumeration of sequences
        P = {}
        for k in range(2, m + 1, 2):
            for seq in permutations(range(m), k):
                if all(beats(base, idx, x, y) for x, y in zip(seq, seq[1:])): P[(seq[0], seq[-1])] = P.get((seq[0], seq[-1]), 0) + 1
        if sum(P.values()) % 2 != 1: fails['evenpaths'] += 1
        for (i, j), u in mu.items():
            if (P.get((i, j), 0) + P.get((j, i), 0)) % 2 != (u // 2) % 2: fails['f44'] += 1; break
        n0 = sum(1 for h in Hx if (h - h0) % 4 == 0); n1 = (1 << m) - n0
        if min(n0, n1) < (1 << m) // 4: fails['cor46a'] += 1
        b = Fraction(n0 - n1, 1 << m)
        sum_b2 += b * b; sum_abs += abs(b); cnt += 1
        sum_D += Fraction(sum(1 if h % 4 == 1 else -1 for h in Hx), 1 << m)
    print("m = %d: %d labelled R ; failures %s ; E_R[b_R^2] = %s ; rho_m = %s ; E_R|b_R| = %s ; D_%d = %s"
          % (m, cnt, fails, sum_b2 / cnt, sum_b2 / cnt * 2 ** m, sum_abs / cnt, m + 1, sum_D / cnt))

def part2(n):
    prs, idx = arcs_index(n)
    tA = perm_tests([x for x in range(n) if x != 0], idx)      # T - v,  v = 0
    tB = perm_tests([x for x in range(n) if x != 1], idx)      # T - v', v' = 1
    cyc = []                                                    # odd cycles through exactly one of 0, 1: (care, want)
    for L in range(3, n + 1, 2):
        for S in combinations(range(n), L):
            if (0 in S) == (1 in S): continue
            for rest in permutations(S[1:]):
                c = (S[0],) + rest; care = want = 0
                for a, b in zip(c, c[1:] + (c[0],)):
                    k = idx[(min(a, b), max(a, b))]; care |= 1 << k
                    if a < b: want |= 1 << k
                cyc.append((care, want))
    s1 = s2 = 0; mism = 0
    for code in range(1 << len(prs)):
        ha = sum(1 for care, want in tA if code & care == want); hb = sum(1 for care, want in tB if code & care == want)
        e1 = -1 if ((ha - hb) // 2) % 2 else 1
        e2 = -1 if sum(1 for care, want in cyc if code & care == want) % 2 else 1
        s1 += e1; s2 += e2; mism += (e1 != e2)
    tot = 1 << len(prs)
    print("n = %d (m = %d): E_T (-1)^((H(T-v)-H(T-v'))/2) = %s ; E_T (-1)^(N_1) = %s ; pointwise mismatches %d of %d"
          % (n, n - 2, Fraction(s1, tot), Fraction(s2, tot), mism, tot))

if __name__ == '__main__':
    for m in (2, 3, 4, 5): part1(m)
    for n in (4, 5, 6): part2(n)
