"""Independent checks of the combinatorics in THM-4819 Lemma B (odd) and Lemma B' (even), and of the cycle-type remark.

(1) Odd case, Z_{2m}, |U| = m.  (a): t in U => t+2 not in U.  (b): for distinct t, t' in U exactly one of 2t-t'+1, 2t'-t+1
    lies in U.  Enumerates ALL m-subsets satisfying (a) (and then (b)), without normalising 0 in U, for m = 2..11.
    Expected from the proof: no (a)-set for odd m >= 3; for even m exactly the four sets (a+4Z) u (b+4Z), a even, b odd;
    none of them satisfies (b) for m >= 4.
(2) Even case, Z_{2m-1}, |U| = m-1, same (a), (b).  Expected: the (a)-sets are exactly the 2m-1 translates of
    2*{0,2,...,2m-4}; (b) has solutions for m = 3, 4 only (lists them with 0 in U); for m = 5..60 the progression has a
    violating pair, and for m >= 7 the pair (0, 4) itself violates, with halved witnesses m-2, m+4.
(3) The remark's cycle types of sigma_t sigma_t' (d = t - t', e = gcd(d, m), l = m/e), by direct composition, m = 3..40,
    and how many gaps are 'settled' if one accepts or rejects the remark's claim that an n-cycle is impossible.
"""
from itertools import combinations
from math import gcd

def sat_a(U, M):
    S = set(U); return all((t + 2) % M not in S for t in U)
def sat_b(U, M):
    S = set(U)
    return all(((2*t - u + 1) % M in S) + ((2*u - t + 1) % M in S) == 1 for t, u in combinations(U, 2))

def a_sets(M, k):
    """all k-subsets of Z_M with no two elements differing by 2 (backtracking on the two +2-orbits is overkill; plain)"""
    res = []
    def rec(start, cur):
        if len(cur) == k: res.append(tuple(cur)); return
        for x in range(start, M):
            if M - x < k - len(cur): break
            if any((x - y) % M in (2, M - 2) for y in cur): continue
            cur.append(x); rec(x + 1, cur); cur.pop()
    rec(0, [])
    return res

print("(1) odd case: Z_{2m}, |U| = m")
for m in range(2, 12):
    M = 2 * m; A = a_sets(M, m); Bs = [U for U in A if sat_b(U, M)]
    expected = None
    if m % 2 == 0:
        expected = sorted(tuple(sorted({(a + 4*i) % M for i in range(m // 2)} | {(b + 4*i) % M for i in range(m // 2)}))
                          for a in (0, 2) for b in (1, 3))
    print(f"   m={m}: (a)-sets {len(A)}" + (f" (= the four mod-4 unions: {sorted(A) == expected})" if expected else "") +
          f"; (a)+(b)-sets {len(Bs)} {Bs[:4]}")

print("(2) even case: Z_{2m-1}, |U| = m-1")
for m in range(3, 13):
    M = 2 * m - 1; A = a_sets(M, m - 1)
    ap = sorted({(2 * (2 * i + s)) % M for i in range(m - 1)} for s in range(M)) if False else None
    translates = sorted(tuple(sorted(((2 * x) + 2 * s) % M for x in range(0, 2 * m - 3, 2))) for s in range(M))
    Bs = [U for U in A if sat_b(U, M)]
    print(f"   m={m}: (a)-sets {len(A)}; equal to the {M} translates of 2*{{0,2,...,2m-4}}: {sorted(A) == sorted(set(translates))}; "
          f"(a)+(b)-sets {len(Bs)}; with 0 in U: {[U for U in Bs if 0 in U]}")
bad = []
for m in range(5, 61):
    M = 2 * m - 1; U = sorted({(2 * x) % M for x in range(0, 2 * m - 3, 2)})
    S = set(U); viol = [(t, u) for t, u in combinations(U, 2) if ((2*t - u + 1) % M in S) + ((2*u - t + 1) % M in S) != 1]
    w1, w2 = (-3) % M, 9 % M; h1, h2 = (w1 * m) % M, (w2 * m) % M
    pair04 = ((w1 in S) + (w2 in S)) != 1
    if not viol or (m >= 7 and (not pair04 or (h1, h2) != (m - 2, (m + 4) % M))): bad.append(m)
    if m <= 8: print(f"   m={m}: progression U = {U}; violating pairs {len(viol)} (first {viol[:3]}); pair (0,4) violates: {pair04}; halved witnesses {h1},{h2}")
print(f"   m = 5..60: progression always has a violating pair, and (0,4) violates with halved witnesses m-2, m+4 for m >= 7: {not bad} {bad}")

print("(3) cycle types of sigma_t sigma_t' (remark 'an alternative route that fails')")
def sigma(t, M):
    p = {'v': t, t: 'v'}
    for x in range(M):
        if x != t: p[x] = (2 * t - x) % M
    return p
def ctype(p):
    seen = set(); L = []
    for s in p:
        if s in seen: continue
        k = 0; x = s
        while x not in seen: seen.add(x); x = p[x]; k += 1
        L.append(k)
    return sorted(L, reverse=True)
def formula(m, d):
    e = gcd(d, m); l = m // e
    return sorted(([1 + 2 * l] + [l] * (2 * e - 2)) if (d // e) % 2 else ([(l + 1) // 2] * 2 + [l] * (2 * e - 1)), reverse=True)
ok = True
for m in range(3, 41):
    M = 2 * m
    for tp in range(M):
        for t in range(M):
            if t == tp: continue
            s_t, s_tp = sigma(t, M), sigma(tp, M)
            g = {x: s_t[s_tp[x]] for x in s_tp}          # apply sigma_t' first
            if ctype(g) != formula(m, (t - tp) % M): ok = False
print(f"   direct composition = closed form for all m = 3..40 and all pairs: {ok}")
def lcm_parity_even(L):
    return any(x % 2 == 0 for x in L)
def clique(M, D, size):
    S = set(D)
    def ext(cur, cands):
        if len(cur) == size: return cur
        for i, c in enumerate(cands):
            r = ext(cur + [c], [d for d in cands[i + 1:] if (d - c) % M in S])
            if r: return r
        return None
    return ext([0], [x for x in range(1, M) if x in S])
uns_with, uns_without = [], []
for m in range(3, 41):
    M = 2 * m
    bad_with = [d for d in range(1, M) if not (lcm_parity_even(formula(m, d)) or formula(m, d) == [2 * m + 1])]
    bad_without = [d for d in range(1, M) if not lcm_parity_even(formula(m, d))]
    if clique(M, bad_with, m): uns_with.append(m)
    if m <= 24 and clique(M, bad_without, m): uns_without.append(m)
print(f"   m <= 40 NOT settled if the n-cycle case counts as settled: {uns_with}")
print(f"   m <= 24 NOT settled if only 'even order is impossible' is used: {uns_without}")
