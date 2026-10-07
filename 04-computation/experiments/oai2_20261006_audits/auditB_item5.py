#!/usr/bin/env python3
"""Audit B, item 5 (note sections 5.3-5.5), all by brute force (no SAT):
 (a) Davenport law: max m such that some a_1..a_m >= 1 have ALL nonempty subset sums in one v_p-level
     equals p - 1 (searched over a_i <= p^2 * small, levels v = 0, 1);
 (b) clique number of the graph x ~ y iff v_p(|x-y|) = v on [0, L) equals p;
 (c) least N forcing a monochromatic {a, b, a+b} (a, b 3-smooth, a+b <= N) in every 2-colouring of [N]:
     with a = b allowed, and with a != b (exhaustive over 2^N colourings of the relevant numbers);
 (d) primitive solutions of x + y = z in 3-smooth numbers up to 10^15 (independent enumeration);
     Omega-parity colouring has no monochromatic Schur triple with a, b, a+b in S (a = b allowed);
     the 2-colouring {S, even Omega} vs rest has no monochromatic quartet {a, b, a+b, ab}, a, b in S;
 (e) book Ramsey: red = disjoint cliques, blue = complete multipartite (all part-size vectors),
     no red B_{n-1} and no blue B_n: max N; also the swapped orientation; compare max(2n, 3n-3).
"""
import itertools, math


def vp(x, p):
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


# (a)
def max_m_level(p, v, A):
    """largest m such that a multiset of m elements from A (a_i with v_p = v) has all subset sums at level v"""
    cand = [a for a in A if vp(a, p) == v]
    best = 0
    # residues mod p of a/p^v suffice; search multisets of residues, then realise
    res = sorted({(a // p ** v) % p for a in cand})
    for m in range(1, p + 2):
        found = False
        for combo in itertools.combinations_with_replacement(res, m):
            ok = True
            for r in range(1, m + 1):
                for sub in itertools.combinations(range(m), r):
                    if sum(combo[i] for i in sub) % p == 0:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                found = True
                break
        if found:
            best = m
        else:
            break
    return best

print("(a) Davenport law: max m (all nonempty subset sums in one v_p level)")
for p in (2, 3, 5, 7, 11):
    A = range(1, 4 * p ** 2)
    r0, r1 = max_m_level(p, 0, A), max_m_level(p, 1, A)
    print(f"    p={p}: level v=0 -> {r0}, v=1 -> {r1}  (p-1 = {p-1})")
    assert r0 == r1 == p - 1
# direct integer check for p = 3: brute force over actual integers a_i <= 30, m = 2, 3
for p in (3, 5):
    for m in (p - 1, p):
        ex = None
        for combo in itertools.combinations_with_replacement(range(1, 31), m):
            sums = {sum(combo[i] for i in sub) for r in range(1, m + 1) for sub in itertools.combinations(range(m), r)}
            if len({vp(s, p) for s in sums}) == 1:
                ex = combo
                break
        print(f"    integer check p={p}, m={m}: example with all subset sums in one level: {ex}")

# (b)
print("(b) clique number of v_p(|x-y|) = v")
def clique_number(p, v, L):
    V = list(range(L))
    adj = {x: {y for y in V if y != x and vp(abs(x - y), p) == v} for x in V}
    best = [0]
    def grow(cl, cand):
        if len(cl) > best[0]:
            best[0] = len(cl)
        for i, y in enumerate(cand):
            if len(cl) + len(cand) - i <= best[0]:
                return
            grow(cl + [y], [z for z in cand[i + 1:] if z in adj[y]])
    grow([], V)
    return best[0]
for p in (2, 3, 5):
    for v in (0, 1):
        L = 3 * p ** (v + 1) + 1
        cn = clique_number(p, v, L)
        print(f"    p={p}, v={v}, on [0,{L}): clique number {cn}")
        assert cn == p

# (c)
def smooth(x):
    while x % 2 == 0:
        x //= 2
    while x % 3 == 0:
        x //= 3
    return x == 1

def forced(N, allow_equal):
    S = [a for a in range(1, N + 1) if smooth(a)]
    trip = set()
    for a in S:
        for b in S:
            if (a < b or (allow_equal and a == b)) and a + b <= N:
                trip.add((a, b, a + b))
    nums = sorted({x for tr in trip for x in tr})
    idx = {x: i for i, x in enumerate(nums)}
    tm = [(idx[a], idx[b], idx[c]) for a, b, c in trip]
    for col in range(1 << len(nums)):
        good = True
        for i, j, k in tm:
            ci, cj, ck = (col >> i) & 1, (col >> j) & 1, (col >> k) & 1
            if ci == cj == ck:
                good = False
                break
        if good:
            return False
    return True

print("(c) least N forcing monochromatic {a, b, a+b}, a, b 3-smooth")
for allow in (True, False):
    N = 2
    while not forced(N, allow):
        N += 1
    print(f"    a = b allowed: {allow} -> least N = {N}")

# (d)
LIM = 10 ** 15
S = sorted(2 ** i * 3 ** j for i in range(60) for j in range(40) if 2 ** i * 3 ** j <= LIM)
Sset = set(S)
prim = set()
for a in S:
    for b in S:
        if a <= b and a + b in Sset:
            g = math.gcd(a, b)
            prim.add((a // g, b // g))
print("(d) primitive x + y = z in 3-smooth numbers <= 1e15:", sorted(prim))
def Om(x):
    c = 0
    while x % 2 == 0:
        x //= 2; c += 1
    while x % 3 == 0:
        x //= 3; c += 1
    return c
mono = [(a, b) for a in S for b in S if a <= b and a + b in Sset and Om(a) % 2 == Om(b) % 2 == Om(a + b) % 2]
print("    monochromatic Schur triples (a <= b, all in S) under Omega mod 2:", mono[:5], len(mono))
# quartet colouring: colour 0 = S with even Omega, colour 1 = everything else
def col(x):
    return 0 if (x in Sset and Om(x) % 2 == 0) else 1
Ssmall = [s for s in S if s <= 10 ** 7]
bad = [(a, b) for a in Ssmall for b in Ssmall if a <= b and len({col(a), col(b), col(a + b), col(a * b)}) == 1]
print("    2-colouring (S with even Omega | rest): monochromatic quartets {a,b,a+b,ab}, a,b in S <= 1e7:", len(bad))

# (e)
print("(e) Turan-type colourings for (B_{n-1}, B_n)")
def partitions_bounded(total_max, kmax, smax):
    # all non-increasing tuples of positive part sizes, k <= kmax, parts <= smax
    def rec(prefix, maxpart, remaining_k):
        yield prefix
        if remaining_k == 0:
            return
        for s in range(min(maxpart, smax), 0, -1):
            yield from rec(prefix + (s,), s, remaining_k - 1)
    yield from rec((), smax, kmax)

def best_N(n, red_book, blue_book):
    """red = disjoint cliques of sizes s_i; blue = complete multipartite on these parts.
    red B_r present iff some clique has s >= r + 2 (edge + r common neighbours);
    blue B_r present iff some blue edge (parts i != j) has >= r common blue neighbours: N - s_i - s_j >= r."""
    best = (0, None)
    for parts in partitions_bounded(None, 6, red_book + 1):
        if not parts:
            continue
        N = sum(parts)
        if any(s >= red_book + 2 for s in parts):
            continue
        k = len(parts)
        if k >= 2:
            srt = sorted(parts)
            if N - srt[0] - srt[1] >= blue_book:  # worst pair = two smallest parts
                continue
        if N > best[0]:
            best = (N, parts)
    return best

for n in list(range(3, 13)) + [16, 20]:
    b1 = best_N(n, n - 1, n)          # red avoids B_{n-1}, blue avoids B_n
    b2 = best_N(n, n, n - 1)          # swapped: red avoids B_n ... i.e. blue cliques avoid B_n, red multipartite avoids B_{n-1}
    print(f"    n={n}: red cliques/blue multipartite max N = {b1[0]} {b1[1]};  max(2n,3n-3) = {max(2*n, 3*n-3)};"
          f"  swapped orientation max N = {b2[0]} {b2[1]}")
