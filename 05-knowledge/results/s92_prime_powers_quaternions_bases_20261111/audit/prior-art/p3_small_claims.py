#!/usr/bin/env python3
"""Auditor prior-art (S92).  Independent recomputation of the small statements the session types NEW or "added"
(note sections 3, 6, 7, 8; THM-5135 (5); THM-5134 (H)), so that the credit findings rest on checked statements.
Own code; nothing imported from the session.

 S1  five known cycles of C(n) = n/2 | 3n+1 on Z: arcs, labels |4n+1|, which labels are path counts.
 S2  orders 4, 5, 6: values of H, signed labels n = (chi_4(H) H - 1)/4; order 4: self-converse iff H = 1 mod 4 iff
     n >= 0 iff c_3 even; strong classes of order 5: H, |Aut|, self-converse; counts of classes by order.
 S3  arc-usage profiles at order 4; classes by strong components: 1/(1 - (x + x^3 + x^4 + 6x^5)).
 S4  fan primes: p < 200 with ord_p(2) even; first p-divisible state against the fan 2^(t-1).
 S5  eight central extensions of F_2^2 by F_2: 16 cocycles, classes, groups, commutativity, triangle types.
 S6  r_2 convolutions: additive = r_4 = 8 sigma(n) (n odd), multiplicative = 0 for n = 3 mod 4.
 S7  base-b one-apex formula against THM-5070 (2); balanced states X = y in base 3 below 2^24.
"""
import itertools
import sys
from math import gcd

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
ok = True


def check(name, cond, detail=""):
    global ok
    ok = ok and bool(cond)
    print("%-5s %s  %s" % (name, "PASS" if cond else "FAIL", detail))


def psi(h):
    if h == 0:
        return 0
    if h % 2 == 0:
        return psi(h // 2) + h // 2
    k = h // 2
    return psi(k) + k // 2 + 1


def X(h):
    return 2 * psi(h) + 1


def Y(h):
    return 2 * h + 1


# ---------------------------------------------------------------- tournaments, small orders
def hp_count(n, adj):
    # adj[u] = bitmask of out-neighbours; number of Hamiltonian paths by subset DP
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if c:
                out = adj[v] & ~mask
                while out:
                    w = (out & -out).bit_length() - 1
                    dp[mask | (1 << w)][w] += c
                    out &= out - 1
    return sum(dp[full])


def tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        adj = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
        yield adj


def canon(n, adj):
    best = None
    for perm in itertools.permutations(range(n)):
        key = 0
        for i in range(n):
            a = adj[perm[i]]
            for j in range(n):
                if a >> perm[j] & 1:
                    key |= 1 << (i * n + j)
        if best is None or key < best:
            best = key
    return best


def converse(n, adj):
    return [((~adj[v]) & ((1 << n) - 1)) & ~(1 << v) for v in range(n)]


def aut_order(n, adj):
    cnt = 0
    for perm in itertools.permutations(range(n)):
        good = True
        for i in range(n):
            a = adj[i]
            for j in range(n):
                if i != j and ((a >> j) & 1) != ((adj[perm[i]] >> perm[j]) & 1):
                    good = False
                    break
            if not good:
                break
        cnt += good
    return cnt


def is_strong(n, adj):
    def reach(start, A):
        seen = 1 << start
        stack = [start]
        while stack:
            v = stack.pop()
            nb = A[v] & ~seen
            while nb:
                w = (nb & -nb).bit_length() - 1
                seen |= 1 << w
                stack.append(w)
                nb &= nb - 1
        return seen
    full = (1 << n) - 1
    return reach(0, adj) == full and reach(0, converse(n, adj)) == full


def c3_count(n, adj):
    c = 0
    for a, b, d in itertools.combinations(range(n), 3):
        if (adj[a] >> b & 1 and adj[b] >> d & 1 and adj[d] >> a & 1) or (adj[a] >> d & 1 and adj[d] >> b & 1 and adj[b] >> a & 1):
            c += 1
    return c


def signed_label(H):
    chi = 1 if H % 4 == 1 else -1
    return (chi * H - 1) // 4


print("=" * 30, "S2/S3: orders 3..6")
classes = {}
for n in (1, 2, 3, 4, 5):
    cl = {}
    for adj in tournaments(n):
        k = canon(n, adj)
        if k not in cl:
            cl[k] = list(adj)
    classes[n] = cl
check("S3a", [len(classes[n]) for n in (1, 2, 3, 4, 5)] == [1, 1, 2, 4, 12], "classes at orders 1..5: %s" % [len(classes[n]) for n in (1, 2, 3, 4, 5)])
strong_counts = [sum(1 for adj in classes[n].values() if is_strong(n, adj)) for n in (1, 2, 3, 4, 5)]
check("S3b", strong_counts == [1, 0, 1, 1, 6], "strong classes at orders 1..5: %s" % strong_counts)
# series 1/(1 - S): coefficients
S = [0] + strong_counts
T = [1] + [0] * 5
for n in range(1, 6):
    T[n] = sum(S[k] * T[n - k] for k in range(1, n + 1))
check("S3c", T[1:] == [1, 1, 2, 4, 12], "1/(1 - (x + x^3 + x^4 + 6x^5)) = 1 + %s (orders 1..5)" % T[1:])

print("order 4 classes: H, c3, |Aut|, self-converse, strong, signed label n, arc-usage profile (non-zero)")
rows4 = []
for k, adj in classes[4].items():
    n = 4
    H = hp_count(n, adj)
    sc = canon(n, converse(n, adj)) == k
    # arc usage
    use = {}
    for perm in itertools.permutations(range(n)):
        if all(adj[perm[i]] >> perm[i + 1] & 1 for i in range(n - 1)):
            for i in range(n - 1):
                use[(perm[i], perm[i + 1])] = use.get((perm[i], perm[i + 1]), 0) + 1
    prof = sorted(use.values(), reverse=True)
    rows4.append((H, c3_count(n, adj), aut_order(n, adj), sc, is_strong(n, adj), signed_label(H), prof))
for r in sorted(rows4):
    print("     H = %d  c3 = %d  |Aut| = %d  self-converse = %s  strong = %s  n = %+d  profile = %s" % r)
profs = sorted(tuple(r[6]) for r in rows4)
check("S3d", profs == sorted([(1, 1, 1), (2, 2, 2, 1, 1, 1), (2, 2, 2, 1, 1, 1), (4, 3, 3, 3, 1, 1)]), "arc-usage profiles at order 4")
check("S2a", all((r[3]) == (r[0] % 4 == 1) == (r[5] >= 0) == (r[1] % 2 == 0) for r in rows4),
      "order 4: self-converse iff H = 1 mod 4 iff n >= 0 iff c3 even")
check("S2b", sorted(r[5] for r in rows4) == [-1, -1, 0, 1], "order 4 labels n: %s" % sorted(r[5] for r in rows4))

print("order 5 classes:")
rows5 = []
for k, adj in classes[5].items():
    n = 5
    H = hp_count(n, adj)
    sc = canon(n, converse(n, adj)) == k
    rows5.append((is_strong(n, adj), H, aut_order(n, adj), sc))
for r in sorted(rows5):
    print("     strong = %-5s H = %2d  |Aut| = %d  self-converse = %s" % r)
st5 = sorted((r[1], r[2], r[3]) for r in rows5 if r[0])
check("S2c", [r[0] for r in st5] == [9, 9, 11, 13, 15, 15], "strong order-5 H values: %s" % [r[0] for r in st5])
check("S2d", sorted(r[1] for r in st5) == [1, 1, 1, 1, 3, 5], "strong order-5 |Aut|: %s" % sorted(r[1] for r in st5))
check("S2e", all(r[2] for r in st5), "all six strong classes of order 5 are self-converse")
nsc = sum(1 for r in rows5 if r[3])
check("S2f", nsc == 8 and all(r[3] for r in rows5 if r[0]) and sum(1 for r in rows5 if not r[3] and not r[0]) == 4,
      "order 5: %d self-converse classes of 12; the 4 others are not strong" % nsc)
H5 = sorted(set(r[1] for r in rows5))
lab5 = sorted(signed_label(H) for H in H5)
check("S2g", lab5 == [-4, -3, -1, 0, 1, 2, 3], "order 5: H values %s, labels %s = [-4, 3] without -2" % (H5, lab5))
H6 = set()
for adj in tournaments(6):
    H6.add(hp_count(6, adj))
lab6 = sorted(signed_label(H) for H in H6)
want6 = [v for v in range(-11, 12) if v not in (-10, -9, -2, 5)]
check("S2h", lab6 == want6, "order 6: H values %s; labels = [-11, 11] without %s" % (sorted(H6), sorted(set(range(-11, 12)) - set(lab6))))

print("=" * 30, "S1: five known cycles of 3n+1 on Z")


def C(n):
    return n // 2 if n % 2 == 0 else 3 * n + 1


cycles = []
for start in (0, 1, -1, -5, -17):
    cyc = [start]
    v = C(start)
    while v != start:
        cyc.append(v)
        v = C(v)
    cycles.append(cyc)
arcs = [(v, C(v)) for cyc in cycles for v in cyc]
labels = [abs(4 * v + 1) for v, _ in arcs]
print("     cycle lengths:", [len(c) for c in cycles], " arcs:", len(arcs))
for cyc in cycles:
    print("     cycle of %4d: labels %s" % (cyc[0], [abs(4 * v + 1) for v in cyc]))
check("S1a", len(arcs) == 29 and len(set(labels)) == 29, "29 arcs, 29 distinct labels")
# which labels are path counts: build a set of attained odd values <= 1200 from one-apex counts, F1 joins and products
LIM = 1200
att = {1}
for h in range(0, 700):
    if X(h) <= LIM:
        att.add(X(h))
states = [(X(h), Y(h)) for h in range(0, 700)]
for (x1, y1) in states:
    if x1 > LIM:
        continue
    for (x2, y2) in states:
        if x1 * x2 > LIM:
            continue
        s = 0
        while True:
            N = x1 * x2 + (1 << s) * y1 * y2 - (1 if s == 0 else 0)
            if N > LIM:
                break
            att.add(N)
            s += 1
changed = True
while changed:
    changed = False
    cur = sorted(v for v in att if v > 1)
    for a in cur:
        for b in cur:
            if a * b > LIM:
                break
            if a * b not in att:
                att.add(a * b)
                changed = True
miss_all = [v for v in range(1, LIM, 2) if v not in att]
print("     odd values <= %d not reached by one-apex counts, F1 joins and their products: %d values, first %s"
      % (LIM, len(miss_all), miss_all[:12]))
print("     (F1 alone is not the whole spectrum; THM-4875 (D) uses F1 and the bank family.)")
# second, larger family: every state of the O/Z/J grammar (THM-4888 setting) from the rooted edge is a path count
att2 = set()
layer = {(1, 0, 2, 1)}
seen_states = set(layer)
for L in range(1, 16):
    new = set()
    for (D, E, U, d) in layer:
        for st in ((D, 0, 2 * U, d), (D - E + U, d, 2 * U - d, d), (D - E + U, D - E, D + U, D)):
            if st[0] <= LIM and st not in seen_states and st[2] <= 8 * LIM:
                new.add(st)
    seen_states |= new
    layer = new
att2 = {st[0] for st in seen_states}
changed = True
while changed:
    changed = False
    cur = sorted(v for v in att2 if v > 1)
    for a in cur:
        for b in cur:
            if a * b > LIM:
                break
            if a * b not in att2:
                att2.add(a * b)
                changed = True
miss2 = [v for v in range(1, LIM, 2) if v not in att2]
print("     odd values <= %d not reached by O/Z/J grammar states (words of length <= 15, U <= %d) and products: %s"
      % (LIM, 8 * LIM, miss2))
att = att | att2
notatt = [l for l in labels if l not in att]
print("     labels of the 29 arcs not in the union of these two families: %s" % notatt)
print("     (these two families are not the whole spectrum; explicit tournaments for every label except 7 are in")
print("      p4_order7_spectrum.out: order-7 exhaustive spectrum, witnesses of orders 8 and 9 for the others)")
check("S1b", 7 in labels and set(notatt) <= {7, 67, 163}, "7 is a label (arc -2 -> -1); the other labels are path counts by p4_order7_spectrum.out")
check("S1c", 21 not in labels and abs(4 * 5 + 1) == 21, "21 = label of the arc 5 -> 16, not on a known cycle")
check("S1d", [abs(4 * ((4 ** k - 1) // 3) + 1) for k in (1, 2, 3)] == [5, 21, 85], "trunk arcs (4^k - 1)/3 -> 4^k have labels 5, 21, 85")

print("=" * 30, "S4: fan primes")


def ordp2(p):
    k, v = 1, 2 % p
    while v != 1:
        v = v * 2 % p
        k += 1
    return k


primes = [p for p in range(3, 200) if all(p % d for d in range(2, int(p ** 0.5) + 1))]
first_div = {}
HMAX = 1 << 20
need = set(primes)
h = 1
xs = None
# psi by doubling, as python ints (2^20 states)
ps = [0] * HMAX
for hh in range(1, HMAX):
    k = hh >> 1
    ps[hh] = ps[k] + k if hh % 2 == 0 else ps[k] + (k >> 1) + 1
for hh in range(1, HMAX):
    g = gcd(2 * ps[hh] + 1, 2 * hh + 1)
    if g > 1:
        for p in list(need):
            if g % p == 0:
                first_div[p] = hh
                need.discard(p)
        if not need:
            break
fanfirst = []
evenord = []
for p in primes:
    o = ordp2(p)
    if o % 2 == 0:
        t = o // 2
        evenord.append(p)
        fan = 1 << (t - 1)
        assert (2 ** t + 1) % p == 0
        if first_div.get(p) == fan:
            fanfirst.append(p)
print("     primes < 200 with ord_p(2) even (fan primes = alternation primes):", evenord)
print("     first p-divisible state (h < 2^20):", {p: first_div.get(p) for p in primes[:25]})
check("S4a", fanfirst == [3, 5, 11, 13, 17, 43, 113], "the fan 2^(t-1) is the first p-divisible state exactly for %s (p < 200)" % fanfirst)
# balanced p-divisible state iff ord even: x = y only at h = 0 and 2^j (check h < 2^20)
bal = [hh for hh in range(0, HMAX) if 2 * ps[hh] + 1 == 2 * hh + 1]
check("S4b", bal == [0] + [1 << j for j in range(20)], "x = y exactly at h = 0 and h = 2^j (h < 2^20)")

print("=" * 30, "S5: central extensions of F_2^2 by F_2")
V = [0, 1, 2, 3]  # xor addition
nz = [1, 2, 3]
cocycles = []
for bits in range(1 << 9):
    f = {}
    k = 0
    for a in V:
        for b in V:
            if a == 0 or b == 0:
                f[(a, b)] = 0
            else:
                f[(a, b)] = bits >> k & 1
                k += 1
    good = all((f[(a, b)] + f[(a ^ b, c)] - f[(b, c)] - f[(a, b ^ c)]) % 2 == 0 for a in V for b in V for c in V)
    if good:
        cocycles.append(f)
check("S5a", len(cocycles) == 16, "normalized 2-cocycles: %d" % len(cocycles))
# classes modulo coboundaries of normalized 1-cochains
cob = set()
for g1, g2, g3 in itertools.product((0, 1), repeat=3):
    g = {0: 0, 1: g1, 2: g2, 3: g3}
    cob.add(tuple((g[a] + g[b] + g[a ^ b]) % 2 for a in V for b in V))
check("S5b", len(cob) == 2, "distinct normalized coboundaries: %d" % len(cob))
cls = {}
for f in cocycles:
    vec = tuple(f[(a, b)] for a in V for b in V)
    reps = sorted(tuple((v + c) % 2 for v, c in zip(vec, cb)) for cb in cob)
    cls.setdefault(reps[0], []).append(f)
check("S5c", len(cls) == 8 and all(len(v) == 2 for v in cls.values()), "classes: %d, each of 2 cocycles" % len(cls))
table = {}
for rep, fs in cls.items():
    f = fs[0]
    Sset = tuple(a for a in nz if f[(a, a)] == 1)
    comm = all(f[(a, b)] == f[(b, a)] for a in V for b in V)
    # group on V x F_2: count elements of order 4
    def mul(u, v):
        return (u[0] ^ v[0], (u[1] + v[1] + f[(u[0], v[0])]) % 2)
    elems = [(a, z) for a in V for z in (0, 1)]
    ord4 = sum(1 for u in elems if mul(u, u) != (0, 0))
    grp = {0: "C2^3", 4: "C4xC2", 2: "D4", 6: "Q8"}[ord4]
    if not comm and grp == "C4xC2":
        grp = "??"
    # sign table on the three points: a -> b iff e_a e_b = + e_(a+b)
    tri = None
    mid = None
    if not comm:
        arcs3 = [(a, b) for a in nz for b in nz if a != b and f[(a, b)] == 0]
        outdeg = {a: sum(1 for (u, v) in arcs3 if u == a) for a in nz}
        tri = "cyclic" if sorted(outdeg.values()) == [1, 1, 1] else "transitive"
        if tri == "transitive":
            mid = [a for a in nz if outdeg[a] == 1][0]
    else:
        # both cocycles of a commutative class are symmetric: signed triangles
        tri = "signed triangle"
    # invariance of S within the class
    assert all(tuple(a for a in nz if g[(a, a)] == 1) == Sset for g in fs)
    table[Sset] = (grp, comm, tri, mid)
for Sset in sorted(table, key=lambda s: (len(s), s)):
    grp, comm, tri, mid = table[Sset]
    print("     S = %-10s group %-6s commutative %-5s table: %s%s" % (Sset, grp, comm, tri, "" if mid is None else ", middle vertex %d" % mid))
check("S5d", len(table) == 8, "the class is determined by S and every subset of the three points occurs")
check("S5e", all((len(S) % 2 == 0) == v[1] for S, v in table.items()), "commutative iff #S even")
check("S5f", all(v[0] == {0: "C2^3", 1: "D4", 2: "C4xC2", 3: "Q8"}[len(S)] for S, v in table.items()), "groups C2^3, D4, C4xC2, Q8 for #S = 0, 1, 2, 3")
check("S5g", all((v[2] == "cyclic") == (len(S) == 3) for S, v in table.items() if not v[1])
      and all(v[3] == S[0] for S, v in table.items() if len(S) == 1), "non-commutative: cyclic iff #S = 3; transitive with the unit of square -1 in the middle iff #S = 1")

print("=" * 30, "S6: r_2 convolutions")
NMAX = 400
r2 = [0] * (NMAX + 1)
B = int(NMAX ** 0.5) + 1
for a in range(-B, B + 1):
    for b in range(-B, B + 1):
        if a * a + b * b <= NMAX:
            r2[a * a + b * b] += 1


def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


addok = all(sum(r2[a] * r2[n - a] for a in range(0, n + 1)) == 8 * sigma(n) for n in range(1, NMAX + 1, 2))
mulok = all(sum(r2[a] * r2[n // a] for a in range(1, n + 1) if n % a == 0) == 0 for n in range(3, NMAX + 1, 4))
mulnz = [n for n in range(1, 60) if sum(r2[a] * r2[n // a] for a in range(1, n + 1) if n % a == 0) > 0]
check("S6a", addok, "sum_{a+b=n} r_2(a) r_2(b) = 8 sigma(n) for odd n <= 400")
check("S6b", mulok, "sum_{ab=n} r_2(a) r_2(b) = 0 for n = 3 mod 4, n <= 400; non-zero n < 60: %s" % mulnz)

print("=" * 30, "S7: base b")


def a_kp(W):  # THM-5070 (2): a(2W+1) = a(W), a(2W) = a(W) + ceil(W/2), a(0) = 0
    if W == 0:
        return 0
    if W % 2 == 1:
        return a_kp(W // 2)
    return a_kp(W // 2) + (W // 2 + 1) // 2


okb = True
for h1 in range(0, 40):
    for h2 in range(0, 40):
        for j1 in range(0, 3):
            for j2 in range(0, 3):
                W = (2 * h1) * (1 << j1) + (1 << j1) - 1
                Wp = (2 * h2) * (1 << j2) + (1 << j2) - 1
                s = j1 + j2
                for b in (2, 3, 4, 5, 6, 10):
                    lhs = (1 + b * a_kp(W)) * (1 + b * a_kp(Wp)) + b * (((W + 1) // 2) * (1 + Wp // 2) + (1 + W // 2) * ((Wp + 1) // 2))
                    cross = (1 << s) * Y(h1) * Y(h2) - (1 if s == 0 else 0)
                    assert cross % 2 == 0
                    rhs = (1 + b * psi(h1)) * (1 + b * psi(h2)) + b * cross // 2
                    if lhs != rhs:
                        okb = False
check("S7a", okb, "THM-5070 (2) with W = 2h 2^j + 2^j - 1 equals (1 + b psi1)(1 + b psi2) + b (2^s y1 y2 - [s=0])/2, s = j1 + j2")
N24 = 1 << 24
bal3 = []
import numpy as np
psarr = np.zeros(N24, dtype=np.int64)
lo = 1
while lo < N24:
    hi = min(2 * lo, N24)
    hh = np.arange(lo, hi, dtype=np.int64)
    k = hh >> 1
    psarr[lo:hi] = np.where(hh % 2 == 0, psarr[k] + k, psarr[k] + (k >> 1) + 1)
    lo = hi
hall = np.arange(N24, dtype=np.int64)
bal3 = hall[(1 + 3 * psarr) == (2 * hall + 1)].tolist()
check("S7b", len(bal3) == 12 and bal3[:6] == [0, 3, 45, 117, 309, 9045], "base 3 balanced states X = y below 2^24: %s" % bal3)
bal4 = int(((1 + 4 * psarr) == (2 * hall + 1)).sum())
minratio = float(((psarr[1:] * 2) > hall[1:]).all())
check("S7c", bal4 == 1 and minratio == 1.0, "base 4: only h = 0; psi(h) > h/2 for 1 <= h < 2^24")
print("ALL PASS" if ok else "SOME FAIL")
