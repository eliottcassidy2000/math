#!/usr/bin/env python3
"""ADVERSARIAL AUDIT of lane berggren_edge_transport (session collatz-mod6-20260917, wave 2026-09-21).

Independent recomputation of every key number of
  04-computation/experiments/collatz_mod6_20260921_berggren_edge_transport.py
with different code paths: Berggren's triple matrices instead of Euclid (m,n) algebra (A1),
direct edge enumeration instead of the numpy pair census (A8), Euclid (m,n) enumeration for the
hypotenuse count (A9), an exact integer linear solve of the S11 'proof sketch' (A11), a component
analysis of the legal sub-forest (A6), and a hunt for the accidental exceptions that the note's
section 4 wording ('never a (3,+-1) edge') overlooked (A7).  Explicit raise only; RAM < 200 MB;
runtime < 1 min.

Reproduce:
  python3 04-computation/experiments/collatz_mod6_20260921_berggren_edge_transport_audit.py \
      > 05-knowledge/results/collatz_mod6_20260921_berggren_edge_transport_audit.out
"""
import math
from collections import defaultdict


def check(cond, msg):
    if not cond:
        raise RuntimeError("AUDIT CHECK FAILED: " + msg)


def v2(n):
    return (n & -n).bit_length() - 1


def is_pow2(q):
    return q >= 1 and (q & (q - 1)) == 0


def read(p, q, a, b):
    """k>=1 with a*p+b = 2^k q, else None."""
    v = a * p + b
    if v <= 0 or v % q:
        return None
    r = v // q
    return v2(r) if (r >= 2 and is_pow2(r)) else None


def legal31(pair):
    """readings of the unordered pair as a (3,+-1) edge, either orientation: list of (src,tgt,b,k)."""
    s, t = pair
    out = []
    for b in (1, -1):
        for (p, q) in ((s, t), (t, s)):
            k = read(p, q, 3, b)
            if k is not None:
                out.append((p, q, b, k))
    return out


def kids(s, t):
    return ((s + 2 * t, t), (2 * s + t, s), (2 * s - t, s))


def edges31(lim):
    """all (3,+-1) edges x->y with max(x,y)<=lim, x!=y: list of (x,y,b,k)."""
    out = []
    for b in (1, -1):
        for x in range(1, 4 * lim + 1, 2):
            v = 3 * x + b
            if v <= 0:
                continue
            k = v2(v)
            y = v >> k
            if y == x or max(x, y) > lim:
                continue
            out.append((x, y, b, k))
    return out


print("=" * 78)
print("AUDIT of lane berggren_edge_transport (independent recomputation)")
print("=" * 78)

# ------------------------------------------------------------------ A1
print("\n### A1  root-coordinate Berggren maps via the TRIPLE matrices (not Euclid (m,n))")
# Berggren matrices on (A odd leg, B even leg, C): standard form.
BERG = (((1, -2, 2), (2, -1, 2), (2, -2, 3)),
        ((1, 2, 2), (2, 1, 2), (2, 2, 3)),
        ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3)))


def isqrt_exact(n):
    r = math.isqrt(n)
    check(r * r == n, "odd square root")
    return r


n1 = 0
matched = defaultdict(int)
for s in range(3, 202, 2):
    for t in range(1, s, 2):
        if math.gcd(s, t) != 1:
            continue
        A, B, C = s * t, (s * s - t * t) // 2, (s * s + t * t) // 2
        check(A * A + B * B == C * C and math.gcd(A, B) == 1, "triple")
        got = set()
        for Mx in BERG:
            A2, B2, C2 = (Mx[i][0] * A + Mx[i][1] * B + Mx[i][2] * C for i in range(3))
            check(A2 * A2 + B2 * B2 == C2 * C2 and A2 > 0 and B2 > 0, "child triple")
            if A2 % 2 == 0:  # make the even leg B2
                A2, B2 = B2, A2
            got.add((isqrt_exact(C2 + B2), isqrt_exact(C2 - B2)))
        want = set(kids(s, t))
        check(got == want, "root form of the triple matrices at %s" % ((s, t),))
        n1 += 1
print("A1 root pairs (s,t) coprime odd, s<=201: %d ; triple-matrix children == {(s+2t,t),(2s+t,s),(2s-t,s)} on all" % n1)
print("    S1: CONFIRMED (independent route through the (A,B,C) matrices).")

# ------------------------------------------------------------------ A2
print("\n### A2  transport table: solve the child's guard FROM the child pair (a<=31, |b|<=31, x<1000)")
n2 = 0
bad = 0
for a in range(1, 32, 2):
    for b in range(-31, 32, 2):
        for x in range(1, 1000, 2):
            v = a * x + b
            if v <= 0:
                continue
            k = v2(v)
            y = v >> k
            if k == 0 or y == x or math.gcd(x, y) != 1:
                continue
            s, t = max(x, y), min(x, y)
            for (p, q) in kids(s, t):
                n2 += 1
                if x in (p, q):  # x-type: solve a' from a' x + b' = 2^k u, b' in {b,-b}, compare with (2)
                    u = q if p == x else p
                    sols = [(bb, (2 ** k * u - bb) // x) for bb in (b, -b) if (2 ** k * u - bb) % x == 0]
                    pred = (2 ** (k + 1) - a, -b) if u == 2 * x - y else (a + 2 ** (k + 1), b)
                    if (pred[1], pred[0]) not in sols:
                        bad += 1
                else:  # y-type: a u + alpha b = M y with u = x+2y (alpha=+1) or 2y-x (alpha=-1); solve M
                    check(y in (p, q), "one endpoint retained")
                    u = q if p == y else p
                    check(u in (x + 2 * y, 2 * y - x), "y-type child form")
                    alpha = +1 if u == x + 2 * y else -1
                    M = a * u + alpha * b
                    if M % y or M // y != alpha * 2 ** k + 2 * a:
                        bad += 1
print("A2 children examined: %d ; guard mismatches: %d" % (n2, bad))
check(bad == 0, "transport table")
print("    S2: CONFIRMED (explorer: 20233 edges a<=15; here a<=31, |b|<=31, wider).")

# ------------------------------------------------------------------ A3
print("\n### A3  multiplier-keeping classification (3) for all odd a<=129, k<=14")
bad = 0
rows = []
for a in range(1, 130, 2):
    for k in range(1, 15):
        for alpha in (1, -1):
            M = alpha * 2 ** k + 2 * a
            keep = M >= 2 and is_pow2(M)
            if k == 1:
                formula = (alpha == 1 and is_pow2(a + 1)) or (alpha == -1 and a >= 3 and is_pow2(a - 1))
            else:
                formula = (alpha == -1 and a == 2 ** (k - 1) + 1)
            if keep != formula:
                bad += 1
            if keep and a <= 9:
                rows.append((a, k, alpha, M))
print("A3 (a,k,alpha,M) with M a power of two, a<=9: %s" % rows)
print("    formula (3) mismatches over a<=129, k<=14: %d" % bad)
check(bad == 0, "classification (3)")
print("    S3: CONFIRMED.")

# ------------------------------------------------------------------ A4
print("\n### A4  orientation = valuation for (3,+-1); smallest |b| witnesses for the two 'possible only if' rows")
viol = 0
for b in (1, -1):
    for x in range(1, 200001, 2):
        v = 3 * x + b
        if v <= 0:
            continue
        k = v2(v)
        y = v >> k
        if y != x and (y > x) != (k == 1):
            viol += 1
print("A4 violations of (y>x <=> k=1), b=+-1, x<=2*10^5: %d" % viol)
check(viol == 0, "orientation")
wit1 = None  # k=1 with y<x
wit2 = None  # k=2 with y>x
for absb in range(1, 40, 2):
    for b in (-absb, absb):
        for x in range(1, 200, 2):
            v = 3 * x + b
            if v <= 0:
                continue
            k = v2(v)
            y = v >> k
            if y == x or k == 0 or math.gcd(x, y) != 1:
                continue
            if k == 1 and y < x and wit1 is None:
                wit1 = (x, y, b, k)
            if k == 2 and y > x and wit2 is None:
                wit2 = (x, y, b, k)
print("    smallest |b| edge with k=1 and y<x: %s ; its B1 child %s reads %s" % (
    wit1, kids(wit1[0], wit1[1])[0], "%d->%d in (3,%d,%d)" % (kids(wit1[0], wit1[1])[0][0], wit1[1], wit1[2], read(kids(wit1[0], wit1[1])[0][0], wit1[1], 3, wit1[2]))))
print("    smallest |b| edge with k=2 and y>x: %s ; its B3 child %s reads %s" % (
    wit2, kids(wit2[1], wit2[0])[2], "%d->%d in (3,%d,%d)" % (kids(wit2[1], wit2[0])[2][0], wit2[1], -wit2[2], read(kids(wit2[1], wit2[0])[2][0], wit2[1], 3, -wit2[2]))))
check(wit1 == (3, 1, -7, 1) and wit2 == (1, 3, 9, 2), "witnesses")
print("    b=-3,-5 have NO k=1 edge with y<x and b=3,5,7 have NO k=2 edge with y>x (the note's 'b<=-3', 'b>=3'")
print("    are necessary conditions only; sharpened in the note).  S4: CONFIRMED.")

# ------------------------------------------------------------------ A5
print("\n### A5  sporadic readings: exact integer linear solve (k,j<=60) + brute force max(x,y)<=20000")
# child pair as integer linear forms in x scaled by 2^k:  P = (c0 + c1*x)/2^k
generic = set()
spor = set()
for b in (1, -1):
    for k in range(1, 61):
        K = 2 ** k
        yf = (b, 3)       # 2^k y = b + 3x
        xf = (0, K)       # 2^k x
        for case in "AB":
            sf, tf = (yf, xf) if case == "A" else (xf, yf)
            fam = {"B1": ((sf[0] + 2 * tf[0], sf[1] + 2 * tf[1]), tf),
                   "B2": ((2 * sf[0] + tf[0], 2 * sf[1] + tf[1]), sf),
                   "B3": ((2 * sf[0] - tf[0], 2 * sf[1] - tf[1]), sf)}
            for name, (pf, qf) in fam.items():
                for j in range(1, 61):
                    for bp in (1, -1):
                        for rd in ("p->q", "q->p"):
                            src, tgt = (pf, qf) if rd == "p->q" else (qf, pf)
                            # 3*src + bp*2^k = 2^j * tgt   (everything times 2^k)
                            c = 3 * src[1] - 2 ** j * tgt[1]
                            d = 2 ** j * tgt[0] - 3 * src[0] - bp * K
                            if c == 0 and d == 0:
                                generic.add((b, k, case, name, rd, bp, j))
                                continue
                            if c == 0 or d % c:
                                continue
                            x = d // c
                            if x <= 0 or x % 2 == 0:
                                continue
                            v = 3 * x + b
                            if v <= 0 or v2(v) != k:
                                continue
                            y = v >> k
                            if y == x or (case == "A") != (y > x):
                                continue
                            p = (pf[0] + pf[1] * x) // K
                            q = (qf[0] + qf[1] * x) // K
                            S, T = (p, q) if rd == "p->q" else (q, p)
                            check(read(S, T, 3, bp) == j, "sporadic recheck")
                            spor.add((x, y, b, k, name, (p, q), S, T, bp, j))
print("A5 generic identities: %d (explorer: 8) ; sporadic readings for k,j<=60: %d (explorer: 5)" % (len(generic), len(spor)))
for r in sorted(spor):
    print("    %s" % (r,))
EXP = {(3, 5, 1, 1, "B3", (7, 5), 5, 7, -1, 1), (7, 5, -1, 2, "B2", (19, 7), 19, 7, -1, 3),
       (7, 5, -1, 2, "B3", (9, 7), 9, 7, 1, 2), (3, 1, -1, 3, "B1", (5, 1), 5, 1, 1, 4),
       (3, 1, -1, 3, "B3", (5, 3), 3, 5, 1, 1)}
check(spor == EXP and len(generic) == 8, "sporadic list")
# brute force, larger window
for LIM in (5000, 20000):
    E = edges31(LIM)
    nk1 = sum(1 for e in E if e[3] == 1)
    gen = 0
    sp = set()
    for (x, y, b, k) in E:
        s, t = max(x, y), min(x, y)
        for idx, pr in enumerate(kids(s, t)):
            for (S, T, bp, j) in legal31(pr):
                if k == 1 and ((idx == 1 and (S, T, bp, j) == (4 * x + b, y, b, 3)) or (idx == 2 and (S, T, bp, j) == (2 * x + b, y, -b, 2))):
                    gen += 1
                else:
                    sp.add((x, y, b, k, "B%d" % (idx + 1), pr, S, T, bp, j))
    print("    brute force max(x,y)<=%d: edges %d, k=1 edges %d, generic legal children %d, sporadic %d" % (LIM, len(E), nk1, gen, len(sp)))
    check(gen == 2 * nk1 and sp == EXP, "brute force %d" % LIM)
print("    S5: CONFIRMED (explorer: 4165 edges, 1666 k=1, 3332 generic, 5 sporadic at 5000).")

# ------------------------------------------------------------------ A6
print("\n### A6  legal sub-forest: degree census at max<=4000 and connected components")
LIM6 = 4000
pairs = {}
for (x, y, b, k) in edges31(LIM6):
    pairs.setdefault((max(x, y), min(x, y)), []).append((x, y, b, k))
deg = defaultdict(int)
parent = {}
for pr in pairs:
    for ch in kids(*pr):
        if ch in pairs:
            parent[ch] = pr
    deg[(min(r[3] for r in pairs[pr]), sum(1 for ch in kids(*pr) if ch in pairs or legal31(ch)))] += 1
print("A6 legal pairs max<=%d: %d ; (min k, legal children) census: %s" % (LIM6, len(pairs), dict(sorted(deg.items()))))
check(all((k == 1 and c == 2) or (k >= 2 and c == 0) or (k, c) == (3, 2) for (k, c) in deg) and deg[(3, 2)] == 1, "degree census")
# components (union-find over parent links inside the window)
root = {}


def find(u):
    while root.get(u, u) != u:
        u = root[u]
    return u


for ch, pr in parent.items():
    a1, b1 = find(ch), find(pr)
    if a1 != b1:
        root[a1] = b1
comp = defaultdict(list)
for pr in pairs:
    comp[find(pr)].append(pr)
sizes = defaultdict(int)
for c in comp.values():
    sizes[len(c)] += 1
big = [sorted(c) for c in comp.values() if len(c) > 3]
print("    component sizes (parent links inside the window): %s" % dict(sorted(sizes.items())))
print("    components with more than 3 pairs: %s" % big)
check(len(big) == 1 and set(big[0]) == {(3, 1), (5, 1), (5, 3), (13, 5), (7, 5), (19, 7), (9, 7)}, "root cluster unique")
print("    S6: CONFIRMED (the size-1 components are k>=4 pairs, and k=1 pairs whose children exceed the window;")
print("    size-2 components are k=1 pairs with one child in the window).")

# ------------------------------------------------------------------ A7
print("\n### A7  parent transport: cones by k at max<=4000, and ACCIDENTAL legal parents (note section 4 wording)")
cnt = defaultdict(int)
acc = []
for (x, y, b, k) in edges31(LIM6):
    s, t = max(x, y), min(x, y)
    if (s, t) == (3, 1):
        cnt[(k, "root")] += 1
        continue
    if 3 * t < s:
        cone, par = "C1", (s - 2 * t, t)
    elif 2 * t < s:
        cone, par = "C2", (t, s - 2 * t)
    else:
        cone, par = "C3", (t, 2 * t - s)
    check(cone == {1: "C3", 2: "C3", 3: "C2"}.get(k, "C1"), "cone by k")
    if k == 1:
        check(par == (x, (x - b) // 2) and read(x, (x - b) // 2, 1, -b) == 1, "k=1 parent halving")
        check(read((x - b) // 2, x, 3, (x + 3 * b) // 2) == 1, "b'' law")
    elif k == 2:
        check(par == (y, (x + b) // 2) and read((x + b) // 2, y, 3, -b) == 1, "k=2 parent")
    elif k == 3:
        check(par == (y, (x - b) // 4) and read((x - b) // 4, y, 3, b) == 1, "k=3 parent")
    else:
        check(par == (x - 2 * y, y) and 3 * (x - 2 * y) + b == (2 ** k - 6) * y, "k>=4 parent")
    cnt[(k, cone)] += 1
    if k in (1,) or k >= 4:
        L = legal31(par)
        if L:
            acc.append(((x, y, b, k), par, L))
print("A7 (k,cone) counts, max<=%d: %s" % (LIM6, dict(sorted(cnt.items()))))
print("    k=1 or k>=4 edges whose PARENT PAIR is itself a (3,+-1) edge (accidental, not by identity): %s" % acc)
check(acc == [((3, 5, 1, 1), (3, 1), [(3, 1, -1, 3)]), ((5, 1, 1, 4), (3, 1), [(3, 1, -1, 3)]),
              ((5, 7, -1, 1), (5, 3), [(3, 5, 1, 1)])], "accidental legal parents")
print("    => the note's sentence 'the parent of a Collatz k=1 edge is never a (3,+-1) edge' is WRONG as written:")
print("       3->5 (k=1) has parent (3,1) = 3->1 in (3,-1,3) (E5 reversed), 5->7 (k=1) has parent (5,3) = 3->5 in (3,1,1)")
print("       (E1 reversed), and 5->1 (k=4) has parent (3,1) (E4 reversed).  Correct statement: never by identity;")
print("       exactly these three accidental cases, which are the sporadic list of S5 read upward (so complete for")
print("       all x, not only this window).  S7: WEAKENED -> fixed in the note.")

# ------------------------------------------------------------------ A8
print("\n### A8  census s<=2000, a<=15, |b|<=15 by DIRECT EDGE ENUMERATION (no numpy pair grid)")
SM = 2000
npairs = sum(1 for s in range(1, SM + 1, 2) for t in range(1, s, 2) if math.gcd(s, t) == 1)
real = defaultdict(list)
per_a = defaultdict(int)
row3 = defaultdict(int)
tot = 0
for a in range(1, 16, 2):
    for b in range(-15, 16, 2):
        for x in range(1, SM + 1, 2):
            v = a * x + b
            if v <= 0:
                continue
            k = v2(v)
            y = v >> k
            if k == 0 or y > SM or y == x or math.gcd(x, y) != 1:
                continue
            real[(max(x, y), min(x, y))].append((x, y, a, b, k))
            tot += 1
            per_a[a] += 1
            if a == 3:
                row3[b] += 1
hist = defaultdict(int)
for L in real.values():
    hist[len(L)] += 1
print("A8 coprime odd pairs t<s<=2000: %d (explorer 405432) ; realizations %d (57837) ; pairs %d (54755) = %.2f%%"
      % (npairs, tot, len(real), 100.0 * len(real) / npairs))
print("    histogram: %s" % dict(sorted(hist.items())))
print("    per a: %s" % dict(sorted(per_a.items())))
print("    a=3 row: %s" % dict(sorted(row3.items())))
check(npairs == 405432 and tot == 57837 and len(real) == 54755 and hist[1] == 52860 and hist[2] == 1415 and max(hist) == 40, "census")
check(row3[-1] == 833 and row3[1] == 832 and per_a[3] == 10582 and per_a[1] == 12685, "a=3 row")
# overlap of b=+1 and b=-1 pairs for a=3
p1 = {pr for pr, L in real.items() if any(r[2] == 3 and r[3] == 1 for r in L)}
pm = {pr for pr, L in real.items() if any(r[2] == 3 and r[3] == -1 for r in L)}
print("    a=3: pairs with a b=+1 reading %d, with a b=-1 reading %d, union %d, intersection %s" % (len(p1), len(pm), len(p1 | pm), sorted(p1 & pm)))
check(len(p1 | pm) == 1664 and not (p1 & pm) and len(pm) == 832, "1664 pairs; (7,5) has two b=-1 readings, NOT a b=+1/b=-1 overlap")
print("    => '1664 pairs: the overlap is (7,5)' is right in count but (7,5) is a double b=-1 reading, not a b=+1/b=-1 overlap.")
print("    S8: CONFIRMED (numbers); wording sharpened.")

# ------------------------------------------------------------------ A9
print("\n### A9  hypotenuse <= 10^6 by Euclid (m,n) enumeration")
X = 10 ** 6
ntri = 0
directed = defaultdict(int)
tri_b = defaultdict(set)
m = 1
while m * m + 1 <= X:
    m += 1
    for n in range(1 + (m % 2), m, 2):
        if m * m + n * n > X:
            break
        if math.gcd(m, n) != 1:
            continue
        ntri += 1
        s, t = m + n, m - n
        for b in (-5, -3, -1, 1, 3, 5):
            for (p, q) in ((s, t), (t, s)):
                if read(p, q, 3, b) is not None:
                    directed[b] += 1
                    tri_b[b].add((s, t))
C_edge = sum(1.0 / math.sqrt(2.0 * (4.0 ** k + 9.0)) for k in range(1, 300))
print("A9 primitive triangles with C<=10^6: %d (explorer 159139)" % ntri)
print("    directed edges by b: %s ; distinct triangles b=+1: %d, b=-1: %d, b=+-1 union: %d, |b|<=5 union: %d (%.3f%%)"
      % (dict(sorted(directed.items())), len(tri_b[1]), len(tri_b[-1]), len(tri_b[1] | tri_b[-1]),
         len(set().union(*tri_b.values())), 100.0 * len(set().union(*tri_b.values())) / ntri))
print("    C_edge = %.10f ; C_edge*sqrt(X) = %.1f ; b=+1 count minus prediction = %.1f" % (C_edge, C_edge * 1000, directed[1] - C_edge * 1000))
check(ntri == 159139 and directed[1] == 507 and directed[-1] == 506 and len(tri_b[-1]) == 505 and len(tri_b[1] | tri_b[-1]) == 1012
      and len(set().union(*tri_b.values())) == 2483 and abs(C_edge - 0.5078191557) < 1e-9, "hypotenuse census")
print("    S9: CONFIRMED.")

# ------------------------------------------------------------------ A10
print("\n### A10 inverse-fibre braid as a Berggren word: closed-form exponents by direct B1 counting")
nf = 0
for b in (1, -1):
    for y in range(1, 1000, 2):
        if y % 3 == 0:
            continue
        for k0 in range(1, 12):
            if (2 ** k0 * y - b) % 3:
                continue
            x0 = (2 ** k0 * y - b) // 3
            if x0 <= 0 or x0 == y:
                continue
            if k0 == 1 and y < x0:
                continue
            for j in range(1, 4):
                kk = k0 + 2 * j
                xk = (2 ** kk * y - b) // 3
                # count B1 steps from the start pair to (xk, y), after one B2 if k0 == 1
                s, t = (y, x0) if k0 == 1 else (x0, y)
                if k0 == 1:
                    s, t = 2 * s + t, s  # B2
                steps = 0
                while (s, t) != (xk, y):
                    check(t == y and s < xk, "on the B1 ray")
                    s, t = s + 2 * t, t
                    steps += 1
                want = (4 ** j - 4) // 3 if k0 == 1 else (2 ** (k0 + 2 * j) - 2 ** k0) // 6
                check(steps == want, "exponent formula k0=%d j=%d" % (k0, j))
                nf += 1
print("A10 fibre words checked (y<1000, k0<=11, j<=3, b=+-1): %d ; exponents (4^j-4)/3 resp. (2^(k0+2j)-2^k0)/6 exact on all" % nf)
check(read(7, 11, 3, 1) == 1 and read(29, 11, 3, 1) == 3 and read(117, 11, 3, 1) == 5 and (2 * 11 + 7, 11) == (29, 11), "example")
print("    note: the k=1 element (y, x_1) of an odd-k fibre has t = x_1 != y, so it is NOT on the B1-ray t=y;")
print("    'the fibre is the subset of the B1-ray' holds for the fibre minus its k=1 element.  S10: CONFIRMED, wording sharpened.")

# ------------------------------------------------------------------ A11
print("\n### A11 consecutive orbit triangles: independent recount + EXACT linear solve of the 'proof sketch'")
LIM = 10 ** 5
np_ = 0
hits = 0
for x0 in range(1, LIM + 1, 2):
    x = x0
    prev = None
    while x != 1:
        v = 3 * x + 1
        y = v >> v2(v)
        cur = (max(x, y), min(x, y))
        if prev is not None:
            np_ += 1
            if cur in kids(*prev) or prev in kids(*cur):
                hits += 1
        prev = cur
        x = y
print("A11 consecutive pairs, odd starts <= 10^5: %d (explorer 1849203) ; parent/child incidences: %d" % (np_, hits))
check(np_ == 1849203 and hits == 0, "orbit recount")
# exact solve: x -> y (k), y -> z (j), b=+1.  Pair1 = roots(x,y), Pair2 = roots(y,z).  Scale by 2^(k+j).
sols = []
for k in range(1, 61):
    for j in range(1, 61):
        D = 2 ** (k + j)
        xf = (0, D)                       # D*x
        yf = (2 ** j, 3 * 2 ** j)         # D*y = 2^j (3x+1)
        zf = (3 + 2 ** k, 9)              # D*z = 2^k (3y+1) = 3(3x+1) + 2^k
        P1 = (yf, xf) if k == 1 else (xf, yf)
        P2 = (zf, yf) if j == 1 else (yf, zf)
        for (par, chd) in ((P1, P2), (P2, P1)):
            s, t = par
            for c in (((s[0] + 2 * t[0], s[1] + 2 * t[1]), t), ((2 * s[0] + t[0], 2 * s[1] + t[1]), s), ((2 * s[0] - t[0], 2 * s[1] - t[1]), s)):
                eqs = [(c[i][1] - chd[i][1], chd[i][0] - c[i][0]) for i in range(2)]  # coef*x = rhs, componentwise
                ok = True
                xs = None
                for (co, rh) in eqs:
                    if co == 0:
                        ok = ok and rh == 0
                        continue
                    if rh % co:
                        ok = False
                        continue
                    xv = rh // co
                    if xs is None:
                        xs = xv
                    elif xs != xv:
                        ok = False
                if not ok or xs is None or xs <= 0 or xs % 2 == 0:
                    continue
                x = xs
                v = 3 * x + 1
                if v2(v) != k:
                    continue
                y = v >> k
                w = 3 * y + 1
                if v2(w) != j:
                    continue
                sols.append((x, y, w >> j, k, j))
print("    exact solutions (x,y,z,k,j) with consecutive triangles parent/child, k,j<=60: %s" % sorted(set(sols)))
check(set(sols) == {(1, 1, 1, 2, 2)}, "S11 exact: only the fixed edge x=y=z=1 (pair (1,1) = its own B3 child), never formed by an orbit")
# ratio bound that makes the exact solve a proof for ALL k,j:  r = max/min of an edge pair, b=+1, x>=3:
#   k=1: r = y/x = (3+1/x)/2 in (3/2, 5/3];  k>=2: r = x/y = 2^k x/(3x+1) in [0.3*2^k, 2^k/3).
# A B2 child has ratio 2+1/rho in (2,3), a B3 child 2-1/rho in (1,2): child exponent j<=3 (j>=4 gives r>=4.8);
# the parent ratio rho = 1/(r2-2) resp. 1/(2-r2) then lies in (3/2,5/2] (B2, j=3) or (2,3] (B3, j=1) or
# [5/4,3/2) (B3, j=2): parent exponent <=3.  A B1 child (s+2t,t) keeps t as its smaller coordinate, so both
# pairs would have the same minimum; with k>=2 that forces z = x+2y, ratio > 2, but j=1 has ratio <= 5/3; with
# k=1 it forces z = x < y, contradicting y -> z with the smaller coordinate x retained.  Hence every candidate
# has k,j <= 3, inside the exact solve.  Numeric check of the ratio intervals:
lo = defaultdict(lambda: 10.0 ** 9)
hi = defaultdict(float)
for x in range(3, 100001, 2):
    v = 3 * x + 1
    k = v2(v)
    y = v >> k
    r = max(x, y) / min(x, y)
    lo[k] = min(lo[k], r)
    hi[k] = max(hi[k], r)
bad = 0
for k in sorted(lo):
    if k == 1:
        ok = 1.5 < lo[k] and hi[k] <= 5 / 3 + 1e-12
    else:
        ok = 0.3 * 2 ** k <= lo[k] and hi[k] < 2 ** k / 3
    bad += 0 if ok else 1
print("    ratio intervals per k (x<=10^5, b=+1): %s" % {k: (round(lo[k], 4), round(hi[k], 4)) for k in sorted(lo) if k <= 8})
print("    violations of k=1: (3/2,5/3], k>=2: [0.3*2^k, 2^k/3): %d" % bad)
check(bad == 0, "ratio intervals")
print("    => B2/B3 parent-child needs (k,j) in {(1,3),(3,3),(3,1),(2,2)} up to direction, B1 never; all inside the")
print("       exact solve, whose only solution is x=y=z=1.  S11: CONFIRMED and upgraded from 'PROVED sketch' to PROVED")
print("       (ratio bound + exact linear solve), with the 10^5 census as the FINITE-EXACT companion.")

print("\nALL AUDIT CHECKS PASSED")
