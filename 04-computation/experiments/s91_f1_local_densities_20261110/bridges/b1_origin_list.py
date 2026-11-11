#!/usr/bin/env python3
"""S91 bridges (1): origin of the owner's list 2, 3, 5, 7, 13, 19, 27, 1039.

Source claimed by the digests: J. Allikvere, "Minimal Trees for the Euler Sombor Index Revisited", MATCH 97 (2026)
1281-1293.  Euler Sombor index EP(G) = sum over edges uv of sqrt(d(u)^2 + d(u) d(v) + d(v)^2).

What this script checks, with its own code (exact integer / interval arithmetic; no floating point decides anything):
  A. the Loeschian norms a^2 + ab + b^2 at (1,1), (1,2), (1,3), (2,3), (3,3) are 3, 7, 13, 19, 27 (and (2,2) gives 12);
  B. eps = 2 sqrt 19 - 5 sqrt 3 = sqrt 76 - sqrt 75 = 1/(sqrt 76 + sqrt 75), delta = sqrt 13 + 2 sqrt 3 - sqrt 7 - sqrt 19,
     and 0 < eps < delta, by certified interval enclosures (integer square roots at 60 digits);
  C. the integers printed in the abstract (read from the pdftotext extraction if that file is present; the result of
     that reading is also hard-coded below so that the script is self-contained);
  D. the number of unicyclic graphs on n = 3..10 vertices, three independent ways:
       (i)   generation by pendant-vertex extension with an exact canonical form (cycle of AHU codes modulo the
             dihedral group),
       (ii)  Polya / Burnside count of bracelets of rooted trees,
       (iii) the mass formula  sum over classes of n!/|Aut| = number of labelled connected unicyclic graphs,
       (iv)  for n <= 7, brute force over all edge sets with n edges; and the networkx graph atlas for n <= 7;
     expected 1, 2, 5, 13, 33, 89, 240, 657 (n = 3..10) and 2 + 5 + 13 + 33 + 89 + 240 + 657 = 1039;
  E. the statements of the abstract about trees, on all trees with 7 <= n <= 14 (networkx.nonisomorphic_trees, exact
     vectors over square roots of squarefree integers): ranks 1-5 of Theorem 1, the unattained value
     4 sqrt 7 + 4 sqrt 19 - 18 sqrt 3 + 2 sqrt 3 n, and ranks 5-13 = base + r delta + c eps for 11 <= n <= 14.

Exit status 0 and a final line ALL CHECKS PASSED iff every check holds.
"""
import itertools
import os
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from math import comb, factorial, gcd, isqrt

import networkx as nx

getcontext().prec = 60
RESULTS = []


def check(name, cond, detail=""):
    RESULTS.append(bool(cond))
    print(("PASS  " if cond else "FAIL  ") + name + (("   [" + detail + "]") if detail else ""))


# ---------------------------------------------------------------- A. norms
def norm(a, b):
    return a * a + a * b + b * b


pairs = [(1, 1), (1, 2), (1, 3), (2,3), (3, 3)]
check("A1 norms at (1,1),(1,2),(1,3),(2,3),(3,3) are 3,7,13,19,27", [norm(*p) for p in pairs] == [3, 7, 13, 19, 27],
      str([norm(*p) for p in pairs]))
alldeg3 = sorted({norm(a, b) for a in range(1, 4) for b in range(a, 4)})
check("A2 radicands for degrees <= 3 are {3,7,12,13,19,27}; the list omits only 12 = N(2,2)", alldeg3 == [3, 7, 12, 13, 19, 27],
      str(alldeg3))
check("A3 (a^3 - b^3) = (a - b) N(a, b) for 1 <= b < a <= 30",
      all(a ** 3 - b ** 3 == (a - b) * norm(a, b) for a in range(1, 31) for b in range(1, a)))
check("A4 7 | 2^3 - 1, 13 | 3^3 - 1, 19 = 3^3 - 2^3, 1039 = N(15, 22) = 26^2 + 3*11^2, 1039 prime",
      (2 ** 3 - 1) % 7 == 0 and (3 ** 3 - 1) % 13 == 0 and 27 - 8 == 19 and norm(15, 22) == 1039
      and 26 ** 2 + 3 * 11 ** 2 == 1039 and all(1039 % d for d in range(2, 33)))
reps1039 = sorted((a, b) for a in range(0, 40) for b in range(a, 40) if norm(a, b) == 1039)
check("A5 1039 has exactly one representation a^2 + ab + b^2 with 0 <= a <= b", reps1039 == [(15, 22)], str(reps1039))
x27 = {p: [(x, y) for x in range(0, 40) for y in range(0, 8) if x * x + 27 * y * y == p] for p in (7, 13, 19, 1039)}
check("A6 none of 7, 13, 19, 1039 is x^2 + 27 y^2 (2 is not a cubic residue)", all(not v for v in x27.values()))
check("A7 cube test: 2^((p-1)/3) != 1 mod p for p = 7, 13, 19, 1039",
      all(pow(2, (p - 1) // 3, p) != 1 for p in (7, 13, 19, 1039)))

# ---------------------------------------------------------------- B. eps and delta, certified
SCALE = 10 ** 60


def sqrt_interval(n):
    """integers (lo, hi) with lo <= sqrt(n) * 10^60 <= hi."""
    lo = isqrt(n * SCALE * SCALE)
    return lo, lo + 1


def lin_interval(terms):
    """terms: list of (coefficient, radicand); returns certified (lo, hi) of sum c sqrt(n), scaled by 10^60."""
    lo = hi = 0
    for c, n in terms:
        a, b = sqrt_interval(n)
        if c >= 0:
            lo += c * a
            hi += c * b
        else:
            lo += c * b
            hi += c * a
    return lo, hi


eps = lin_interval([(2, 19), (-5, 3)])
eps2 = lin_interval([(1, 76), (-1, 75)])
delta = lin_interval([(1, 13), (2, 3), (-1, 7), (-1, 19)])
dme = lin_interval([(1, 13), (7, 3), (-1, 7), (-3, 19)])
check("B1 4*19 = 76 = 75 + 1 = 3*5^2 + 1, so eps = 2 sqrt 19 - 5 sqrt 3 = sqrt 76 - sqrt 75 = 1/(sqrt 76 + sqrt 75)",
      4 * 19 == 76 and 25 * 3 == 75 and eps[0] <= eps2[1] and eps2[0] <= eps[1])
check("B2 eps > 0 (certified)", eps[0] > 0, "eps in [%.12f, %.12f]" % (eps[0] / SCALE, eps[1] / SCALE))
check("B3 delta - eps = sqrt 13 + 7 sqrt 3 - sqrt 7 - 3 sqrt 19 > 0 (certified)", dme[0] > 0,
      "delta = %.12f, delta - eps = %.12f" % (delta[0] / SCALE, dme[0] / SCALE))
# the same sign by pure integer squaring:  sqrt13 + 7 sqrt3 > sqrt7 + 3 sqrt19
#   LHS^2 = 160 + 14 sqrt 39,  RHS^2 = 178 + 6 sqrt 133;  need 14 sqrt39 - 6 sqrt133 > 18, i.e. 7 sqrt39 > 9 + 3 sqrt133
#   squares: 1911 > 81 + 1197 + 54 sqrt133  <=>  633 > 54 sqrt 133  <=>  633^2 = 400689 > 54^2 * 133 = 387828
check("B4 the same inequality by integer squaring (400689 > 387828)", 633 ** 2 == 400689 and 54 ** 2 * 133 == 387828 and 400689 > 387828
      and 13 + 147 == 160 and 7 + 9 * 19 == 178 and 49 * 39 == 1911 and 9 * 133 == 1197)
# reduced forms of discriminant -75
red = [(a, b, c) for a in range(1, 6) for b in range(-a + 1, a + 1) for c in range(a, 40)
       if b * b - 4 * a * c == -75 and (abs(b) <= a <= c) and not (b < 0 and (abs(b) == a or a == c))]
prim = sorted(f for f in red if gcd(gcd(f[0], f[1]), f[2]) == 1)
check("B5 reduced forms of discriminant -75: primitive (1,1,19), (3,3,7); and the imprimitive 5*(1,1,1) "
      "(the digest's list omits the qualifier 'primitive')", prim == [(1, 1, 19), (3, 3, 7)]
      and sorted(red) == [(1, 1, 19), (3, 3, 7), (5, 5, 5)], str(sorted(red)))

# ---------------------------------------------------------------- C. integers of the abstract
OWNER = [2, 3, 5, 7, 13, 19, 27, 1039]
# Hard-coded reading of the abstract (numerals only; the exponent 2 of d(u)^2 included; the citation
# "94 (2025) 525-548" excluded): formulas  4 sqrt7 + 4 sqrt19 - 18 sqrt3 + 2 sqrt3 n,  2 sqrt3, (2,2), sqrt27, S(3,3),
# delta = sqrt13 + 2 sqrt3 - sqrt7 - sqrt19, eps = 2 sqrt19 - 5 sqrt3, 0 < eps < delta; "Theorem 3", "Theorem 4";
# ranges n >= 7, n >= 11, n = 4, 7 <= n <= 18, 4 <= n <= 10; census sizes 204,990 and 1,039.
ABSTRACT_NUMERALS = {0, 2, 3, 4, 5, 7, 10, 11, 13, 18, 19, 27, 1039, 204990}
check("C1 every integer of the owner's list is printed in the abstract", set(OWNER) <= ABSTRACT_NUMERALS)
extra = sorted(ABSTRACT_NUMERALS - set(OWNER))
print("      abstract numerals NOT in the owner's list:", extra,
      " (so 'exactly the integers of the abstract' is an overstatement)")
in_delta_eps = {13, 2, 3, 7, 19} | {2, 19, 5, 3} | {27}
check("C2 the list = integers of delta, eps and sqrt 27, together with the census size 1039",
      sorted(in_delta_eps | {1039}) == OWNER)
txt_path = os.path.join(os.environ.get("S91_SCRATCH", ""), "es_revisited_out", "layout_utf8.txt")
if os.path.isfile(txt_path):
    raw = open(txt_path, encoding="utf-8").read()
    body = raw[raw.index("Abstract"):raw.index("1 Introduction")]
    body = body[:body.index("This work is licensed")] + body[body.index("unicyclic graphs, which"):]
    body = body.replace("d(u)2", "d(u) 2").replace("d(v)2", "d(v) 2")
    nums = {int(t.replace(",", "")) for t in re.findall(r"\d[\d,]*", body)}
    cite = {94, 2025, 525, 548}
    found = nums - cite
    check("C3 regex reading of the extracted abstract agrees with the hard-coded set", found == ABSTRACT_NUMERALS,
          "found %s" % sorted(found))
else:
    print("      (C3 skipped: extraction file not present; set S91_SCRATCH to the scratchpad folder to run it)")

# ---------------------------------------------------------------- D. unicyclic graphs
NMAX = 10


def tree_code(adj, v, parent, banned):
    kids = sorted(tree_code(adj, c, v, banned) for c in adj[v] if c != parent and c not in banned)
    return "(" + "".join(kids) + ")"


def tree_aut(adj, v, parent, banned):
    """order of the automorphism group of the rooted tree hanging at v (children avoid 'banned' and the parent)."""
    kids = [c for c in adj[v] if c != parent and c not in banned]
    codes = {}
    a = 1
    for c in kids:
        a *= tree_aut(adj, c, v, banned)
        k = tree_code(adj, c, v, banned)
        codes[k] = codes.get(k, 0) + 1
    for m in codes.values():
        a *= factorial(m)
    return a


def cycle_of(adj):
    n = len(adj)
    deg = [len(a) for a in adj]
    alive = [True] * n
    stack = [v for v in range(n) if deg[v] == 1]
    while stack:
        v = stack.pop()
        alive[v] = False
        for w in adj[v]:
            if alive[w]:
                deg[w] -= 1
                if deg[w] == 1:
                    stack.append(w)
    cyc = [v for v in range(n) if alive[v]]
    start = cyc[0]
    order = [start]
    prev = None
    cur = start
    while True:
        nxt = [w for w in adj[cur] if alive[w] and w != prev]
        if len(order) == len(cyc):
            break
        w = nxt[0]
        order.append(w)
        prev, cur = cur, w
    return order


def canon_unicyclic(adj):
    """(canonical form, |Aut|) of a connected unicyclic graph given by adjacency sets."""
    cyc = cycle_of(adj)
    k = len(cyc)
    cs = set(cyc)
    codes = [tree_code(adj, v, None, cs - {v}) for v in cyc]
    auts = 1
    for v in cyc:
        auts *= tree_aut(adj, v, None, cs - {v})
    variants = []
    for seq in (codes, codes[::-1]):
        for r in range(k):
            variants.append(tuple(seq[r:] + seq[:r]))
    best = min(variants)
    sym = sum(1 for t in variants if t == best)          # stabiliser of the coloured cycle in the dihedral group D_k
    return best, auts * sym


def is_connected(adj):
    n = len(adj)
    seen = {0}
    st = [0]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                st.append(w)
    return len(seen) == n


# (i) generation
classes = {3: {}}
c3 = [{1, 2}, {0, 2}, {0, 1}]
classes[3][canon_unicyclic(c3)[0]] = (c3, canon_unicyclic(c3)[1])
for n in range(4, NMAX + 1):
    cur = {}
    cyc = [{(i - 1) % n, (i + 1) % n} for i in range(n)]
    cf, au = canon_unicyclic(cyc)
    cur[cf] = (cyc, au)
    for adj, _ in classes[n - 1].values():
        for v in range(n - 1):
            new = [set(a) for a in adj] + [{v}]
            new[v].add(n - 1)
            cf, au = canon_unicyclic(new)
            if cf not in cur:
                cur[cf] = (new, au)
    classes[n] = cur
gen_counts = [len(classes[n]) for n in range(3, NMAX + 1)]
EXPECT = [1, 2, 5, 13, 33, 89, 240, 657]
check("D1 generation with canonical forms: unicyclic counts for n = 3..10", gen_counts == EXPECT, str(gen_counts))

# (ii) Polya count: bracelets of rooted trees
r = [0, 1]
for n in range(1, NMAX):
    s = 0
    for k in range(1, n + 1):
        s += sum(d * r[d] for d in range(1, k + 1) if k % d == 0) * r[n - k + 1]
    r.append(s // n)
check("D2 rooted trees r_1..r_10 = 1,1,2,4,9,20,48,115,286,719", r[1:] == [1, 1, 2, 4, 9, 20, 48, 115, 286, 719], str(r[1:]))


def pmul(a, b):
    out = [Fraction(0)] * (NMAX + 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y and i + j <= NMAX:
                    out[i + j] += x * y
    return out


def ppow(a, e):
    out = [Fraction(0)] * (NMAX + 1)
    out[0] = Fraction(1)
    for _ in range(e):
        out = pmul(out, a)
    return out


def subst(d):
    out = [Fraction(0)] * (NMAX + 1)
    for j in range(1, NMAX + 1):
        if j * d <= NMAX:
            out[j * d] = Fraction(r[j])
    return out


def phi(n):
    return sum(1 for i in range(1, n + 1) if gcd(i, n) == 1)


polya = [Fraction(0)] * (NMAX + 1)
for k in range(3, NMAX + 1):
    tot = [Fraction(0)] * (NMAX + 1)
    for d in range(1, k + 1):
        if k % d == 0:
            term = ppow(subst(d), k // d)
            tot = [t + phi(d) * u / (2 * k) for t, u in zip(tot, term)]
    if k % 2:
        term = pmul(subst(1), ppow(subst(2), (k - 1) // 2))
        tot = [t + u / 2 for t, u in zip(tot, term)]
    else:
        t1 = pmul(ppow(subst(1), 2), ppow(subst(2), (k - 2) // 2))
        t2 = ppow(subst(2), k // 2)
        tot = [t + (u + v) / 4 for t, u, v in zip(tot, t1, t2)]
    polya = [a + b for a, b in zip(polya, tot)]
polya_counts = [polya[n] for n in range(3, NMAX + 1)]
check("D3 Polya (dihedral cycle index of rooted trees): the same counts", polya_counts == EXPECT, str([int(x) for x in polya_counts]))

# (iii) mass formula
mass_ok = True
labelled = []
for n in range(3, NMAX + 1):
    lhs = sum(Fraction(factorial(n), au) for _, au in classes[n].values())
    rhs = sum(Fraction(comb(n, k) * factorial(k - 1) * k * n ** (n - k), 2 * n) for k in range(3, n + 1))
    labelled.append(int(rhs))
    mass_ok &= (lhs == rhs and rhs.denominator == 1)
check("D4 mass formula sum n!/|Aut| = labelled connected unicyclic graphs (n = 3..10)", mass_ok, str(labelled))
check("D5 labelled counts agree with (n-1)!/2 * sum_{k<=n-3} n^k/k!",
      all(Fraction(factorial(n - 1), 2) * sum(Fraction(n ** k, factorial(k)) for k in range(0, n - 2)) == labelled[n - 3]
          for n in range(3, NMAX + 1)))

# (iv) brute force for n <= 7 and the graph atlas
bf_counts, bf_labelled = [], []
for n in range(3, 8):
    edges = list(itertools.combinations(range(n), 2))
    seen = set()
    lab = 0
    for sub in itertools.combinations(edges, n):
        adj = [set() for _ in range(n)]
        for a, b in sub:
            adj[a].add(b)
            adj[b].add(a)
        if min(len(a) for a in adj) == 0 or not is_connected(adj):
            continue
        lab += 1
        seen.add(canon_unicyclic(adj)[0])
    bf_counts.append(len(seen))
    bf_labelled.append(lab)
check("D6 brute force over all n-edge graphs, n = 3..7: classes 1,2,5,13,33 and labelled counts", bf_counts == EXPECT[:5]
      and bf_labelled == labelled[:5], "%s %s" % (bf_counts, bf_labelled))
atlas = [0] * 8
for g in nx.graph_atlas_g():
    n = g.number_of_nodes()
    if 3 <= n <= 7 and g.number_of_edges() == n and nx.is_connected(g):
        atlas[n] += 1
check("D7 networkx graph atlas: connected graphs with n vertices and n edges, n = 3..7", atlas[3:8] == EXPECT[:5], str(atlas[3:8]))
# VF2 cross-check of the canonical form at n = 8: the 89 class representatives are pairwise non-isomorphic
reps8 = [nx.Graph([(a, b) for a in range(8) for b in adj[a] if a < b]) for adj, _ in classes[8].values()]
noniso = all(not nx.is_isomorphic(reps8[i], reps8[j]) for i in range(len(reps8)) for j in range(i))
check("D8 the 89 representatives at n = 8 are pairwise non-isomorphic (networkx VF2)", noniso)
check("D9 2 + 5 + 13 + 33 + 89 + 240 + 657 = 1039 (unicyclic graphs with 4 <= n <= 10)", sum(gen_counts[1:]) == 1039,
      "sum = %d; with n = 3: %d" % (sum(gen_counts[1:]), sum(gen_counts)))


# ---------------------------------------------------------------- E. trees
def squarefree_split(n):
    m, q = 1, n
    d = 2
    while d * d <= q:
        while q % (d * d) == 0:
            q //= d * d
            m *= d
        d += 1
    return m, q


SQ = {}


def dsqrt(q):
    if q not in SQ:
        SQ[q] = Decimal(q).sqrt()
    return SQ[q]


def ep_vector(g):
    vec = {}
    for u, v in g.edges():
        a, b = g.degree(u), g.degree(v)
        m, q = squarefree_split(a * a + a * b + b * b)
        vec[q] = vec.get(q, 0) + m
    return tuple(sorted(vec.items()))


def val(vec):
    return sum(Decimal(c) * dsqrt(q) for q, c in vec)


def vec_of(coeffs, n):
    """coeffs: dict radicand -> coefficient, plus the term 2 sqrt 3 * n; returns the canonical vector."""
    d = dict(coeffs)
    d[3] = d.get(3, 0) + 2 * n
    return tuple(sorted((q, c) for q, c in d.items() if c))


tree_counts = []
ok_rank, ok_unatt, ok_513 = True, True, True
detail = []
for n in range(7, 15):
    vecs = {}
    cnt = 0
    for t in nx.nonisomorphic_trees(n):
        cnt += 1
        v = ep_vector(t)
        vecs[v] = vecs.get(v, 0) + 1
    tree_counts.append(cnt)
    order = sorted(vecs, key=val)
    # distinct vectors have distinct values (square roots of distinct squarefree integers are Q-independent);
    # numerically the consecutive gaps must be far above the working precision
    gaps = [val(order[i + 1]) - val(order[i]) for i in range(len(order) - 1)]
    ok_rank &= min(gaps) > Decimal(10) ** -30
    th = [vec_of({7: 2, 3: -6}, n), vec_of({19: 3, 7: 3, 3: -14}, n), vec_of({13: 1, 19: 2, 7: 2, 3: -12}, n),
          vec_of({13: 2, 19: 1, 7: 1, 3: -10}, n)]
    ok_rank &= (order[:4] == th)
    printed = vec_of({7: 4, 19: 4, 3: -18}, n)
    ok_unatt &= (printed not in vecs)
    if n >= 10:
        ok_rank &= (order[4] == vec_of({7: 4, 19: 4, 3: -17}, n))
    if n >= 11:
        base = {7: 4, 19: 4, 3: -17}
        fam = []
        for rr, cc in [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1), (4, 1)]:
            d = dict(base)
            # r delta = r (sqrt13 + 2 sqrt3 - sqrt7 - sqrt19),  c eps = c (2 sqrt19 - 5 sqrt3)
            d[13] = d.get(13, 0) + rr
            d[3] += 2 * rr - 5 * cc
            d[7] -= rr
            d[19] += -rr + 2 * cc
            fam.append(vec_of(d, n))
        ok_513 &= (order[4:13] == fam)
    detail.append("n=%d: %d trees, %d values, min gap %.3e" % (n, cnt, len(order), min(gaps)))
for line in detail:
    print("      " + line)
check("E1 numbers of trees for n = 7..14 are 11, 23, 47, 106, 235, 551, 1301, 3159", tree_counts == [11, 23, 47, 106, 235, 551, 1301, 3159],
      str(tree_counts))
check("E2 ranks 1-4 (n = 7..14) and rank 5 (n = 10..14) are the values of Theorem 1", ok_rank)
check("E3 the value 4 sqrt7 + 4 sqrt19 - 18 sqrt3 + 2 sqrt3 n is attained by no tree with 7 <= n <= 14", ok_unatt)
check("E4 ranks 5-13 for n = 11..14 are base + r delta + c eps in the order (0,0),(0,1),(1,0),(1,1),(2,0),(2,1),(3,0),(3,1),(4,1)", ok_513)

if all(RESULTS):
    print("ALL CHECKS PASSED (%d checks)" % len(RESULTS))
    sys.exit(0)
print("CHECK FAILED (%d of %d)" % (RESULTS.count(False), len(RESULTS)))
sys.exit(1)
