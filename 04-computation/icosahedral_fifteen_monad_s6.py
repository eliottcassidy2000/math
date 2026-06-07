#!/usr/bin/env python3
"""
icosahedral_fifteen_monad_s6.py  -- monad-explorer-2026-06-07-S6

Pushes THM-436 sec.5 (the UNUSED Klein/icosahedral handle) into verified
combinatorics. Central object: the "15" of the 5-point cause.

Claims verified here:
 (A) #{overlapping cyclic-triangle pairs on [n]} = 15 * C(n,5)
     (unordered pairs of 3-subsets {X,Y} with |X cap Y| = 1).
     Per-5-set count = 15.  Matches THM-436's 0,0,15,90,315,840 (n=3..8).
 (B) The ORIENTED count = 60 * C(n,5); per-5-set = 60 = |A_5|.
 (C) CANONICAL BIJECTION (no choices): overlapping pair on a 5-set
        {X,Y}  <->  involution (ab)(cd) of A_5
     via (shared vertex v = X cap Y) <-> (fixed point f),
         ({X\v, Y\v} the two off-pairs) <-> ({a,b},{c,d} the two transpositions).
     Both sides have 15 elements; the map is a bijection.
 (D) A_5 conjugacy-class / icosahedral-axis dictionary:
        15 involutions   = 15 two-fold axes  (30 edge-points)   -> deg-15 invariant
        20 three-cycles  = 10 three-fold axes (20 face-points)  -> deg-10 invariant
        24 five-cycles   =  6 five-fold axes  (12 vertices)      -> deg-6  invariant
     i.e. the three nontrivial icosahedral invariant degrees {6,10,15}
     are exactly the AXIS COUNTS, and {6*4,10*2,15*1}={24,20,15} the class sizes.
 (E) Commutator map on ORIENTED overlapping pairs: tabulate [sigma_X, sigma_Y].
 (F) Round tournament C_5: among its 5 cyclic triangles, how many pairs overlap
     in exactly one vertex (the tournament realisation of the 5-point cause).
"""
from itertools import combinations, permutations
from math import comb

# ---------- permutation helpers (on a small ground set, as dicts) ----------
def cyc_to_perm(cycle, ground):
    """cycle = tuple (a,b,c...) meaning a->b->c->...->a; ground = full point set."""
    p = {x: x for x in ground}
    L = len(cycle)
    for i in range(L):
        p[cycle[i]] = cycle[(i + 1) % L]
    return p

def compose(p, q):
    """(p after q): x -> p[q[x]]."""
    return {x: p[q[x]] for x in q}

def inverse(p):
    return {v: k for k, v in p.items()}

def commutator(a, b):
    # [a,b] = a b a^-1 b^-1
    return compose(compose(a, b), compose(inverse(a), inverse(b)))

def cycle_type(p):
    seen = set()
    ct = []
    for x in p:
        if x in seen:
            continue
        L = 0
        y = x
        while y not in seen:
            seen.add(y)
            y = p[y]
            L += 1
        if L > 1:
            ct.append(L)
    return tuple(sorted(ct, reverse=True))

# ---------- (A)/(B) overlap counts ----------
def count_overlap_pairs(n, oriented=False):
    pts = list(range(n))
    if not oriented:
        triples = list(combinations(pts, 3))
        cnt = 0
        for i in range(len(triples)):
            for j in range(i + 1, len(triples)):
                if len(set(triples[i]) & set(triples[j])) == 1:
                    cnt += 1
        return cnt
    else:
        # oriented = pick a cyclic orientation of each 3-subset (2 each)
        # represent an oriented triangle by frozenset of its directed arcs
        def oriented_tris(s):
            a, b, c = s
            return [frozenset([(a, b), (b, c), (c, a)]),
                    frozenset([(a, c), (c, b), (b, a)])]
        triples = list(combinations(pts, 3))
        otris = []
        for t in triples:
            for o in oriented_tris(t):
                otris.append((set(t), o))
        cnt = 0
        for i in range(len(otris)):
            for j in range(i + 1, len(otris)):
                if len(otris[i][0] & otris[j][0]) == 1:
                    cnt += 1
        return cnt

print("=" * 70)
print("(A)/(B)  overlap-pair counts vs 15*C(n,5) and 60*C(n,5)")
print("=" * 70)
print(f"{'n':>2} {'unoriented':>11} {'15*C(n,5)':>10} {'oriented':>9} {'60*C(n,5)':>10}")
for n in range(3, 9):
    u = count_overlap_pairs(n, oriented=False)
    o = count_overlap_pairs(n, oriented=True) if n <= 7 else None
    pred_u = 15 * comb(n, 5)
    pred_o = 60 * comb(n, 5)
    os = f"{o:>9}" if o is not None else f"{'(skip)':>9}"
    print(f"{n:>2} {u:>11} {pred_u:>10} {os} {pred_o:>10}  "
          f"{'OK' if u == pred_u and (o is None or o == pred_o) else 'MISMATCH'}")

# ---------- (C) canonical bijection on a 5-set ----------
print()
print("=" * 70)
print("(C)  canonical bijection: overlapping pair <-> involution of A_5")
print("=" * 70)
S = [0, 1, 2, 3, 4]
# enumerate overlapping unordered pairs {X,Y} on S
triples = list(combinations(S, 3))
overlap_pairs = []
for i in range(len(triples)):
    for j in range(i + 1, len(triples)):
        X, Y = set(triples[i]), set(triples[j])
        if len(X & Y) == 1:
            overlap_pairs.append((triples[i], triples[j]))

# enumerate involutions of A_5 = (ab)(cd)
involutions = []
for f in S:                                   # fixed point
    rest = [x for x in S if x != f]
    # perfect matchings of 4 elements into 2 pairs: 3 of them
    a = rest[0]
    for b in rest[1:]:
        cd = [x for x in rest if x not in (a, b)]
        c, d = cd
        inv = cyc_to_perm((a, b), S)
        inv = compose(inv, cyc_to_perm((c, d), S))
        involutions.append((f, frozenset([frozenset((a, b)), frozenset((c, d))]), inv))

print(f"#overlapping pairs on 5-set = {len(overlap_pairs)}")
print(f"#involutions of A_5         = {len(involutions)}")

# build the map and check bijectivity
def pair_signature(X, Y):
    v = (set(X) & set(Y)).pop()
    px = frozenset(set(X) - {v})
    py = frozenset(set(Y) - {v})
    return (v, frozenset([px, py]))

pair_sigs = set(pair_signature(X, Y) for (X, Y) in overlap_pairs)
inv_sigs = set((f, m) for (f, m, _) in involutions)
print(f"distinct pair signatures   = {len(pair_sigs)}")
print(f"distinct inv signatures    = {len(inv_sigs)}")
print(f"BIJECTION (signatures equal as sets): {pair_sigs == inv_sigs}")

# ---------- (D) class / axis dictionary ----------
print()
print("=" * 70)
print("(D)  A_5 conjugacy classes <-> icosahedral axes <-> invariant degrees")
print("=" * 70)
# build all of A_5 and bucket by cycle type
A5 = []
for perm in permutations(S):
    p = {S[i]: perm[i] for i in range(5)}
    # even permutations only
    # parity via cycle type: even iff (5 - #cycles) is even
    seen = set(); ncyc = 0
    for x in S:
        if x in seen: continue
        ncyc += 1; y = x
        while y not in seen:
            seen.add(y); y = p[y]
    if (5 - ncyc) % 2 == 0:
        A5.append(p)
from collections import Counter
ctc = Counter(cycle_type(p) for p in A5)
print(f"|A_5| = {len(A5)}")
print("cycle-type class sizes:", dict(ctc))
print()
print(" class            size  axis-type      #axes  rots/axis  #points  inv-deg")
rows = [
    ("identity ()",        ctc[()],    "-",          0, 0, 0,  0),
    ("involution (ab)(cd)", ctc[(2,2)], "2-fold",    15, 1, 30, 15),
    ("3-cycle (abc)",       ctc[(3,)],  "3-fold",    10, 2, 20, 10),
    ("5-cycle (abcde)",     ctc[(5,)],  "5-fold",     6, 4, 12,  6),
]
for name, size, axt, naxes, rpa, npts, deg in rows:
    chk = "" if name.startswith("identity") else \
          ("OK" if size == naxes * rpa else "??")
    print(f" {name:18} {size:4}  {axt:10} {naxes:6} {rpa:9} {npts:7} "
          f"{deg if deg else '-':>6}  {chk}")
print("note: class size = (#axes)*(rots/axis); inv-degree = #axes;")
print("      Euler: V-E+F = 12-30+20 = 2; deg-set {2,6,10,15} (2 = the sphere).")

# ---------- (E) commutator map on oriented overlapping pairs ----------
print()
print("=" * 70)
print("(E)  commutators [sigma_X, sigma_Y] over oriented overlapping pairs (5-set)")
print("=" * 70)
def oriented_tris(s):
    a, b, c = s
    return [(a, b, c), (a, c, b)]
comm_types = Counter()
total = 0
nonid = 0
generate_full = set()
for (X, Y) in overlap_pairs:
    for ox in oriented_tris(X):
        for oy in oriented_tris(Y):
            sx = cyc_to_perm(ox, S)
            sy = cyc_to_perm(oy, S)
            c = commutator(sx, sy)
            comm_types[cycle_type(c)] += 1
            total += 1
            if cycle_type(c) != ():
                nonid += 1
print(f"total oriented overlapping pairs = {total} (= 60? {total==60})")
print("commutator cycle-type distribution:", dict(comm_types))
print(f"non-identity commutators: {nonid}/{total}")

# subgroup generated by one canonical overlapping oriented pair
def subgroup(gens):
    grp = {tuple(sorted({x:x for x in S}.items()))}
    # represent perms as frozenset of items
    elems = {frozenset({x:x for x in S}.items())}
    frontier = list(elems)
    gp = [g for g in gens]
    changed = True
    allp = set(frozenset(g.items()) for g in gens)
    allp.add(frozenset({x:x for x in S}.items()))
    while changed:
        changed = False
        cur = list(allp)
        for e in cur:
            ed = dict(e)
            for g in gens:
                ne = frozenset(compose(g, ed).items())
                if ne not in allp:
                    allp.add(ne); changed = True
    return allp
g1 = cyc_to_perm((0,1,2), S)
g2 = cyc_to_perm((2,3,4), S)
sub = subgroup([g1, g2])
print(f"<(012),(234)> order = {len(sub)}  (A_5 has order 60: "
      f"{'A_5' if len(sub)==60 else 'NOT A_5'})")

# ---------- (F) round tournament C_5 ----------
print()
print("=" * 70)
print("(F)  round tournament C_5: overlapping cyclic-triangle pairs")
print("=" * 70)
# C_5: arc i->j iff (j-i) mod 5 in {1,2}
n = 5
def arc(i, j):
    return (j - i) % n in (1, 2)
cyc_tris = []
for t in combinations(range(n), 3):
    # a 3-subset is a cyclic triangle iff it is a directed 3-cycle
    a, b, c = t
    # check the 3 possible cyclic orders for a directed 3-cycle
    is_cyclic = False
    for perm in [(a, b, c), (a, c, b)]:
        x, y, z = perm
        if arc(x, y) and arc(y, z) and arc(z, x):
            is_cyclic = True
    if is_cyclic:
        cyc_tris.append(t)
print(f"#cyclic triangles in C_5 = {len(cyc_tris)} (expect 5): {cyc_tris}")
ov = 0; share0 = 0; share2 = 0
for i in range(len(cyc_tris)):
    for j in range(i + 1, len(cyc_tris)):
        s = len(set(cyc_tris[i]) & set(cyc_tris[j]))
        if s == 1: ov += 1
        elif s == 0: share0 += 1
        elif s == 2: share2 += 1
print(f"pairs of cyclic triangles sharing exactly 1 vertex = {ov}")
print(f"   sharing 0 vertices = {share0}, sharing 2 vertices = {share2}, "
      f"total pairs C(5,2)={comb(5,2)}")
print(f"=> the C_5 tournament realises {ov} overlapping cyclic-triangle pairs"
      f" (the on-tournament shadow of the 15 abstract ones).")

# ---------- (G) Sylow / axis-count refinement ----------
print()
print("=" * 70)
print("(G)  icosahedral axis-counts {6,10,15} as Sylow / involution counts of A_5")
print("=" * 70)
def subgroup_from(gens):
    allp = {tuple(sorted({x:x for x in S}.items())): {x:x for x in S}}
    changed = True
    while changed:
        changed = False
        for e in list(allp.values()):
            for g in gens:
                ne = compose(g, e); k = tuple(sorted(ne.items()))
                if k not in allp:
                    allp[k] = ne; changed = True
    return frozenset(allp.keys())
five  = [p for p in A5 if cycle_type(p) == (5,)]
three = [p for p in A5 if cycle_type(p) == (3,)]
invs  = [p for p in A5 if cycle_type(p) == (2, 2)]
syl5 = set(subgroup_from([g]) for g in five)
syl3 = set(subgroup_from([g]) for g in three)
syl2 = set()
for i in range(len(invs)):
    for j in range(len(invs)):
        if i != j and compose(invs[i], invs[j]) == compose(invs[j], invs[i]):
            syl2.add(subgroup_from([invs[i], invs[j]]))
print(f"  n_5 = #Sylow-5 subgroups = {len(syl5)}  (icosa 5-fold axes = 6)  "
      f"{'OK' if len(syl5)==6 else '??'}")
print(f"  n_3 = #Sylow-3 subgroups = {len(syl3)}  (icosa 3-fold axes = 10) "
      f"{'OK' if len(syl3)==10 else '??'}")
print(f"  n_2 = #Sylow-2 (V4) subgroups = {len(syl2)}  (icosa 2-fold axes = 15 = #involutions, NOT n_2)")
print(f"  #involutions = {len(invs)}  (icosa 2-fold axes = 15) {'OK' if len(invs)==15 else '??'}")
print("  => {6,10,15} = {n_5, n_3, #involutions}; cyclic Sylows (p=3,5): axes = #subgroups;")
print("     p=2 deviates because Sylow-2 = V4 is non-cyclic (5 subgroups x 3 involutions = 15).")
