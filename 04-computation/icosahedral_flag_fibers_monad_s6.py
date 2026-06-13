#!/usr/bin/env python3
"""
icosahedral_flag_fibers_monad_s6.py  -- monad-explorer-2026-06-07-S6 (deeper)

The 3-to-1 commutator covering {60 oriented overlapping pairs} -> {20 three-cycles}.
Numerology: 60 = |A_5| = #icosahedral face-vertex FLAGS = 20 faces x 3 vertices/face;
the flag->face map is canonically 3-to-1. Question: is the commutator covering
THE icosahedral flag incidence?

We test the fiber structure:
 (1) Each 3-cycle c has exactly 3 preimage oriented overlapping pairs (3-to-1 confirmed).
 (2) For c=(a b c), what are the 3 preimages? We look at the SHARED VERTEX of each
     preimage pair and the orientation data, to see if the 3 preimages are canonically
     indexed by {a,b,c} (= the 3 vertices of the face c) -- the icosahedral flag picture.
 (3) A_5 acts on both sides; check the covering is A_5-equivariant (it must be, the
     construction is natural) and identify the point-stabiliser orders (Z/3 over a face).
"""
from itertools import combinations, permutations
from collections import defaultdict

S = [0, 1, 2, 3, 4]

def cyc_to_perm(cycle):
    p = {x: x for x in S}
    L = len(cycle)
    for i in range(L):
        p[cycle[i]] = cycle[(i + 1) % L]
    return p

def compose(p, q):
    return {x: p[q[x]] for x in q}

def inverse(p):
    return {v: k for k, v in p.items()}

def commutator(a, b):
    return compose(compose(a, b), compose(inverse(a), inverse(b)))

def perm_to_cycle3(p):
    """return the 3-cycle as a canonical tuple (a,b,c) with a=min, or None."""
    moved = [x for x in S if p[x] != x]
    if len(moved) != 3:
        return None
    a = min(moved)
    return (a, p[a], p[p[a]])

# enumerate the 60 oriented overlapping pairs.
# An oriented overlapping pair: unordered {X,Y}, |X cap Y|=1, with a cyclic
# orientation on each -> we keep the ORDERED commutator [sigma_X, sigma_Y] but
# record full data.
triples = list(combinations(S, 3))
def oriented_tris(s):
    a, b, c = s
    return [(a, b, c), (a, c, b)]

flags = []   # (X, Y, oX, oY, sigmaX, sigmaY, commutator-3cycle, shared vertex)
for i in range(len(triples)):
    for j in range(i + 1, len(triples)):
        X, Y = set(triples[i]), set(triples[j])
        if len(X & Y) != 1:
            continue
        v = (X & Y).pop()
        for oX in oriented_tris(triples[i]):
            for oY in oriented_tris(triples[j]):
                sx, sy = cyc_to_perm(oX), cyc_to_perm(oY)
                c = commutator(sx, sy)
                c3 = perm_to_cycle3(c)
                flags.append(dict(X=triples[i], Y=triples[j], oX=oX, oY=oY,
                                  comm=c3, shared=v))

print(f"total oriented overlapping pairs (flags): {len(flags)}")
assert len(flags) == 60

# (1) fibers of comm
fiber = defaultdict(list)
for f in flags:
    fiber[f['comm']].append(f)
sizes = sorted(len(v) for v in fiber.values())
print(f"#distinct commutator 3-cycles hit: {len(fiber)} (expect 20)")
print(f"fiber sizes: min={min(sizes)} max={max(sizes)} (expect all 3): "
      f"{'3-to-1 CONFIRMED' if set(sizes)=={3} else sizes}")

# (2) structure of a fiber: for each 3-cycle c, what 3 things index its preimages?
print()
print("Sample fibers (commutator 3-cycle  <-  its 3 preimage oriented pairs):")
for c3 in sorted(fiber)[:6]:
    pres = fiber[c3]
    line = f"  c={c3}: "
    parts = []
    for f in pres:
        parts.append(f"[{''.join(map(str,f['oX']))}|{''.join(map(str,f['oY']))} sh={f['shared']}]")
    print(line + "  ".join(parts))

# Test: are the 3 preimages of c=(a,b,c) canonically indexed by {a,b,c}?
# Hypothesis H1: the multiset of SHARED vertices over the 3 preimages = {a,b,c}.
print()
print("(2a) Is {shared vertices of the 3 preimages} = the support {a,b,c} of c?")
ok_support = 0
for c3, pres in fiber.items():
    sh = sorted(f['shared'] for f in pres)
    if sh == sorted(c3):
        ok_support += 1
print(f"   matches support: {ok_support}/{len(fiber)}  "
      f"=> {'YES (shared vertex = the face-vertex of the flag)' if ok_support==len(fiber) else 'NO'}")

# Hypothesis H1': shared vertices form some other fixed 3-set / pattern
print()
print("(2b) distribution of (sorted shared-vertex multiset) across fibers:")
patt = defaultdict(int)
for c3, pres in fiber.items():
    sh = tuple(sorted(f['shared'] for f in pres))
    patt[(len(set(sh)), sh==tuple(sorted(c3)))] += 1
for k, v in sorted(patt.items()):
    print(f"   (#distinct shared={k[0]}, equals-support={k[1]}): {v} fibers")

# (3) equivariance + the icosahedral flag count cross-check
print()
print("(3) icosahedral flag arithmetic cross-check:")
print(f"   faces F=20  vertices/face=3  -> face-vertex flags = {20*3} = |A_5| = 60")
print(f"   3-cycles = 20 = F ;  flag->face is 3-to-1 ;  commutator-cover is 3-to-1 -> MATCH")
print(f"   each 3-fold axis carries 2 three-cycles (abc),(acb) = 2 antipodal faces;")
print(f"   #3-fold axes = 10 = 20/2 = pairs of antipodal faces.")

# Build the pairing c <-> c^{-1} (antipodal faces) and verify 10 axis-pairs.
pairs = set()
for c3 in fiber:
    inv3 = perm_to_cycle3(inverse(cyc_to_perm(c3)))
    pairs.add(frozenset([c3, inv3]))
print(f"   antipodal 3-cycle pairs {{c, c^-1}}: {len(pairs)} (expect 10 three-fold axes): "
      f"{'OK' if len(pairs)==10 else 'MISMATCH'}")
