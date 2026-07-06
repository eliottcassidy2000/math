#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S35: ADDITIVE ISOLATION = the hidden link between
mac-mini's kissing number (HYP-4542), opus's height-independence (HYP-4476),
and my boundary-seating picture (S34).

CLAIM (small, checkable): in a gap member = base AP + defects, the FAR OUTLIER is
the runner that participates in NO additive triple (v_i + v_j = v_l) -- it is
"additively isolated", contributing 0 to the relation lattice's kissing number
(mac-mini HYP-4542: kissing = # additive triples).  Two consequences we test:

  (A) BINDING: the additively-isolated runner is always among the M-binding
      runners (achieves the min residue-distance) at the witness t.
  (B) HEIGHT-INDEPENDENCE (opus HYP-4476, tested structurally): replacing the
      isolated outlier by ANY other additively-isolated value leaves M unchanged.
      => M depends on the DEFECT PATTERN, not the height, BECAUSE the outlier is
      additively isolated (no relation to constrain it, so only its residue counts,
      and a generic isolated residue lands at the same boundary seat).

If both hold, the synthesis is: additive isolation (kissing=0) => height-independence
(opus) => the outlier sits at a boundary seat (my S34).  One fact, three faces.
"""
from fractions import Fraction
from itertools import combinations
import numpy as np

def M_and_witness(v):
    """Return (M as Fraction, best (c,q), residue-distances at best)."""
    v = [x for x in v if x != 0]
    S = int(sum(abs(x) for x in v)); Q = min(4*S, 2*max(abs(x) for x in v)+2)
    va = np.array(v, dtype=np.int64); bn, bd, bc = 0, 1, 1
    for q in range(2, Q+1):
        for c in range(1, q):
            r = (va*c) % q; d = np.minimum(r, q-r); bq = int(d.min())
            if bq*bd > bn*q: bn, bd, bc = bq, q, c
    dists = [int(min((x*bc) % bd, bd-(x*bc) % bd)) for x in v]
    return Fraction(bn, bd), (bc, bd), dists

def additive_triples(v):
    s = set(v); tri = []
    for i in range(len(v)):
        for j in range(i, len(v)):
            if v[i]+v[j] in s and v[i]+v[j] not in (v[i], v[j]):
                tri.append((v[i], v[j], v[i]+v[j]))
    return tri

def isolated(v):
    """runners in NO additive triple (as summand or sum)."""
    tri = additive_triples(v)
    used = set()
    for a, b, c in tri:
        used.add(a); used.add(b); used.add(c)
    return [x for x in v if x not in used], tri

MEMBERS = {
    "n=8 A": [1,2,3,4,5,7,18],
    "n=8 B": [1,3,4,5,7,13,18],
    "n=7 (5/33)": None,  # filled below by search
}

# find an n=7 gap member (M in (1/7,2/13)) by quick search among base-AP+outlier
def find_n7():
    lo, hi = Fraction(1,7), Fraction(2,13)
    for base_d in (1,2):
        base = [1+i*base_d for i in range(5)]
        for out in combinations(range(base[-1]+1, 34), 1):
            v = sorted(set(base+list(out)))
            if len(v) != 6: continue
            import math
            g = 0
            for x in v: g = math.gcd(g, x)
            if g != 1: continue
            M,_,_ = M_and_witness(v)
            if lo < M < hi:
                return v, M
    return None, None

v7, M7 = find_n7()
if v7: MEMBERS["n=7 (5/33)"] = v7

print("=== ADDITIVE ISOLATION on the gap members (kps-S35) ===\n", flush=True)
for name, v in MEMBERS.items():
    if v is None:
        print(f"{name}: (no member found)\n"); continue
    M, (c,q), dists = M_and_witness(v)
    iso, tri = isolated(v)
    mind = min(dists)
    binders = [v[i] for i in range(len(v)) if dists[i] == mind]
    print(f"{name}: v={v}", flush=True)
    print(f"  M={M}  witness t={c}/{q}  residue-dists={dict(zip(v,dists))}", flush=True)
    print(f"  additive triples ({len(tri)}): {tri}", flush=True)
    print(f"  ADDITIVELY ISOLATED (kissing-0): {iso}", flush=True)
    print(f"  BINDING (min dist {mind}): {binders}", flush=True)
    print(f"  (A) isolated subset of binders? {set(iso).issubset(set(binders))}", flush=True)
    print(flush=True)

# (B) HEIGHT-INDEPENDENCE test: vary the outlier of n=8 A over isolated vs non-isolated values
print("=== (B) HEIGHT-INDEPENDENCE: {1,2,3,4,5,7,x}, x=13..40 -- M vs whether x is isolated ===", flush=True)
base = [1,2,3,4,5,7]
gapM = Fraction(3,23)
same = []; diff = []
for x in range(13, 41):
    v = sorted(base+[x])
    if len(set(v)) != 7: continue
    M,_,_ = M_and_witness(v)
    iso,_ = isolated(v)
    is_iso = x in iso
    tag = "ISO" if is_iso else "   "
    mark = "== 3/23" if M == gapM else "!= 3/23"
    print(f"  x={x:>2} {tag}  M={str(M):>7}  {mark}", flush=True)
    (same if M == gapM else diff).append((x, is_iso))
iso_same = [x for x,i in same if i]; noniso_same = [x for x,i in same if not i]
print(f"\n  x giving M=3/23: {[x for x,_ in same]}", flush=True)
print(f"    of these, additively isolated: {iso_same}", flush=True)
print(f"    non-isolated: {noniso_same}", flush=True)
print(f"  x giving M!=3/23: {[x for x,_ in diff]}", flush=True)
print("\nREADING: if 'M=3/23' coincides with 'x additively isolated', then additive", flush=True)
print("isolation (kissing=0) IS opus's height-independence, and the outlier's only", flush=True)
print("role is to occupy a boundary seat (my S34). One fact, three faces.", flush=True)
