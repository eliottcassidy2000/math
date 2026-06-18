#!/usr/bin/env python3
"""
ADVERSARIAL REFUTATION of the M5 half-grid parity-XOR forbidden-class claim.

CLAIM (theme pinning-tightlocus, method M5 half-grid parity-XOR):
  Map: vertices = n runners with speeds v. Fix half-unit set {1,3,5}.
  For i<j, arc i->j iff  #{a in {1,3,5}: (v_i*a)%14 < (v_j*a)%14}  is ODD,
  else j->i. (GF(2)/XOR aggregation of three half-grid section orders.)
  CLAIMED FORBIDDEN at N=14: the regular T_5 score (2,2,2,2,2), c3=5, H=15.
  i.e. NO primitive 5-speed LRC input realizes the regular tournament.

We try to REFUTE: realize the regular T_5 iso class even once with this EXACT
map over a much broader exact search of LRC-constrained inputs:
  - all 5-subsets of large speed ranges (vmax up to several hundred)
  - 5-subsets drawn from covering / tight / sporadic 13-speed sets
  - residues mod 14 varied freely (incl. multiples of 14, all unit residues)
  - random large primitive 5-sets, huge speeds
  - sanity: confirm map is order-independent (relabel-invariant) so "iso class"
    is well defined, and confirm ties (v_i==v_j mod 14) handling.

All arithmetic is exact integer. No floats used for any decision.

instance: kps-S2-wf
"""
from math import gcd
from itertools import combinations, permutations
from functools import reduce
import random, sys, functools
print=functools.partial(print, flush=True)

HALF = (1, 3, 5)
N = 14

def is_primitive(S):
    return reduce(gcd, S) == 1

# ---- the M5 map, defined on an ARBITRARY ordered tuple of speeds ----
def m5_adj(speeds):
    """Return kxk adjacency (adj[i][j]=True means i->j), or None if not a
    valid tournament (a tie occurs: v_i==v_j mod 14 with c2 even/odd ambiguity).
    Orientation rule: arc i->j iff #{a in HALF: (v_i a)%N < (v_j a)%N} is ODD.
    NB: this is ANTISYMMETRIC by construction PROVIDED no pair ties in all three
    comparisons in a way that makes both directions even/odd. We detect ties."""
    k = len(speeds)
    adj = [[None]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            # count a where (v_i a)%N < (v_j a)%N
            cij = sum(1 for a in HALF if (speeds[i]*a) % N < (speeds[j]*a) % N)
            cji = sum(1 for a in HALF if (speeds[j]*a) % N < (speeds[i]*a) % N)
            # if equal in some a, that a contributes to neither cij nor cji.
            oij = (cij % 2 == 1)
            oji = (cji % 2 == 1)
            if oij == oji:
                # ambiguous: either both claim the arc or neither does.
                # The published code uses the i<j convention: arc i->j iff cij odd,
                # else j->i. Replicate that exactly (deterministic) but FLAG it.
                if cij % 2 == 1:
                    adj[i][j] = True; adj[j][i] = False
                else:
                    adj[i][j] = False; adj[j][i] = True
            else:
                if oij:
                    adj[i][j] = True; adj[j][i] = False
                else:
                    adj[i][j] = False; adj[j][i] = True
    for i in range(k):
        adj[i][i] = False
    return adj

def m5_adj_strict(speeds):
    """Strict replica of the published m5(): arc i->j iff cij ODD else j->i,
    using sorted-index convention. Returns adj or None if degenerate
    (some pair has v_i==v_j mod 14 making the comparison vacuous => still
    deterministic by the rule, but we want to know)."""
    k = len(speeds)
    adj = [[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            c2 = sum(1 for a in HALF if (speeds[i]*a) % N < (speeds[j]*a) % N)
            if c2 % 2 == 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
    # validity: tournament iff adj[i][j]!=adj[j][i] for all i<j (always true here
    # since we set exactly one). But check no self-contradiction.
    return adj

def is_tournament(adj, k):
    for i in range(k):
        for j in range(i+1, k):
            if adj[i][j] == adj[j][i]:
                return False
    return True

def score(adj, k):
    return tuple(sorted(sum(1 for j in range(k) if j != i and adj[i][j]) for i in range(k)))

def num3(adj, k):
    c = 0
    for i, j, l in combinations(range(k), 3):
        if adj[i][j] and adj[j][l] and adj[l][i]:
            c += 1
        if adj[i][l] and adj[l][j] and adj[j][i]:
            c += 1
    return c

# canonical iso-class signature for 5-vertex tournaments: (score, c3) is enough to
# separate all 12 classes EXCEPT we double-check by full canonical form for the
# regular target.
def is_regular_T5(adj):
    return is_tournament(adj, 5) and score(adj, 5) == (2,2,2,2,2)

def canon5(adj):
    """Full canonical adjacency under relabeling (brute over 120 perms),
    as a frozenset of directed edges on canonical labels -> a hashable key."""
    best = None
    for p in permutations(range(5)):
        edges = frozenset((p[i], p[j]) for i in range(5) for j in range(5)
                          if i != j and adj[i][j])
        key = tuple(sorted(edges))
        if best is None or key < best:
            best = key
    return best

# the unique regular T_5 canonical key (build one explicitly: rotational tournament
# on Z_5 with i->i+1,i+2)
def regular_T5_canon():
    adj = [[False]*5 for _ in range(5)]
    for i in range(5):
        adj[i][(i+1) % 5] = True
        adj[i][(i+2) % 5] = True
    return canon5(adj)

REG_KEY = regular_T5_canon()

# FAST EXACT test: among the 12 iso classes of 5-vertex tournaments, the regular
# tournament is the UNIQUE one with score (2,2,2,2,2). So score==(2,2,2,2,2) is a
# complete, exact identifier (no canon needed). We confirm canon agrees once.
def is_reg_fast(adj):
    return is_tournament(adj, 5) and score(adj, 5) == (2,2,2,2,2)
# sanity: among ALL 12 iso classes of 5-vertex tournaments, EXACTLY ONE has score
# (2,2,2,2,2) and it is the regular T5 (c3=5). Verified by full enumeration here.
_pairs = list(combinations(range(5),2))
_seen = {}
for _mask in range(1024):
    _a = [[False]*5 for _ in range(5)]
    for _b,(_i,_j) in enumerate(_pairs):
        if _mask>>_b & 1: _a[_i][_j]=True
        else: _a[_j][_i]=True
    _k = canon5(_a)
    if _k not in _seen:
        _seen[_k] = (score(_a,5), num3(_a,5))
_regs = [c for c in _seen if _seen[c][0]==(2,2,2,2,2)]
assert len(_seen)==12, f"expected 12 iso classes, got {len(_seen)}"
assert len(_regs)==1 and _seen[_regs[0]]==((2,2,2,2,2),5), "score (2,2,2,2,2) not unique-regular!"
assert _regs[0]==REG_KEY, "regular reference mismatch"
print(f"[sanity] 12 iso classes; score (2,2,2,2,2) is UNIQUELY the regular T5 (c3=5): CONFIRMED")

# ===================================================================
# SANITY 0: order-independence of the iso class produced by m5_adj_strict.
# If the map's iso class depends on input ordering, "forbidden class" is
# ill-posed. We test on random 5-sets whether all 120 orderings give the
# same iso class.
# ===================================================================
def iso_key_of_set(speedset):
    """Apply strict map to a given ordering and return canonical key + valid."""
    adj = m5_adj_strict(list(speedset))
    if not is_tournament(adj, len(speedset)):
        return None
    return canon5(adj)

print("="*72)
print("SANITY 0: is the m5 iso class independent of input ordering?")
print("="*72)
random.seed(1)
ord_dep = 0
ord_tested = 0
for _ in range(300):
    base = sorted(random.sample(range(1, 200), 5))
    if not is_primitive(base):
        continue
    keys = set()
    for p in permutations(base):
        adj = m5_adj_strict(list(p))
        if not is_tournament(adj, 5):
            keys.add(("INVALID",))
            continue
        keys.add(canon5(adj))
    ord_tested += 1
    if len(keys) > 1:
        ord_dep += 1
        if ord_dep <= 3:
            print(f"  ORDER-DEPENDENT set {base}: {len(keys)} distinct iso classes")
print(f"  tested {ord_tested} primitive 5-sets; order-dependent: {ord_dep}")
print(f"  (if >0, the map depends on labeling; published claim uses SORTED order)")
print()

# Going forward we use SORTED order (matching the published evidence which uses
# combinations -> sorted tuples).

def realizes_regular(speeds_sorted):
    adj = m5_adj_strict(list(speeds_sorted))
    if not is_tournament(adj, len(speeds_sorted)):
        return False
    return canon5(adj) == REG_KEY

# ===================================================================
# ATTACK 1: brute force ALL primitive sorted 5-sets up to large vmax.
# Published evidence: vmax=40 -> NOT realized. We push much further.
# ===================================================================
print("="*72)
print("ATTACK 1: exhaustive primitive sorted 5-sets, escalating vmax")
print("="*72)
witnesses = []
for vmax in (40, 55):
    found = None
    count = 0
    classes = set()
    for S in combinations(range(1, vmax+1), 5):
        if not is_primitive(S):
            continue
        adj = m5_adj_strict(list(S))
        if not is_tournament(adj, 5):
            continue
        count += 1
        classes.add((score(adj,5), num3(adj,5)))
        if is_reg_fast(adj):
            found = S
            witnesses.append(("attack1", S))
            break
    print(f"  vmax={vmax:3d}: {count} valid tournaments, "
          f"{len(classes)} distinct (score,c3) classes, regular T5 found={found}")
    if found:
        break
print()

# ===================================================================
# ATTACK 2: residues mod 14 are what actually matter for the comparisons
# (v*a)%14 depends only on v mod 14. So the iso class depends ONLY on the
# multiset of residues (v_i mod 14) AND the sorted ORDER of the actual speeds.
# To break the claim we should sweep ALL residue tuples mod 14 in all sorted
# realizations. Enumerate every sorted 5-set of residues r_i in 0..13 and ask:
# does ANY assignment of distinct increasing integer speeds with those residues
# realize the regular class? Since the map uses (v_i mod 14) for comparisons but
# index order = sorted speed order, the iso class depends on (residue-in-sorted-
# -speed-order). We can realize ANY ordering of residues by choosing speeds, as
# long as residues need not be distinct. So: enumerate ALL ordered residue
# 5-tuples (14^5) and test the map directly on the residues as the "speeds"
# in that fixed index order. This is the COMPLETE residue-level search.
# ===================================================================
print("="*72)
print("ATTACK 2: COMPLETE residue-level search (all 14^5 ordered residue tuples)")
print("  The map output depends only on residues mod 14 in index order.")
print("  Any ordered residue tuple is realizable by integer speeds (add 14*k).")
print("="*72)
from itertools import product
total = 0
valid = 0
reg_hits = []
classes_seen = set()
for r in product(range(14), repeat=5):
    total += 1
    adj = m5_adj_strict(list(r))
    if not is_tournament(adj, 5):
        continue
    valid += 1
    sc = score(adj, 5); c3 = num3(adj, 5)
    classes_seen.add((sc, c3))
    if is_reg_fast(adj):
        reg_hits.append(r)
print(f"  total ordered residue tuples: {total}")
print(f"  valid tournaments: {valid}")
print(f"  distinct (score,c3) classes realized: {len(classes_seen)}")
print(f"  regular T5 realizations: {len(reg_hits)}")
if reg_hits:
    print(f"  FIRST regular-T5 residue witness: {reg_hits[0]}")
print()

# But ATTACK 2 ignores the constraint that index order = SORTED SPEED order.
# In a real LRC input the speeds are positive integers; if residues are distinct
# we can pick speeds to realize ANY ordering of residues (e.g. residue order
# need not match value order: speed 15 (res 1) > speed 3 (res 3) means a smaller
# residue at a larger index). So actually arbitrary ordered residue tuples ARE
# realizable -> ATTACK 2 IS the faithful complete search. We make this rigorous:

print("="*72)
print("ATTACK 2b: confirm every ordered residue tuple is LRC-realizable by speeds")
print("  Construct explicit primitive integer 5-sets for any regular-T5 residue hit")
print("="*72)
def realize_residue_order(r):
    """Given ordered residues r[0..4], produce a sorted integer 5-set whose
    SORTED order has those residues in order, i.e. speeds s0<s1<...<s4 with
    s_t % 14 == r[t]. Choose s_t = r[t] + 14*t (+offset to keep positive &
    strictly increasing). r[t]+14t is strictly increasing in t since 14>13>=r."""
    speeds = [r[t] + 14*t for t in range(5)]
    # ensure all positive & strictly increasing & distinct
    # r[t]+14t : for t<t', r[t']+14t' - (r[t]+14t) = 14(t'-t)+(r[t']-r[t]) >= 14-13>0
    speeds = [s if s > 0 else s+14 for s in speeds]
    # make primitive
    g = reduce(gcd, speeds)
    if g > 1:
        speeds = [s//g for s in speeds]  # note: dividing changes residues! so instead:
    return speeds

# Better: verify by direct construction that the map on the constructed speeds
# (sorted) reproduces the regular class for each residue hit.
confirmed = 0
example_sets = []
for r in reg_hits[:50]:
    speeds = [r[t] + 14*t for t in range(5)]
    # speeds strictly increasing already; residues match r exactly
    assert speeds == sorted(speeds)
    assert all(speeds[t] % 14 == r[t] for t in range(5))
    if realizes_regular(tuple(speeds)):
        confirmed += 1
        if len(example_sets) < 5:
            example_sets.append((r, speeds, is_primitive(speeds)))
print(f"  residue hits checked: {min(len(reg_hits),50)}; "
      f"reproduced as integer-speed sets: {confirmed}")
for r, sp, prim in example_sets:
    print(f"    residues {r} -> speeds {sp}  primitive={prim}")
print()

# ===================================================================
# ATTACK 3: but are these LRC-VALID inputs? LRC wants PRIMITIVE sets and the
# forbidden-class claim is about ANY primitive 5-speed set. Find an explicit
# PRIMITIVE integer 5-set realizing the regular class (if any residue hit exists).
# ===================================================================
print("="*72)
print("ATTACK 3: explicit PRIMITIVE integer 5-sets realizing regular T5 (if any)")
print("="*72)
prim_witnesses = []
if reg_hits:
    # search small primitive sets whose residue tuple (in sorted order) is a hit
    reg_hit_set = set(reg_hits)
    for vmax in (100, 200, 400):
        for S in combinations(range(1, vmax+1), 5):
            if not is_primitive(S):
                continue
            r = tuple(v % 14 for v in S)  # S sorted -> residues in sorted-speed order
            if r in reg_hit_set:
                # double check the actual map
                if realizes_regular(S):
                    prim_witnesses.append(S)
                    if len(prim_witnesses) >= 10:
                        break
        if prim_witnesses:
            print(f"  vmax={vmax}: found {len(prim_witnesses)} primitive witnesses")
            for w in prim_witnesses[:10]:
                adj = m5_adj_strict(list(w))
                print(f"    S={w} residues={tuple(v%14 for v in w)} "
                      f"score={score(adj,5)} c3={num3(adj,5)}")
            break
        else:
            print(f"  vmax={vmax}: no primitive witness yet among residue hits")
else:
    print("  no residue hits => no primitive witness possible at residue level")
print()

# ===================================================================
# ATTACK 4: 5-subsets of canonical LRC 13-speed sets (covering, tight, sporadic)
# ===================================================================
print("="*72)
print("ATTACK 4: 5-subsets drawn from canonical LRC 13-speed sets")
print("="*72)
families = {
    "AP {1..13}": list(range(1, 14)),
    "Goddyn-Wong {1..11,13,24}": [1,2,3,4,5,6,7,8,9,10,11,13,24],
    "covering {1..11,13,84}": [1,2,3,4,5,6,7,8,9,10,11,13,84],
    "covering {1..12,84}": list(range(1,13))+[84],  # has multiple of 14? 84=14*6 yes
    "AP-drop-6 core+w {1,2,3,4,5,7,8,9,10,11,12,13,84}": [1,2,3,4,5,7,8,9,10,11,12,13,84],
}
for name, big in families.items():
    hit = None
    cnt = 0
    cls = set()
    for S in combinations(sorted(big), 5):
        if not is_primitive(S):
            continue
        adj = m5_adj_strict(list(S))
        if not is_tournament(adj, 5):
            continue
        cnt += 1
        cls.add((score(adj,5), num3(adj,5)))
        if is_reg_fast(adj):
            hit = S
            break
    print(f"  {name}: {cnt} 5-subsets, {len(cls)} classes, regular T5={hit}")
print()

# ===================================================================
# ATTACK 5: random huge primitive 5-sets (sporadic, large speeds)
# ===================================================================
print("="*72)
print("ATTACK 5: random large primitive 5-sets (sporadic)")
print("="*72)
random.seed(12345)
huge_hit = None
trials = 200000
seen_cls = set()
for _ in range(trials):
    S = tuple(sorted(random.sample(range(1, 100000), 5)))
    if not is_primitive(S):
        continue
    adj = m5_adj_strict(list(S))
    if not is_tournament(adj, 5):
        continue
    seen_cls.add((score(adj,5), num3(adj,5)))
    if is_reg_fast(adj):
        huge_hit = S
        break
print(f"  {trials} random trials, {len(seen_cls)} (score,c3) classes seen, "
      f"regular T5 found={huge_hit}")
print()

# ===================================================================
# VERDICT
# ===================================================================
print("="*72)
print("VERDICT SUMMARY")
print("="*72)
print(f"  REG_KEY (canonical regular T5) = {REG_KEY}")
print(f"  ATTACK1 witnesses: {[w for t,w in witnesses if t=='attack1']}")
print(f"  ATTACK2 residue-level regular-T5 realizations: {len(reg_hits)}")
print(f"  ATTACK3 primitive integer witnesses: {prim_witnesses[:5]}")
print(f"  ATTACK5 huge random witness: {huge_hit}")
refuted = bool(witnesses or reg_hits or prim_witnesses or huge_hit)
print()
print(f"  CLAIM REFUTED: {refuted}")
if reg_hits and not prim_witnesses:
    print("  NOTE: residue-level realization exists but check primitive realizability.")
print("="*72)
