#!/usr/bin/env python3
"""
LRC(14) creative tournament-generation, theme = "pinning-tightlocus".

We build several DIFFERENT arc-definition maps that take LRC data
(speeds S, lonely time tau, sections, binding/pinning structure) and
produce a tournament on a small vertex set. For each we:
  (b) check NON-TRIVIALITY (does it ever produce a 3-cycle, H>1?)
  (c) enumerate realized iso classes vs the FULL free set (A000568)
  (d) report forbidden classes with search extent.

EXACT rationals only. M tool copied verbatim from the task spec.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

# ---------- EXACT M TOOL (verbatim) ----------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t):
    return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = g(S, t)
        if v > b:
            b = v; at = t
    return b, at
# All optima (there may be several tau achieving the max)
def M_all(S):
    b, _ = M(S)
    ats = sorted(t for t in cand(S) if g(S, t) == b)
    return b, ats

# ---------- tournament canonicalization ----------
# A tournament on k labeled vertices: adj[i][j] = True iff arc i->j.
# We need iso-class canonical form. Use brute force over k! relabelings
# of the 0/1 adjacency upper-triangle bit-string; canonical = min bitstring.

def tour_key(adj, k):
    """Canonical key = min over all relabelings of the edge-bit integer."""
    best = None
    verts = list(range(k))
    for perm in permutations(verts):
        bits = 0; idx = 0
        for i in range(k):
            for j in range(k):
                if i == j: continue
                # arc perm-image: original edge perm[i]->perm[j]
                bit = 1 if adj[perm[i]][perm[j]] else 0
                bits = (bits << 1) | bit
        if best is None or bits < best:
            best = bits
    return best

def is_tournament(adj, k):
    for i in range(k):
        for j in range(i+1, k):
            if adj[i][j] == adj[j][i]:
                return False  # need exactly one arc
    return True

def num_3cycles(adj, k):
    c = 0
    for i, j, l in combinations(range(k), 3):
        # count cyclic triangles among the three
        # gather arcs
        def arc(a, b): return adj[a][b]
        # cyclic if i->j->l->i or i->l->j->i
        if arc(i,j) and arc(j,l) and arc(l,i): c += 1
        if arc(i,l) and arc(l,j) and arc(j,i): c += 1
    return c

def num_hampaths(adj, k):
    """H(T) = number of Hamiltonian paths (directed)."""
    h = 0
    for perm in permutations(range(k)):
        ok = True
        for a in range(k-1):
            if not adj[perm[a]][perm[a+1]]:
                ok = False; break
        if ok: h += 1
    return h

def score_seq(adj, k):
    return tuple(sorted(sum(1 for j in range(k) if i != j and adj[i][j]) for i in range(k)))

# A000568 full free set counts
A000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456}

def enumerate_free_set(k):
    """All iso classes of tournaments on k vertices: dict canonkey -> rep info."""
    classes = {}
    # iterate over all orientations of C(k,2) edges
    edges = list(combinations(range(k), 2))
    m = len(edges)
    for bits in range(1 << m):
        adj = [[False]*k for _ in range(k)]
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
        key = tour_key(adj, k)
        if key not in classes:
            classes[key] = {
                'score': score_seq(adj, k),
                'c3': num_3cycles(adj, k),
                'H': num_hampaths(adj, k),
            }
    return classes

FREE = {}
for k in (3, 4, 5):
    FREE[k] = enumerate_free_set(k)
    assert len(FREE[k]) == A000568[k], f"free set k={k}: {len(FREE[k])} != {A000568[k]}"

def class_label(k, key):
    info = FREE[k][key]
    return f"[score={info['score']},c3={info['c3']},H={info['H']}]"

# ================================================================
# Helper LRC machinery
# ================================================================
def is_primitive(S):
    from functools import reduce
    return reduce(gcd, S) == 1

def sections(S, a, N=14):
    """grid time tau=a/N; section r_i = v_i*a mod N."""
    return [ (v*a) % N for v in S ]

# binding runners at the optimum: runners i with ||v_i*tau|| == M
def binding_set(S, tau, gap):
    return [i for i, v in enumerate(S) if nrm(v*tau) == gap]

# frac position of runner at tau (in [0,1))
def fracpos(v, tau):
    r = v*tau - int(v*tau)
    return r + 1 if r < 0 else r

# ================================================================
# METHODS
# Each method: given an LRC input, produce (k, adj) or None (if undefined).
# We then range over a constrained family of inputs and bucket iso classes.
# ================================================================

# ----------------------------------------------------------------
# Input family generators
# ----------------------------------------------------------------
def small_speed_sets(nspeed, vmax):
    """All primitive sorted speed sets of size nspeed with speeds in 1..vmax."""
    out = []
    for S in combinations(range(1, vmax+1), nspeed):
        if is_primitive(S):
            out.append(S)
    return out

# ================================================================
# METHOD 1: PLATEAU-PINNING TOURNAMENT (pairs as vertices)
# ----------------------------------------------------------------
# The tight locus {1..13} pins gap 1/14 at tau in {1/14,3/14,5/14},
# each plateau bound by a disjoint extreme PAIR. Generalize:
# VERTICES = the binding PAIRS {a,b} of speeds that co-pin some optimum
# plateau. We need a tournament so we pick a fixed small set of candidate
# pairs and orient by WHICH pin (tau value) each pair binds first.
#
# Cleaner concrete version (Method 1): VERTICES = runners that are BINDING
# at the optimum tau. Arc i->j iff frac(v_i*tau) < frac(v_j*tau) AND they
# are on OPPOSITE sides (one near 0+gap, other near 1-gap). This is the
# "endpoint protection" order: who hugs the lower edge vs upper edge.
# But a single tau snapshot frac-order is transitive -> likely dead.
# We instead orient by the SIGNED pin: i->j iff runner i binds from BELOW
# (frac = gap) and j binds from ABOVE (frac = 1-gap)  -> bipartite-ish.
# Test non-triviality.
# ================================================================
def method1(S):
    """Vertices = binding runners at THE (first) optimum tau.
    Arc i->j iff fracpos(i) < fracpos(j) (snapshot order among binders).
    SNAPSHOT order -> expected transitive. Included to confirm dead."""
    gap, ats = M_all(S)
    tau = ats[0]
    B = binding_set(S, tau, gap)
    k = len(B)
    if k < 3: return None
    adj = [[False]*k for _ in range(k)]
    fr = [fracpos(S[i], tau) for i in B]
    for a in range(k):
        for b in range(k):
            if a == b: continue
            if fr[a] < fr[b]:
                adj[a][b] = True
    if not is_tournament(adj, k): return None
    return k, adj

# ================================================================
# METHOD 2: CO-PINNING INCIDENCE / CYCLIC PLATEAU ORDER
# ----------------------------------------------------------------
# VERTICES = a chosen set of speeds (runners). For an ordered pair (i,j),
# look at ALL optima tau (the plateaus). Arc i->j iff, scanning the sorted
# list of optimal taus, the FIRST plateau where exactly one of {i,j} is
# binding has i binding (i pins earlier). If both/neither always co-occur,
# break by speed. This uses the MULTI-plateau pin structure (not a single
# snapshot) -> can be cyclic.
# ================================================================
def method2(S):
    gap, ats = M_all(S)
    k = len(S)
    if k < 3: return None
    # binding indicator per tau
    bind = [ set(binding_set(S, t, gap)) for t in ats ]
    adj = [[False]*k for _ in range(k)]
    ok = True
    for i in range(k):
        for j in range(k):
            if i == j: continue
            decided = False
            for t_idx in range(len(ats)):
                bi = i in bind[t_idx]; bj = j in bind[t_idx]
                if bi != bj:
                    adj[i][j] = bi  # i pins at this earliest distinguishing plateau
                    decided = True
                    break
            if not decided:
                # tie -> orient by speed (smaller speed -> )
                adj[i][j] = S[i] < S[j]
    if not is_tournament(adj, k): return None
    return k, adj

# ================================================================
# METHOD 3: SECTION-CROSSING PARITY TOURNAMENT (on-grid sections)
# ----------------------------------------------------------------
# VERTICES = runners. Look at the grid times tau=a/14, a in (Z/14)*.
# Section r_i(a) = v_i*a mod 14. Arc i->j by the SIGNED SUM over a of
# sign(section_i - section_j) compared across the unit group, i.e. a
# "majority who sits higher in section order across the 6 grid columns".
# Majority votes can produce a Condorcet-cycle -> potentially cyclic.
# Ties / zero-sum broken by speed.
# ================================================================
UNITS14 = [a for a in range(1,14) if gcd(a,14)==1]  # {1,3,5,9,11,13}
def method3(S, N=14):
    k = len(S)
    if k < 3: return None
    adj = [[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            tally = 0
            for a in UNITS14:
                ri = (S[i]*a) % N
                rj = (S[j]*a) % N
                # compare distance-to-0 in the section circle (loneliness sense):
                # use centered residue |.| with sign by which is "more lonely"
                di = min(ri, N-ri); dj = min(rj, N-rj)
                if di > dj: tally += 1
                elif di < dj: tally -= 1
            if tally > 0: adj[i][j] = True
            elif tally < 0: adj[i][j] = False
            else: adj[i][j] = S[i] < S[j]
    if not is_tournament(adj, k): return None
    return k, adj

# ================================================================
# METHOD 4: PAIR-BINDING CYCLIC TOURNAMENT (pairs as vertices)
# ----------------------------------------------------------------
# The KEY pinning fact: each plateau is bound by a binding PAIR {a,b}
# with crossing tau = k/(v_a +- v_b). VERTICES = the binding pairs that
# achieve gap M (the pins). Arc P->Q (pairs P,Q) iff the pin-tau of P is
# LESS than the pin-tau of Q. With multiple co-optimal pins at the SAME
# gap, several pairs share the plateau-time structure; we orient by
# (pin-tau, then by the cyclic difference of section midpoints).
# To get a cyclic tournament we orient pairs by the CYCLIC order of their
# crossing taus on the circle R/Z (rotational), which is genuinely cyclic.
# ================================================================
def binding_pairs(S, tau, gap):
    """pairs (i,j) such that tau is a crossing k/(vi+-vj) AND both bind."""
    B = set(binding_set(S, tau, gap))
    pairs = []
    for i, j in combinations(range(len(S)), 2):
        if i in B and j in B:
            pairs.append((i, j))
    return pairs

def method4(S):
    gap, ats = M_all(S)
    # collect ALL binding pairs across all optimal plateaus, tag with tau
    pinned = []  # (pair, tau)
    seen = set()
    for t in ats:
        for p in binding_pairs(S, t, gap):
            if p not in seen:
                seen.add(p); pinned.append((p, t))
    k = len(pinned)
    if k < 3: return None
    if k > 7: return ('TOOBIG', k)  # skip canonicalization for big ones, but record
    adj = [[False]*k for _ in range(k)]
    # cyclic order of taus on circle: pick reference = min tau, orient by
    # rotational position; arc P->Q iff (tauQ - tauP) mod 1 in (0, 1/2)
    for a in range(k):
        for b in range(k):
            if a == b: continue
            ta = pinned[a][1]; tb = pinned[b][1]
            diff = (tb - ta) % 1
            if diff == 0:
                # same tau -> orient by section midpoint of the pair
                pa = pinned[a][0]; pb = pinned[b][0]
                ma = (S[pa[0]] + S[pa[1]]); mb = (S[pb[0]] + S[pb[1]])
                adj[a][b] = ma < mb
            elif diff < F(1,2):
                adj[a][b] = True
            else:
                adj[a][b] = False
    if not is_tournament(adj, k): return None
    return k, adj

# ================================================================
# METHOD 5: COMPLEMENT/REVERSAL-TWISTED SECTION TOURNAMENT
# ----------------------------------------------------------------
# a=13=-1 sends section r -> 14-r (reversal = T->T^op). VERTICES = runners.
# For each unit a, runner i has section r_i(a). Define arc i->j by the
# PARITY of the number of grid columns a in {1,3,5} (a half-set of units,
# the "first half" before the reversal a=13 mirrors them) where
# section_i(a) < section_j(a). This is a parity/XOR aggregation of
# half-grid orders -> NOT a single snapshot, can be cyclic.
# ================================================================
HALF_UNITS = [1, 3, 5]  # the half not including the -1 reversal partners
def method5(S, N=14):
    k = len(S)
    if k < 3: return None
    adj = [[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            cnt = 0
            for a in HALF_UNITS:
                ri = (S[i]*a) % N; rj = (S[j]*a) % N
                if ri < rj: cnt += 1
            # parity orientation
            if cnt % 2 == 1: adj[i][j] = True
            else: adj[i][j] = False
            # ensure antisymmetry: check i->j and j->i won't both be set
    # Fix: this parity rule may NOT be antisymmetric. Rebuild properly:
    adj = [[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            cnt = 0
            for a in HALF_UNITS:
                ri = (S[i]*a) % N; rj = (S[j]*a) % N
                if ri < rj: cnt += 1
                elif ri > rj: cnt -= 1
            # orient by sign of parity-weighted... use signed count parity:
            # arc by: if cnt>0 i->j else j->i, tie by speed; but that's a vote (method3-like).
            # Make it genuinely XOR: count columns where ri<rj, take parity.
            c2 = sum(1 for a in HALF_UNITS if (S[i]*a)%N < (S[j]*a)%N)
            if c2 % 2 == 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
    if not is_tournament(adj, k): return None
    return k, adj

# ================================================================
# METHOD 6: PIN-RANK TOURNAMENT (who reaches gap M at the smallest tau)
# ----------------------------------------------------------------
# VERTICES = runners. For each runner i define t_i = the SMALLEST tau in
# (0,1/2] at which ||v_i*tau|| == gap M(S) (first time runner i is pinned
# at the global gap). Arc i->j iff t_i < t_j (who gets pinned first as tau
# sweeps up). Ties broken by speed. This sweeps the WHOLE pin history, not
# a single snapshot. Each t_i = (2k+1)/(2 v_i)-ish; could be cyclic? A pure
# numeric order is transitive, BUT we instead use CYCLIC sweep order from a
# random phase... Actually pin-time is a real number -> transitive. We add a
# twist: arc i->j iff t_i < t_j measured on the CIRCLE from the optimum tau*
# (rotate so tau* = 0). Rotational cut can break transitivity.
# ================================================================
def first_pin_time(v, gap, N_den_max=200):
    """smallest tau in (0,1/2] with ||v*tau||==gap. gap is a Fraction.
    ||v*tau||=gap means v*tau = k +- gap. tau=(k+-gap)/v. smallest >0."""
    best = None
    # k from 0.. up to v
    for k in range(0, v+1):
        for s in (gap, -gap):
            tau = (F(k) + s) / v
            if F(0) < tau <= F(1,2):
                if best is None or tau < best:
                    best = tau
    return best

def method6(S):
    gap, ats = M_all(S)
    taustar = ats[0]
    k = len(S)
    if k < 3: return None
    t = []
    for v in S:
        fp = first_pin_time(v, gap)
        t.append(fp if fp is not None else F(1))
    adj = [[False]*k for _ in range(k)]
    # rotational: position on circle measured from taustar
    pos = [ (t[i] - taustar) % 1 for i in range(k) ]
    for i in range(k):
        for j in range(k):
            if i == j: continue
            d = (pos[j] - pos[i]) % 1
            if d == 0:
                adj[i][j] = S[i] < S[j]
            elif d < F(1,2):
                adj[i][j] = True
            else:
                adj[i][j] = False
    if not is_tournament(adj, k): return None
    return k, adj

# ================================================================
# DRIVER: range over constrained families, bucket iso classes
# ================================================================
METHODS = {
    'M1_snapshot_binder_order': method1,
    'M2_plateau_pin_order': method2,
    'M3_section_majority_vote': method3,
    'M4_pair_binding_cyclic': method4,
    'M5_halfgrid_parity_xor': method5,
    'M6_rotational_pinrank': method6,
}

def run_method(name, fn, families):
    """families: list of (label, list_of_speed_sets). Returns per-k realized classes."""
    realized = {}  # k -> set of canonkeys
    nontrivial_examples = {}  # k -> example
    total_inputs = 0
    toobig = 0
    for label, sets in families:
        for S in sets:
            total_inputs += 1
            res = fn(S)
            if res is None: continue
            if isinstance(res, tuple) and res[0] == 'TOOBIG':
                toobig += 1; continue
            k, adj = res
            if not is_tournament(adj, k): continue
            if k not in (3,4,5): continue
            key = tour_key(adj, k)
            realized.setdefault(k, set()).add(key)
            if num_3cycles(adj, k) > 0 and k not in nontrivial_examples:
                nontrivial_examples[k] = (S, k, num_3cycles(adj,k), num_hampaths(adj,k))
    return realized, nontrivial_examples, total_inputs, toobig

def report(name, realized, examples, total, toobig):
    print(f"\n{'='*70}\nMETHOD: {name}")
    print(f"  inputs scanned: {total}  (skipped-toobig: {toobig})")
    any_nontrivial = bool(examples)
    print(f"  NON-TRIVIAL (produces a 3-cycle anywhere): {any_nontrivial}")
    if examples:
        for k, ex in sorted(examples.items()):
            print(f"    e.g. k={k}: S={ex[0]} c3={ex[2]} H={ex[3]}")
    for k in (3,4,5):
        rs = realized.get(k, set())
        if not rs: continue
        free = set(FREE[k].keys())
        forbidden = free - rs
        print(f"  k={k}: realized {len(rs)}/{A000568[k]} free classes.")
        if forbidden:
            print(f"    FORBIDDEN ({len(forbidden)}):")
            for key in sorted(forbidden):
                print(f"       {class_label(k, key)}")
        else:
            print(f"    (realizes ALL free classes)")
    return any_nontrivial

# ----------------------------------------------------------------
# Build input families. We want LRC-CONSTRAINED inputs at small vertex
# counts so canonicalization (k<=5) is exhaustive.
#  - n=3 speeds (k may be 3): all primitive 3-subsets of 1..14
#  - n=4 speeds: all primitive 4-subsets of 1..12
#  - n=5 speeds: all primitive 5-subsets of 1..10 (limited for runtime)
# Methods that use runners-as-vertices: k = nspeed.
# Methods that use pairs/binders: k varies.
# ----------------------------------------------------------------
def main():
    fam3 = [('3of14', small_speed_sets(3, 14))]
    fam4 = [('4of12', small_speed_sets(4, 12))]
    fam5 = [('5of10', small_speed_sets(5, 10))]
    print("Family sizes:",
          "3of14=",len(fam3[0][1]),
          "4of12=",len(fam4[0][1]),
          "5of10=",len(fam5[0][1]))

    # methods with runners-as-vertices want each nspeed separately
    runner_families = fam3 + fam4 + fam5

    summary = {}
    for name, fn in METHODS.items():
        realized, examples, total, toobig = run_method(name, fn, runner_families)
        nt = report(name, realized, examples, total, toobig)
        summary[name] = (nt, realized)

    # Special deeper dive for any NON-TRIVIAL method with a forbidden class
    print("\n" + "#"*70)
    print("SUMMARY")
    for name, (nt, realized) in summary.items():
        flags = []
        for k in (3,4,5):
            rs = realized.get(k,set())
            if rs:
                fb = set(FREE[k].keys()) - rs
                flags.append(f"k{k}:{len(rs)}/{A000568[k]}" + (f" FORBIDS{len(fb)}" if fb else ""))
        print(f"  {name}: nontrivial={nt}  {' '.join(flags)}")

if __name__ == '__main__':
    main()
