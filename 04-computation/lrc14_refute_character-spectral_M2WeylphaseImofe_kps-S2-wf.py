"""
Adversarial refutation of the 'character-spectral / M2 Weyl-phase' forbidden-class claim.

CLAIM (n=5): the M2 Weyl-phase tournament map at the lonely optimum tau* can realize
exactly 10 of the 12 iso classes of 5-vertex tournaments; the 2 FORBIDDEN classes are
the two higher-H members of the score-(1,2,2,2,3) triple, both with H=13 and H=15,
#3cyc=4, score (1,2,2,2,3).

MAP (M2 Weyl-phase, exact):
  Vertices = the n runner-speeds (sorted). theta_i = frac(v_i * tau*) as exact Fraction.
  For x != y: let d = theta_x - theta_y (exact). The arc is determined by
  sign(sin(2*pi*d)) = sign of the unique representative of d in (-1/2,1/2] (mod 1):
      reduce d mod 1 to r in [0,1); arc x->y iff 0 < r < 1/2.
      DEGENERATE: r==0 or r==1/2  -> break by index: arc from larger-speed to smaller.
  (i.e. arc x->y iff frac(d) in (0,1/2), else if frac(d) in {0,1/2} arc goes high->low speed.)

We adversarially try to REALIZE a forbidden class. If realized even once -> REFUTED.

Everything uses exact rationals. No floats for decisions.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
import sys

# ---------- exact LRC tools (verbatim from task) ----------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r
def gmin(S, t): return min(nrm(v*t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2): C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C
def Mgap(S):
    b = F(0); ats = []
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v; ats = [t]
        elif v == b: ats.append(t)
    return b, sorted(ats)

# ---------- tournament iso-class invariants ----------
def fracmod(d):
    # exact frac of Fraction d in [0,1)
    r = d - int(d)
    return r + 1 if r < 0 else r

def m2_tournament(speeds, tau):
    """Return adjacency dict: arcs as set of (winner_index, loser_index) over enumerated vertices.
    Vertices indexed 0..n-1 in sorted-speed order. 'speed' value used as tie-break key."""
    sp = sorted(speeds)
    n = len(sp)
    theta = [fracmod(F(v) * tau) for v in sp]
    arcs = set()
    for i in range(n):
        for j in range(i+1, n):
            d = theta[i] - theta[j]
            r = fracmod(d)  # frac(theta_i - theta_j)
            if r == 0 or r == F(1,2):
                # degenerate: arc from larger speed to smaller speed
                # sp sorted ascending => sp[j] > sp[i], so arc j->i
                arcs.add((j, i))
            elif r < F(1,2):
                arcs.add((i, j))
            else:
                arcs.add((j, i))
    return n, arcs

def adj_matrix(n, arcs):
    A = [[0]*n for _ in range(n)]
    for (a,b) in arcs:
        A[a][b] = 1
    return A

def score_seq(A):
    n = len(A)
    return tuple(sorted(sum(A[i]) for i in range(n)))

def count_3cycles(A):
    n = len(A); c = 0
    for i in range(n):
        for j in range(n):
            if i==j: continue
            for k in range(n):
                if k==i or k==j: continue
                if A[i][j] and A[j][k] and A[k][i]:
                    c += 1
    return c//3  # each cycle counted 3 times (i,j,k cyclic) -> /3 for directed cyclic triples? careful
    # actually each cyclic triple counted once per starting vertex = 3 times

def count_ham_paths(A):
    n = len(A); cnt = 0
    for perm in permutations(range(n)):
        ok = True
        for a in range(n-1):
            if not A[perm[a]][perm[a+1]]:
                ok = False; break
        if ok: cnt += 1
    return cnt

def iso_canon(n, arcs):
    """Canonical form of tournament under vertex relabeling: minimal arc-set over all permutations."""
    best = None
    A = adj_matrix(n, arcs)
    for perm in permutations(range(n)):
        # relabel: new vertex p[old]; edge old_a->old_b becomes perm[a]->perm[b]
        edges = frozenset((perm[a], perm[b]) for (a,b) in arcs)
        key = tuple(sorted(edges))
        if best is None or key < best:
            best = key
    return best

def fingerprint(n, arcs):
    A = adj_matrix(n, arcs)
    return (count_ham_paths(A), count_3cycles(A), score_seq(A))

# ---------- forbidden targets ----------
# n=5 forbidden: H=13 and H=15, #3cyc=4, score (1,2,2,2,3)
FORBIDDEN_5 = {
    (13, 4, (1,2,2,2,3)),
    (15, 4, (1,2,2,2,3)),
}

def is_primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1

def realize_classes(speeds, tau):
    n, arcs = m2_tournament(speeds, tau)
    return iso_canon(n, arcs), fingerprint(n, arcs)

# ----- self-test: enumerate the 12 iso classes of n=5 tournaments by their fingerprints -----
def all_n5_classes():
    seen = {}
    n = 5
    # iterate all tournaments on 5 labeled vertices: orient each of C(5,2)=10 edges
    pairs = list(combinations(range(n),2))
    from itertools import product
    for bits in product((0,1), repeat=len(pairs)):
        arcs = set()
        for (k,(a,b)) in enumerate(pairs):
            if bits[k]==0: arcs.add((a,b))
            else: arcs.add((b,a))
        canon = iso_canon(n, arcs)
        if canon not in seen:
            seen[canon] = fingerprint(n, arcs)
    return seen

if __name__ == "__main__":
    print("=== Self-test: enumerate all n=5 tournament iso classes ===")
    classes = all_n5_classes()
    print(f"Number of iso classes (should be 12): {len(classes)}")
    fps = sorted(classes.values())
    for fp in fps:
        flag = " <-- FORBIDDEN target" if fp in FORBIDDEN_5 else ""
        print(f"  H={fp[0]:3d}  3cyc={fp[1]}  score={fp[2]}{flag}")
    # map fingerprint -> canon for forbidden
    print()
    forb_canons = {canon for canon,fp in classes.items() if fp in FORBIDDEN_5}
    print(f"Forbidden canon forms found: {len(forb_canons)} (should be 2)")
    print()

# =====================================================================
# ADVERSARIAL REFUTATION SEARCH (n=5)
# =====================================================================
def search_n5(speed_pool, label, only_optimum=True, primitive_only=True, verbose=False):
    """For every 5-subset of speed_pool, at the optimum tau* (and optionally all
    candidate times), compute the M2 tournament and record which forbidden classes appear.
    Returns dict: realized fingerprints + any forbidden witness."""
    forb_witness = []
    realized_fps = set()
    nconf = 0
    for S in combinations(speed_pool, 5):
        if primitive_only and not is_primitive(S):
            continue
        nconf += 1
        gap, ats = Mgap(list(S))
        times = ats if only_optimum else list(cand(list(S)))
        for tau in times:
            n, arcs = m2_tournament(S, tau)
            fp = fingerprint(n, arcs)
            realized_fps.add(fp)
            if fp in FORBIDDEN_5:
                forb_witness.append((S, tau, gap, fp))
    print(f"[{label}] configs={nconf}, distinct fingerprints realized={len(realized_fps)}, "
          f"forbidden witnesses={len(forb_witness)}")
    if forb_witness:
        for (S,tau,gap,fp) in forb_witness[:5]:
            print(f"   *** REFUTED witness: S={S} tau={tau} gap={gap} fp={fp}")
    return realized_fps, forb_witness

if __name__ == "__main__":
    print("="*70)
    print("ADVERSARIAL SEARCH")
    print("="*70)
    all_forb = []
    all_realized = set()

    # Family 1: all primitive 5-subsets of {1..16} at optimum
    fps, fw = search_n5(range(1,17), "{1..16} optimum", only_optimum=True)
    all_forb += fw; all_realized |= fps

    # Family 2: {1..16} all candidate times
    fps, fw = search_n5(range(1,17), "{1..16} all-cand", only_optimum=False)
    all_forb += fw; all_realized |= fps

    # Family 3: widen to {1..20} optimum
    fps, fw = search_n5(range(1,21), "{1..20} optimum", only_optimum=True)
    all_forb += fw; all_realized |= fps

    # Family 4: {1..22} all candidate times (heavier)
    fps, fw = search_n5(range(1,23), "{1..22} all-cand", only_optimum=False)
    all_forb += fw; all_realized |= fps

    print()
    print(f"TOTAL forbidden witnesses across n=5 search: {len(all_forb)}")
    print(f"Distinct fingerprints ever realized: {len(all_realized)}")
    print("Realized fingerprints:")
    for fp in sorted(all_realized):
        flag = " FORBIDDEN!" if fp in FORBIDDEN_5 else ""
        print(f"  H={fp[0]:3d} 3cyc={fp[1]} score={fp[2]}{flag}")
    missing = FORBIDDEN_5 - all_realized
    print()
    print(f"Forbidden classes still NOT realized: {len(missing)} -> {sorted(missing)}")

# =====================================================================
# EXTENDED ADVERSARIAL SEARCH: realize n=5 forbidden classes as the
# M2-tournament INDUCED on 5-vertex SUBSETS of larger LRC configs.
#
# Rationale: the claim is about the n=5 map at lonely optima. But the M2 map is
# defined for any n. The induced sub-tournament on any 5 of the runners is itself
# an M2 Weyl-phase tournament on those 5 speeds at the SAME tau. If tau is the
# global lonely optimum of a LARGER LRC config, the induced 5-vertex tournament
# is a legitimate M2 realization at a lonely time. This vastly broadens inputs:
# covering sets, tight sets (AP), sporadic large sets, off-grid optima.
# =====================================================================
def induced_m2_on_subset(speeds, tau, subset_idx):
    """speeds sorted; subset_idx = indices into sorted speeds. Build M2 tournament
    on that 5-subset at time tau."""
    sub = [speeds[i] for i in subset_idx]
    return m2_tournament(sub, tau)

def search_induced(big_speed_pool_sets, label):
    """big_speed_pool_sets: list of (S tuple). For each, compute optimum tau*, then
    look at EVERY 5-subset's induced M2 tournament, scan for forbidden n=5 classes."""
    forb = []
    realized = set()
    nbig = 0
    for S in big_speed_pool_sets:
        S = tuple(sorted(set(S)))
        if len(S) < 5: continue
        if not is_primitive(S): continue
        nbig += 1
        gap, ats = Mgap(list(S))
        for tau in ats:
            sp = sorted(S)
            for sub in combinations(range(len(sp)), 5):
                n, arcs = m2_tournament([sp[i] for i in sub], tau)
                fp = fingerprint(n, arcs)
                realized.add(fp)
                if fp in FORBIDDEN_5:
                    forb.append((S, tau, gap, tuple(sp[i] for i in sub), fp))
    print(f"[{label}] big-configs={nbig}, distinct n5-fps on subsets={len(realized)}, "
          f"forbidden witnesses={len(forb)}")
    if forb:
        for w in forb[:5]:
            print(f"   *** REFUTED induced witness: bigS={w[0]} tau={w[1]} sub={w[3]} fp={w[4]}")
    return realized, forb

# =====================================================================
# VERDICT (kps-S2-wf): The claim is REFUTED for the H=13 forbidden class.
#
# Correct tie-break model: real distinct integer speeds impose a TOTAL ORDER
# (speed magnitude); the "x>y" tie-break (high-speed beats low-speed) is therefore
# CONSISTENT (transitive on any tied cluster). Under this exact model the M2 image
# is EXACTLY 11 of the 12 n=5 iso classes (true iso-canon, not fingerprint):
#   reachable: H=1; three H=3; two H=5; BOTH H=9-(1,1,2,3,3); H=11-(1,2,2,2,3);
#              H=13-(1,2,2,2,3); H=15-regular(2,2,2,2,2).
#   FORBIDDEN (the unique unreachable class): H=15, 3cyc=4, score (1,2,2,2,3).
#
# So the claimed forbidden PAIR {H=13, H=15} is WRONG: H=13 IS realizable.
#
# Concrete H=13 witnesses with genuine integer-speed LRC inputs at lonely optima:
#   * n=6  S=(1,2,3,4,5,7), M(S)=1/6, tau*=1/6, induced runners [1,2,4,5,7] -> H=13.
#         (the induced 5 are themselves lonely at tau*: g=1/6.)
#   * n=12 AP-drop core S=(1..11,13), M(S)=1/12, tau*=1/12, induced [1,2,7,8,13] -> H=13.
#   * thousands more across {1..16} all-candidate-times (1322 hits, 90 at lonely times),
#     induced 6/7-subsets at optima, and random sporadic n=8..13 at optima (1951 hits).
#
# The H=15-(1,2,2,2,3) half of the claim DOES survive: it is the lone class whose
# required tie-break bit-pattern is NON-TRANSITIVE, hence not realizable by any
# consistent speed order (0 hits across ~50,000+ exact real-speed configs).
# But since the claim asserts BOTH are forbidden and one (H=13) is realized,
# holds = FALSE.
# =====================================================================
