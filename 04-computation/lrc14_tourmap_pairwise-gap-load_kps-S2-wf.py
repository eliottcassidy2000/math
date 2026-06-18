"""
lrc14_tourmap_pairwise-gap-load_kps-S2-wf.py

THEME: "pairwise-gap-load" tournament generators from Lonely Runner data.

We map LRC inputs (speed set S, lonely time tau, sections on the grid) to
tournaments via NON-OBVIOUS arc rules, then ask which tournament ISO CLASSES
are realized as (S, tau) range over LRC-constrained inputs. A class in the FREE
set (all A000568 iso classes) but NEVER realized = a FORBIDDEN CLASS.

All decisions use EXACT rationals (fractions.Fraction). No floats.

Methods (theme = pairwise-gap-load), each must pass a NON-TRIVIALITY gate
(does it ever produce a 3-cycle / H>1?) before we enumerate realized classes.

Validated M tool copied verbatim from the task prompt.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product
import sys

# ----------------------------------------------------------------------
# Validated M tool (verbatim)
# ----------------------------------------------------------------------
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

# Signed normalized residue (for direction-sensitive rules): in (-1/2, 1/2]
def snrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r  # r in [0,1)
    # map to (-1/2, 1/2]
    return r if r <= F(1, 2) else r - 1

# ----------------------------------------------------------------------
# Tournament canonicalization utilities
# ----------------------------------------------------------------------
def adj_from_arcfn(verts, arcfn):
    """arcfn(i,j) returns True if arc i->j (i beats j). Must be a tournament:
    for i!=j exactly one of arcfn(i,j), arcfn(j,i) true. Returns adjacency
    dict-of-frozenset of out-neighbors, plus a validity flag."""
    n = len(verts)
    A = [[0]*n for _ in range(n)]
    valid = True
    for a in range(n):
        for b in range(a+1, n):
            i, j = verts[a], verts[b]
            ij = arcfn(i, j)
            ji = arcfn(j, i)
            if ij == ji:
                valid = False  # not a tournament (tie or double)
            if ij:
                A[a][b] = 1
            else:
                A[b][a] = 1
    return A, valid

def score_seq(A):
    n = len(A)
    return tuple(sorted(sum(A[i]) for i in range(n)))

def count_3cycles(A):
    n = len(A); c = 0
    for i in range(n):
        for j in range(n):
            if i == j or not A[i][j]:
                continue
            for k in range(n):
                if k == i or k == j:
                    continue
                if A[j][k] and A[k][i]:
                    c += 1
    return c // 3  # each directed 3-cycle counted 3x

def canon(A):
    """Brute canonical form via all n! relabelings: lexicographically
    smallest flattened upper-triangular-ish full matrix tuple."""
    n = len(A)
    best = None
    for p in permutations(range(n)):
        flat = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or flat < best:
            best = flat
    return best

def H_hampaths(A):
    """Number of Hamiltonian paths (directed) = Redei H."""
    n = len(A)
    cnt = 0
    for p in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if not A[p[i]][p[i+1]]:
                ok = False; break
        if ok:
            cnt += 1
    return cnt

# Precompute the full free set of iso classes for n=3,4,5 (by exhaustive
# enumeration of all tournaments). Keyed by canon().
def free_set(n):
    classes = {}
    pairs = list(combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = [[0]*n for _ in range(n)]
        for idx, (a, b) in enumerate(pairs):
            if (bits >> idx) & 1:
                A[a][b] = 1
            else:
                A[b][a] = 1
        c = canon(A)
        if c not in classes:
            classes[c] = (score_seq(A), count_3cycles(A), H_hampaths(A))
    return classes

FREE = {n: free_set(n) for n in (3, 4, 5)}
print("FREE SET sizes (should be A000568 = 2,4,12):")
for n in (3, 4, 5):
    print(f"  n={n}: {len(FREE[n])} classes")
print()

# ----------------------------------------------------------------------
# LRC input generation: small primitive speed sets and their lonely times.
# We enumerate primitive speed sets of size n with speeds in a bounded range,
# compute M(S) and the optimal tau, and feed (S, tau) into each method.
# ----------------------------------------------------------------------
def is_primitive(S):
    g0 = 0
    for v in S:
        g0 = gcd(g0, v)
    return g0 == 1

def gen_speedsets(n, vmax):
    """All primitive n-subsets of {1..vmax} (positive speeds)."""
    out = []
    for combo in combinations(range(1, vmax+1), n):
        if is_primitive(combo):
            out.append(combo)
    return out

# ----------------------------------------------------------------------
# METHOD DEFINITIONS
# Each method: given (S, tau) produce an arcfn over a vertex set.
# We return (verts, arcfn) or None if undefined for this input.
# ----------------------------------------------------------------------

# --- METHOD 1: Gap-margin tournament on RUNNERS ---
# Vertex = runner i (speed v_i). At the lonely time tau, each runner has
# distance d_i = ||v_i tau||. Arc i->j by "who is more binding / closer to 0":
# i->j iff d_i < d_j (i is tighter/closer to the danger boundary). Ties broken
# by the SIGNED residue (smaller signed residue wins). This is NOT the frac
# total-order (which is the dead overtaking map): it folds frac to distance,
# then compares distances. Folding destroys the total order -> can be cyclic.
def method1(S, tau):
    S = sorted(set(S))
    verts = list(S)
    d = {v: nrm(v*tau) for v in S}
    sgn = {v: snrm(v*tau) for v in S}
    def arcfn(i, j):
        if d[i] != d[j]:
            return d[i] < d[j]
        # tie on distance: compare signed residue, then speed (total deterministic)
        if sgn[i] != sgn[j]:
            return sgn[i] < sgn[j]
        return i < j
    return verts, arcfn

# --- METHOD 2: Difference-speed danger tournament on RUNNERS ---
# Arc i->j by the DIFFERENCE speed danger: how close is ||(v_i - v_j) tau||
# combined with orientation. Concretely we orient i->j iff the signed residue
# of (v_i - v_j)*tau is positive (i pulls ahead of j at tau in the circle
# direction). This is the SIGN of the pairwise gap on the difference speed.
# Sign of snrm((v_i-v_j)tau): if positive i->j else j->i. snrm is antisymmetric
# (snrm(-x) = -snrm(x) except at the 1/2 boundary), so this is a valid
# tournament except when snrm == 0 (then degenerate -> we break by index).
def method2(S, tau):
    S = sorted(set(S))
    verts = list(S)
    def arcfn(i, j):
        s = snrm((i - j) * tau)
        if s > 0:
            return True
        if s < 0:
            return False
        return i < j  # boundary tie-break
    return verts, arcfn

# --- METHOD 3: Section-distance tournament on RUNNERS (off-grid analog) ---
# Use the on-grid sections at a NEARBY grid time. But to stay general (S of any
# size, tau off-grid), define a "load" L_i = numerator of the binding fraction:
# L_i = the closest half-integer count, i.e. round(2 v_i tau). Arc i->j iff
# (L_i mod something) ... Instead use a CLEAN 3-runner-free rule:
# i->j iff ||v_i tau|| + ||v_j tau|| is closer to v-weighted...
# We define: i->j iff frac((v_i+v_j) tau) has signed residue > 0 (SUM-speed
# orientation). Pairwise SUM analog of method 2. Sum speeds are exactly the
# OTHER binding-pair family (THM-524 uses both v_i+v_j and v_i-v_j).
def method3(S, tau):
    S = sorted(set(S))
    verts = list(S)
    def arcfn(i, j):
        s = snrm((i + j) * tau)
        if s > 0:
            return True
        if s < 0:
            return False
        return i < j
    return verts, arcfn

# --- METHOD 4: 3-runner MAJORITY / Condorcet tournament on RUNNERS ---
# Arc i->j = majority vote over all other runners k of the relation
# "i is closer than j to the danger boundary as seen relative to k", i.e.
# count k for which ||(v_i - v_k) tau|| < ||(v_j - v_k) tau||. i->j iff i wins
# more head-to-head difference-distance comparisons. This is a genuine
# Condorcet aggregation -> classically can be cyclic. Tie-break by index.
def method4(S, tau):
    S = sorted(set(S))
    verts = list(S)
    dd = {}  # dd[(i,k)] = ||(i-k) tau||
    for a in S:
        for b in S:
            dd[(a, b)] = nrm((a - b) * tau)
    def arcfn(i, j):
        wins_i = 0; wins_j = 0
        for k in S:
            if k == i or k == j:
                continue
            di = dd[(i, k)]; dj = dd[(j, k)]
            if di < dj:
                wins_i += 1
            elif dj < di:
                wins_j += 1
        if wins_i != wins_j:
            return wins_i > wins_j
        return i < j
    return verts, arcfn

# --- METHOD 5: Pairwise-load tournament on PAIRS (vertex = unordered pair) ---
# For a base speed set, vertices are the C(n,2) pairs {i,j}. The "load" of a
# pair is its binding distance ||(v_i+v_j) tau|| ... we want a small vertex
# count so use n=3 base (3 pairs -> tournament on 3 vertices), n=4 base
# (6 pairs -> too big for n5 free set), so restrict pair-tournaments to base
# sizes giving <=5 pairs: base size 3 (3 pairs) only for clean iso classes.
# Arc between pair P=(a,b) and pair Q=(c,d): compare the SUM-speed binding
# distances: P->Q iff ||(a+b)tau|| < ||(c+d)tau|| (P is more binding). Tie ->
# lexicographic. Folding via nrm again -> potentially cyclic.
def method5(S, tau):
    S = sorted(set(S))
    prs = list(combinations(S, 2))
    verts = list(range(len(prs)))
    load = [nrm((a+b)*tau) for (a, b) in prs]
    def arcfn(i, j):
        if load[i] != load[j]:
            return load[i] < load[j]
        return prs[i] < prs[j]
    return verts, arcfn

# --- METHOD 6: Parity/switch tournament on RUNNERS ---
# Built from a SWITCH count, not a snapshot. For the optimal tau (a rational
# p/q), the number of times runner i has "crossed the danger boundary" up to
# tau is floor(2 v_i tau + 1/2) = round-ish count = number of half-integers
# <= v_i tau. Define for each runner the integer
#   c_i = number of integers m>=1 with m/(2) ...
# Concretely c_i = floor(v_i * tau * 2) (count of half-laps completed).
# Arc i->j iff (c_i + c_j) is ODD  ... that's symmetric, not a tournament.
# Instead: arc i->j iff c_i > c_j (more half-laps); tie -> by sign residue;
# tie -> index. c_i is a near-linear monotone in v_i for fixed tau, so this is
# likely transitive -> we test the gate honestly.
def method6(S, tau):
    S = sorted(set(S))
    verts = list(S)
    c = {v: int(2 * v * tau) for v in S}  # floor of 2 v tau
    sgn = {v: snrm(v*tau) for v in S}
    def arcfn(i, j):
        if c[i] != c[j]:
            return c[i] > c[j]
        if sgn[i] != sgn[j]:
            return sgn[i] < sgn[j]
        return i < j
    return verts, arcfn

METHODS = {
    "M1_gap_margin_runner": method1,
    "M2_diff_speed_sign": method2,
    "M3_sum_speed_sign": method3,
    "M4_condorcet_diffdist": method4,
    "M5_pair_sumload": method5,
    "M6_halflap_count": method6,
}

# ----------------------------------------------------------------------
# Driver: for each method and each base size n_base, range over speed sets,
# build the tournament at the OPTIMAL lonely tau, canonicalize, bucket.
# ----------------------------------------------------------------------
def run_method(name, fn, n_base, vmax, also_all_cand=False, verbose=False):
    """n_base = number of runners (so vertex count depends on method:
    runner-methods -> n_base vertices; pair-method -> C(n_base,2) vertices)."""
    speedsets = gen_speedsets(n_base, vmax)
    realized = {}      # canon -> (score_seq, c3, H, example)
    nonvalid = 0
    total = 0
    # determine vertex count for this method by probing one example
    if not speedsets:
        return None
    for S in speedsets:
        Mval, tau = M(S)
        taus = [tau]
        if also_all_cand:
            taus = list(cand(S))
        for t in taus:
            res = fn(S, t)
            if res is None:
                continue
            verts, arcfn = res
            A, valid = adj_from_arcfn(verts, arcfn)
            total += 1
            if not valid:
                nonvalid += 1
                continue
            c = canon(A)
            if c not in realized:
                realized[c] = (score_seq(A), count_3cycles(A), H_hampaths(A), tuple(S), t)
    return {
        "n_vertices": len(verts),
        "realized": realized,
        "nonvalid": nonvalid,
        "total": total,
        "n_speedsets": len(speedsets),
    }

# Map n_base -> vertex count per method to compare against the right FREE set
def vertex_count_for(name, n_base):
    if name == "M5_pair_sumload":
        return n_base*(n_base-1)//2
    return n_base

# ----------------------------------------------------------------------
# RUN: non-triviality gate + class enumeration
# ----------------------------------------------------------------------
print("="*70)
print("NON-TRIVIALITY GATE + ISO-CLASS ENUMERATION (at optimal lonely tau)")
print("="*70)

# Choose base sizes so the resulting vertex count is in {3,4,5} (free set known).
# Runner methods: n_base in {3,4,5}. Pair method: n_base=3 -> 3 verts.
configs = {
    "M1_gap_margin_runner": [(3, 12), (4, 12), (5, 11)],
    "M2_diff_speed_sign":   [(3, 14), (4, 12), (5, 11)],
    "M3_sum_speed_sign":    [(3, 14), (4, 12), (5, 11)],
    "M4_condorcet_diffdist":[(3, 14), (4, 12), (5, 11)],
    "M5_pair_sumload":      [(3, 16)],   # 3 pairs -> 3 verts; base 4 -> 6 verts (skip)
    "M6_halflap_count":     [(3, 14), (4, 12), (5, 11)],
}

summary = {}
for name, fn in METHODS.items():
    print(f"\n##### {name} #####")
    summary[name] = {}
    for (n_base, vmax) in configs[name]:
        vc = vertex_count_for(name, n_base)
        if vc not in FREE:
            print(f"  base n={n_base}: vertex count {vc} not in known free set, skip")
            continue
        r = run_method(name, fn, n_base, vmax)
        if r is None:
            print(f"  base n={n_base}: no primitive speed sets")
            continue
        realized = r["realized"]
        free = FREE[vc]
        # non-trivial = any realized class with c3>0 (H>1)
        nontrivial = any(v[1] > 0 for v in realized.values())
        forbidden = [(c, free[c]) for c in free if c not in realized]
        print(f"  base n={n_base} (vmax={vmax}): vertices={vc}, "
              f"speedsets={r['n_speedsets']}, valid-tournaments={r['total']-r['nonvalid']}, "
              f"non-valid(ties)={r['nonvalid']}")
        print(f"    realized {len(realized)}/{len(free)} free classes; "
              f"NON-TRIVIAL={nontrivial}")
        # report realized class signatures (score_seq, c3, H)
        sigs = sorted(set((v[0], v[1], v[2]) for v in realized.values()))
        print(f"    realized signatures (score,c3,H): {sigs}")
        if forbidden:
            fsigs = sorted(set((v[0], v[1], v[2]) for (_, v) in forbidden))
            print(f"    FORBIDDEN {len(forbidden)} classes; signatures: {fsigs}")
        summary[name][n_base] = {
            "vertices": vc,
            "realized_count": len(realized),
            "free_count": len(free),
            "nontrivial": nontrivial,
            "forbidden_count": len(forbidden),
            "forbidden_sigs": sorted(set((v[0], v[1], v[2]) for (_, v) in forbidden)) if forbidden else [],
            "realized_sigs": sigs,
        }

# ----------------------------------------------------------------------
# DEEPER PASS: for the most promising non-trivial methods, ALSO range tau over
# the WHOLE candidate set (not just the optimal lonely tau) to see if loneliness
# (restricting to the optimal tau) genuinely forbids classes that arbitrary
# crossing-times would reach. This isolates the loneliness constraint.
# ----------------------------------------------------------------------
print("\n" + "="*70)
print("LONELINESS-vs-ALL-CROSSINGS COMPARISON (does optimal-tau forbid more?)")
print("="*70)
for name in ("M1_gap_margin_runner", "M4_condorcet_diffdist", "M2_diff_speed_sign", "M3_sum_speed_sign"):
    fn = METHODS[name]
    for (n_base, vmax) in [(4, 10), (5, 9)]:
        vc = vertex_count_for(name, n_base)
        if vc not in FREE:
            continue
        r_opt = run_method(name, fn, n_base, vmax, also_all_cand=False)
        r_all = run_method(name, fn, n_base, vmax, also_all_cand=True)
        free = FREE[vc]
        ro = set(r_opt["realized"].keys())
        ra = set(r_all["realized"].keys())
        print(f"\n  {name} base n={n_base} vmax={vmax} (verts={vc}):")
        print(f"    optimal-tau realized: {len(ro)}/{len(free)}")
        print(f"    all-crossings realized: {len(ra)}/{len(free)}")
        only_all = ra - ro
        never = set(free.keys()) - ra
        print(f"    classes reached by SOME crossing but NOT at optimal lonely tau: {len(only_all)}")
        if only_all:
            os = sorted(set((free[c][0], free[c][1], free[c][2]) for c in only_all))
            print(f"      -> their signatures (score,c3,H): {os}")
        print(f"    classes NEVER reached even over all crossings: {len(never)}")
        if never:
            ns = sorted(set((free[c][0], free[c][1], free[c][2]) for c in never))
            print(f"      -> signatures: {ns}")

print("\nDONE.")
