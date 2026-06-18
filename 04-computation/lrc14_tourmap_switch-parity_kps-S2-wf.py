#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_tourmap_switch-parity_kps-S2-wf.py

THEME: build a tournament from the LRC BINDING-PAIR SWITCHES / CROSSINGS /
PARITIES (DYNAMICS), not a single-time snapshot order. The snapshot
"overtaking" order (i->j iff frac(v_i tau)>frac(v_j tau)) is ALWAYS transitive
(frac is a total order) so it is DEAD (H=1). We avoid that.

We define >=4 distinct arc-rules from switch/crossing/parity data, gate each on
NON-TRIVIALITY (does it ever produce a 3-cycle / H>1?), and for the non-trivial
ones we enumerate which tournament ISO CLASSES are realized as (S,tau) range over
LRC-constrained inputs, comparing against the FULL free set (A000568).

A000568 (iso classes of tournaments): n=1..7 -> 1,1,2,4,12,56,456.

All exact rationals via fractions.Fraction. stdlib only.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
from collections import Counter, defaultdict

# ----------------------------------------------------------------------
# exact M tool (validated, verbatim from prompt)
# ----------------------------------------------------------------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2): C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v; at = t
    return b, at
def all_optima(S):
    Mv, _ = M(S)
    return Mv, sorted(t for t in cand(S) if gmin(S, t) == Mv)

# ----------------------------------------------------------------------
# tournament tools
# adj is a dict (i,j)->bool over an ordered vertex list V; exactly one of
# (i,j),(j,i) is True for i!=j.
# ----------------------------------------------------------------------
def ham_path_count(adj, V):
    n = len(V)
    dp = [[0] * n for _ in range(1 << n)]
    for k in range(n): dp[1 << k][k] = 1
    full = (1 << n) - 1
    for mask in range(1 << n):
        for last in range(n):
            c = dp[mask][last]
            if c == 0: continue
            for nxt in range(n):
                if mask & (1 << nxt): continue
                if adj[(V[last], V[nxt])]:
                    dp[mask | (1 << nxt)][nxt] += c
    return sum(dp[full][k] for k in range(n))

def scores(adj, V):
    return tuple(sorted(sum(1 for u in V if u != v and adj[(v, u)]) for v in V))

def count_3cycles(adj, V):
    c = 0
    for a, b, cc in combinations(V, 3):
        # 3-cycle iff not transitive on the triple; count cyclic triangles
        s = sum(1 for (x, y) in [(a, b), (b, cc), (cc, a)] if adj[(x, y)])
        # cyclic if all-forward or all-backward in this fixed cyclic order
        if s == 3 or s == 0:
            c += 1
    return c

def matrix(adj, V):
    n = len(V)
    return tuple(tuple(1 if (i != j and adj[(V[i], V[j])]) else 0
                       for j in range(n)) for i in range(n))

def canon(adj, V):
    """Canonical key: minimal adjacency matrix over all relabelings (n<=5 ok)."""
    n = len(V)
    base = [[1 if (i != j and adj[(V[i], V[j])]) else 0 for j in range(n)]
            for i in range(n)]
    best = None
    for p in permutations(range(n)):
        m = tuple(tuple(base[p[i]][p[j]] for j in range(n)) for i in range(n))
        if best is None or m < best:
            best = m
    return best

def is_tournament(adj, V):
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            a, b = V[i], V[j]
            if adj[(a, b)] == adj[(b, a)]:
                return False
    return True

# Full free set (iso classes) by exhaustive build over all 2^C(n,2) tournaments.
def free_classes(n):
    V = list(range(n))
    pairs = list(combinations(range(n), 2))
    seen = {}
    for bits in range(1 << len(pairs)):
        adj = {}
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1:
                adj[(i, j)] = True; adj[(j, i)] = False
            else:
                adj[(i, j)] = False; adj[(j, i)] = True
        k = canon(adj, V)
        if k not in seen:
            seen[k] = {"score": scores(adj, V),
                       "c3": count_3cycles(adj, V),
                       "H": ham_path_count(adj, V)}
    return seen

# ----------------------------------------------------------------------
# LRC input families. For small n (3..5 runners) we want to range over many
# primitive speed sets S with a positive gap, and use the optimal tau.
# ----------------------------------------------------------------------
def primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1

def speed_sets(n, maxspeed):
    """All n-subsets of {1..maxspeed}, primitive, distinct, with positive gap."""
    out = []
    for S in combinations(range(1, maxspeed + 1), n):
        if not primitive(S):
            continue
        out.append(S)
    return out

# crossing data helpers ------------------------------------------------
def frac(x):
    r = x - int(x)
    return r + 1 if r < 0 else r

def signed_dev(v, t):
    """signed deviation of runner v from observer at 0 on circle, in (-1/2,1/2].
       i.e. the representative of v*t mod 1 in (-1/2,1/2]."""
    r = frac(v * t)
    return r - 1 if r > F(1, 2) else r

def triangle_slope(v, t):
    """sign of d/dt ||v t|| at t (the triangle wave's local slope), +/-1.
       ||v t|| rises when frac(v t) in (0,1/2), falls when in (1/2,1).
       At exact peaks/troughs (frac=0 or 1/2) slope is ill-defined; return 0."""
    r = frac(v * t)
    if r == 0 or r == F(1, 2):
        return 0
    return v if r < F(1, 2) else -v  # scaled by speed: rate is v * (+/-1)

def crossing_count_diff(vi, vj, t):
    """# of times runner i has lapped/crossed runner j on (0,t], by difference
       speed: floor((vi-vj)*t) style. Use signed count = floor((vi-vj)*t)."""
    d = (vi - vj) * t
    return d.__floor__()

def crossing_count_sum(vi, vj, t):
    """# of head-on meetings (sum-speed crossings) on (0,t]: floor((vi+vj)*t)."""
    d = (vi + vj) * t
    return d.__floor__()

# ======================================================================
# METHOD 1: DIFFERENCE-CROSSING PARITY tournament.
#   Vertices = runners. Arc i->j iff parity of N_ij := floor((v_i - v_j)*tau*)
#   ... but that is antisymmetric-ish; we need a real orientation rule.
#   RULE: for i<j (as speeds), let c = floor((v_i+v_j)*tau*) (# head-on meetings)
#   XOR floor(|v_i-v_j|*tau*) (# overtakings). Orient i->j if c is EVEN, else j->i.
#   This blends BOTH crossing families and is NOT a snapshot order.
# ======================================================================
def method1_adj(S, t):
    V = list(range(len(S)))
    sp = list(S)
    adj = {}
    for i, j in combinations(V, 2):
        vi, vj = sp[i], sp[j]
        c_sum = crossing_count_sum(vi, vj, t)
        c_dif = crossing_count_diff(max(vi, vj), min(vi, vj), t)
        par = (c_sum ^ c_dif) & 1
        if par == 0:
            adj[(i, j)] = True; adj[(j, i)] = False
        else:
            adj[(i, j)] = False; adj[(j, i)] = True
    return V, adj

# ======================================================================
# METHOD 2: BINDING-PAIR SLOPE tournament (rising triangle-wave).
#   Vertices = runners. Arc i->j iff at the optimal tau*, runner i's distance
#   ||v_i tau|| is RISING and j's is FALLING, OR (both same slope sign) the one
#   with larger |slope*sign| ... we need a strict tournament. RULE:
#   define key(v) = ( sign of slope of ||v t|| , signed_dev(v,t) ). Orient i->j
#   iff key(i) > key(j) lexicographically. The PRIMARY sort is by RISING vs
#   FALLING (a switch/dynamics quantity), tie-broken by position.
#   NOTE: lexicographic total order => transitive. So this is a GATE-FAIL test
#   by design (we include it to confirm the dead pattern and contrast Method 3).
# ======================================================================
def method2_adj(S, t):
    V = list(range(len(S)))
    sp = list(S)
    key = {}
    for i in V:
        sl = triangle_slope(sp[i], t)
        sgn = 1 if sl > 0 else (-1 if sl < 0 else 0)
        key[i] = (sgn, signed_dev(sp[i], t), sp[i])
    adj = {}
    for i, j in combinations(V, 2):
        if key[i] > key[j]:
            adj[(i, j)] = True; adj[(j, i)] = False
        else:
            adj[(i, j)] = False; adj[(j, i)] = True
    return V, adj

# ======================================================================
# METHOD 3: PAIRWISE-MEETING-PARITY tournament (the real candidate).
#   Vertices = runners. For each unordered pair we look at the LAST crossing
#   (binding-pair switch) of i and j strictly before/at tau*, and orient by the
#   PARITY OF THE CROSSING INDEX k where their relative phase wrapped.
#   Concretely: the pair (i,j) meet (||v_i t||=||v_j t|| via sum OR difference)
#   at times k/(v_i+v_j) and k/(v_i-v_j). Let k* = the largest k with
#   k/(v_i+v_j) <= tau* among SUM crossings (the # of head-on meetings so far),
#   and orient i->j iff (k* is even) == (v_i < v_j). This couples the crossing
#   PARITY with the speed order, producing a non-snapshot, potentially cyclic map.
# ======================================================================
def method3_adj(S, t):
    V = list(range(len(S)))
    sp = list(S)
    adj = {}
    for i, j in combinations(V, 2):
        vi, vj = sp[i], sp[j]
        # number of head-on (sum) meetings strictly within (0, t]
        kstar = ((vi + vj) * t).__floor__()
        cond = ((kstar & 1) == 0) == (vi < vj)
        if cond:
            adj[(i, j)] = True; adj[(j, i)] = False
        else:
            adj[(i, j)] = False; adj[(j, i)] = True
    return V, adj

# ======================================================================
# METHOD 4: SECTION-WINDING tournament (off-grid generalization of sections).
#   Vertices = runners. At tau*, runner i sits at signed position p_i in
#   (-1/2,1/2]. Define the WINDING count w_i = floor(v_i * tau*) (# full laps).
#   Orient i->j iff (w_i + [p_i>0]) and (w_j + [p_j>0]) ... we use a genuinely
#   2D rule: orient i->j iff ( (w_i mod 2, sign(p_i)) beats (w_j mod 2, sign(p_j)) )
#   under the cyclic (rock-paper-scissors-like) order on the 4 states
#   (0,+)->(1,+)->(1,-)->(0,-)->(0,+) ... with ties (same state) broken by |p|.
#   A CYCLIC order on states can produce 3-cycles among runners in different states.
# ======================================================================
# 4 states encoded 0..3: state(v) = 2*(w mod 2) + (1 if p>0 else 0)
#   mapping: (w even,p<=0)=0 ; (w even,p>0)=1 ; (w odd,p<=0)=2 ; (w odd,p>0)=3
# cyclic dominance: define a tournament on the 4 states (a "rps4" oriented 4-cycle
#   plus diagonals): state x beats state y iff (x - y) mod 4 in {1,2}? That gives
#   a regular tournament on 4 states is impossible (4 even). Use a fixed
#   NON-TRANSITIVE base tournament on the 4 states (a near-regular one with a
#   3-cycle) and lift it to runners; same-state ties broken by |p| (then by speed).
# ======================================================================
# Base tournament on 4 states: choose the one with a 3-cycle 1->2->3->1 and
#   0 dominated/dominating to taste. We'll use the "almost transitive with one
#   back-edge" tournament: standard order 0<1<2<3 EXCEPT 3->1 (back edge),
#   creating 3-cycle 1->2->3->1. Wait need 1->2,2->3,3->1 => yes cycle.
STATE_BEATS = {}
def _build_state_tour():
    # default transitive a->b iff a<b
    for a in range(4):
        for b in range(4):
            if a == b: continue
            STATE_BEATS[(a, b)] = (a < b)
    # introduce back edges to make a 3-cycle on {1,2,3}: keep 1->2,2->3, flip 3 vs1
    STATE_BEATS[(3, 1)] = True
    STATE_BEATS[(1, 3)] = False
_build_state_tour()

def _state(v, t):
    w = (v * t).__floor__()
    p = signed_dev(v, t)
    return 2 * (w & 1) + (1 if p > 0 else 0), p

def method4_adj(S, t):
    V = list(range(len(S)))
    sp = list(S)
    st = {i: _state(sp[i], t) for i in V}
    adj = {}
    for i, j in combinations(V, 2):
        si, pi = st[i]; sj, pj = st[j]
        if si != sj:
            iwins = STATE_BEATS[(si, sj)]
        else:
            # same state: break by |p| then by speed (total order within state)
            iwins = (abs(pi), sp[i]) > (abs(pj), sp[j])
        if iwins:
            adj[(i, j)] = True; adj[(j, i)] = False
        else:
            adj[(i, j)] = False; adj[(j, i)] = True
    return V, adj

# ======================================================================
# METHOD 5: PAIR-SECTION-COLLISION tournament on PAIRS (vertices = pairs of runners).
#   Different vertex set! Vertices = the C(n,2) runner-pairs. For two pairs P,Q
#   sharing... too dense. Instead: vertices = runners, arc i->j iff the number of
#   OTHER runners k strictly between i and j on the circle at tau* is EVEN.
#   "Between" = in the open arc from p_i to p_j going counterclockwise. This is a
#   PARITY-OF-SEPARATION rule (a genuine combinatorial switch count), not a
#   snapshot order. Antisymmetry: #between(i,j ccw) + #between(j,i ccw) = n-2
#   (constant), so parities are equal iff n-2 even => need care; define
#   orient i->j iff #strictly-between-ccw(i->j) is EVEN; if both even (n even...)
#   tie-break: this can fail antisymmetry, so we add the cyclic-position sign.
# ======================================================================
def method5_adj(S, t):
    V = list(range(len(S)))
    sp = list(S)
    pos = {i: frac(sp[i] * t) for i in V}  # in [0,1)
    adj = {}
    order = sorted(V, key=lambda i: (pos[i], sp[i]))
    rank = {v: r for r, v in enumerate(order)}
    n = len(V)
    for i, j in combinations(V, 2):
        # number of runners strictly between i and j going ccw from i to j
        ri, rj = rank[i], rank[j]
        # ccw from i to j: ranks strictly between ri and rj in cyclic order
        if ri < rj:
            between = rj - ri - 1
        else:
            between = (n - ri - 1) + rj
        par = between & 1
        # orient i->j iff par==0, but ensure antisymmetry: between(i->j ccw) +
        #  between(j->i ccw) = n-2; if n-2 even both have same parity -> tie.
        # break ties by rank order.
        if (n - 2) % 2 == 0:
            # parity identical both ways: degenerate. Use rank+parity blend:
            iwins = (par == 0) ^ (rank[i] > rank[j])
        else:
            iwins = (par == 0)
        if iwins:
            adj[(i, j)] = True; adj[(j, i)] = False
        else:
            adj[(i, j)] = False; adj[(j, i)] = True
    return V, adj

# ======================================================================
# METHOD 6: SLOPE-PRODUCT / RELATIVE-RISING tournament.
#   Vertices = runners. For pair (i,j) consider the RELATIVE distance function
#   D_ij(t) = ||v_i t|| - ||v_j t||. Orient i->j iff D_ij is RISING at tau*
#   (d/dt (||v_i t|| - ||v_j t||) > 0), i.e. i is pulling away from observer
#   faster than j. This is the EXACT "which triangle-wave is rising" seed from
#   the prompt, done PAIRWISE (not as a global key). slope_i - slope_j with
#   slope = +/- v. Ties (equal slope incl. sign*speed) broken by signed_dev.
#   Because slope depends on the section (frac<1/2 rising, >1/2 falling) and is
#   scaled by speed, the relation can be cyclic.
# ======================================================================
def method6_adj(S, t):
    V = list(range(len(S)))
    sp = list(S)
    sl = {i: triangle_slope(sp[i], t) for i in V}  # signed, scaled by speed
    adj = {}
    for i, j in combinations(V, 2):
        d = sl[i] - sl[j]
        if d > 0:
            iwins = True
        elif d < 0:
            iwins = False
        else:
            iwins = signed_dev(sp[i], t) > signed_dev(sp[j], t) or \
                    (signed_dev(sp[i], t) == signed_dev(sp[j], t) and sp[i] > sp[j])
        if iwins:
            adj[(i, j)] = True; adj[(j, i)] = False
        else:
            adj[(i, j)] = False; adj[(j, i)] = True
    return V, adj

# ======================================================================
# METHOD 7: GAP-MARGIN x SLOPE x WINDING-PARITY tournament. EXPLICITLY uses the
#   gap value M as input (the LRC quantity itself). For pair (i,j):
#     mu_v = ||v t|| - M (margin above the gap; =0 exactly for binding runners),
#     s_v  = sign of triangle-wave slope of ||v t|| at t (rising +1 / flat 0 / fall -1),
#     w_v  = floor(v t) (winding count).
#   base order key(v) = (s_v, -mu_v) [rising & near-binding ranks high]; orient
#   i->j iff (key_i > key_j) XOR (parity of (w_i+w_j) is odd). The winding XOR
#   injects a non-snapshot, potentially cyclic twist coupling dynamics to the gap.
# ======================================================================
def method7_adj(S, t):
    V = list(range(len(S))); sp = list(S)
    Mval = gmin(S, t)
    d = {i: nrm(sp[i] * t) for i in V}
    mu = {i: d[i] - Mval for i in V}
    sl = {i: (1 if triangle_slope(sp[i], t) > 0 else
              (-1 if triangle_slope(sp[i], t) < 0 else 0)) for i in V}
    adj = {}
    for i, j in combinations(V, 2):
        wi = (sp[i] * t).__floor__(); wj = (sp[j] * t).__floor__()
        base = (sl[i], -mu[i]) > (sl[j], -mu[j])
        flip = ((wi + wj) & 1) == 1
        iw = base ^ flip
        adj[(i, j)] = iw; adj[(j, i)] = not iw
    return V, adj

METHODS = {
    "M1_diff_sum_crossing_parity": method1_adj,
    "M2_binding_slope_lex (control)": method2_adj,
    "M3_meeting_parity_x_speed": method3_adj,
    "M4_section_winding_rps4": method4_adj,
    "M5_separation_parity": method5_adj,
    "M6_relative_rising_pairwise": method6_adj,
    "M7_gapmargin_slope_winding": method7_adj,
}

# ----------------------------------------------------------------------
# Driver: for each n in {3,4,5}, range S over primitive speed sets up to a cap,
# use ALL optimal tau (there can be several), build each method's tournament,
# bucket by iso class, compare to free set.
# ----------------------------------------------------------------------
def run(n, maxspeed, use_all_optima=True):
    V = list(range(n))
    free = free_classes(n)
    free_keys = set(free.keys())
    # map canon-key -> nice label
    def label(k):
        # recompute score/c3/H from the matrix k
        adjk = {}
        for i in range(n):
            for j in range(n):
                if i != j:
                    adjk[(i, j)] = bool(k[i][j])
        return (scores(adjk, V), count_3cycles(adjk, V), ham_path_count(adjk, V))

    results = {}
    Ssets = speed_sets(n, maxspeed)
    for name, fn in METHODS.items():
        realized = defaultdict(int)
        nontrivial_example = None
        bad = 0
        for S in Ssets:
            Mv, opts = all_optima(S)
            if Mv == 0:
                continue
            taus = opts if use_all_optima else opts[:1]
            for t in taus:
                Vv, adj = fn(S, t)
                if not is_tournament(adj, Vv):
                    bad += 1
                    continue
                k = canon(adj, Vv)
                realized[k] += 1
                if nontrivial_example is None and ham_path_count(adj, Vv) > 1:
                    nontrivial_example = (S, t, label(k))
        results[name] = {
            "realized": dict(realized),
            "n_classes_realized": len(realized),
            "nontrivial_example": nontrivial_example,
            "bad_nontournament": bad,
        }
    return free, free_keys, label, results

def report(n, maxspeed):
    free, free_keys, label, results = run(n, maxspeed)
    print("=" * 72)
    print(f"n = {n} runners, speeds in 1..{maxspeed}, "
          f"#primitive speed-sets = {len(speed_sets(n, maxspeed))}")
    print(f"FREE set (A000568): {len(free_keys)} iso classes")
    free_labels = sorted(set(label(k) for k in free_keys))
    print(f"  free (score,c3,H) signatures: {free_labels}")
    for name, r in results.items():
        rk = set(r["realized"].keys())
        realized_labels = sorted(set(label(k) for k in rk))
        Hset = sorted(set(l[2] for l in realized_labels))
        nontriv = any(l[2] > 1 for l in realized_labels)
        forbidden = free_keys - rk
        forb_labels = sorted(set(label(k) for k in forbidden))
        print("-" * 60)
        print(f"[{name}]")
        print(f"  realized iso classes: {r['n_classes_realized']} / {len(free_keys)}")
        print(f"  H values realized: {Hset}   NON-TRIVIAL (H>1 seen)? {nontriv}")
        print(f"  realized signatures: {realized_labels}")
        if r["bad_nontournament"]:
            print(f"  WARNING bad (non-tournament outputs): {r['bad_nontournament']}")
        if forbidden:
            print(f"  FORBIDDEN classes (in free, never realized): "
                  f"{len(forbidden)}")
            print(f"    forbidden signatures: {forb_labels}")
        if r["nontrivial_example"]:
            S, t, lab = r["nontrivial_example"]
            print(f"  nontrivial example: S={S} tau={t} -> (score,c3,H)={lab}")
    return results

if __name__ == "__main__":
    for n, ms in [(3, 14), (4, 14), (5, 12)]:
        report(n, ms)
        print()
    # Also a targeted run at n=5 restricted to COVERING-like / hard cores would go
    # here, but for class enumeration the broad primitive sweep is the population.

# ======================================================================
# DEEPER PROBE: is the M4 / M5 forbidden set due to LONELINESS, or just due to
# the MAP's algebra? Control = range tau over ALL grid times (not just optima),
# and over a dense rational grid, ignoring loneliness. If the same classes are
# forbidden WITHOUT the loneliness constraint, then loneliness is NOT the cause
# (the map itself can't reach them). A genuine LRC-forbidden class is one that
# free arc-flips reach AND the map reaches at SOME tau, but NEVER at a LONELY tau.
# ======================================================================
def probe_map(n, maxspeed, fn, dense_den=60):
    """Return (realized_at_optima, realized_at_anytau) iso-class signature sets.
       'anytau' ranges tau over k/d for d up to dense_den (all rationals), AND
       over each speed-set's own candidate set — ignoring whether it's lonely."""
    V = list(range(n))
    free = free_classes(n)
    def label(k):
        adjk = {}
        for i in range(n):
            for j in range(n):
                if i != j: adjk[(i, j)] = bool(k[i][j])
        return (scores(adjk, V), count_3cycles(adjk, V), ham_path_count(adjk, V))
    Ssets = speed_sets(n, maxspeed)
    at_opt = set(); at_any = set()
    # dense tau grid (shared, independent of S)
    grid = set()
    for d in range(2, dense_den + 1):
        for k in range(1, d):
            if gcd(k, d) == 1:
                grid.add(F(k, d))
    grid = sorted(grid)
    for S in Ssets:
        Mv, opts = all_optima(S)
        if Mv > 0:
            for t in opts:
                Vv, adj = fn(S, t)
                if is_tournament(adj, Vv):
                    at_opt.add(canon(adj, Vv))
        # any tau: use the candidate set of THIS S (the natural crossings) + grid
        Cset = cand(S) | set(grid)
        for t in Cset:
            if t <= 0 or t >= 1: continue
            Vv, adj = fn(S, t)
            if is_tournament(adj, Vv):
                at_any.add(canon(adj, Vv))
    return ({label(k) for k in at_opt}, {label(k) for k in at_any},
            {label(k) for k in free})

def report_probe(n, maxspeed):
    print("#" * 72)
    print(f"PROBE n={n}, speeds 1..{maxspeed}: optima(LONELY) vs any-tau vs free")
    for name in ["M4_section_winding_rps4", "M5_separation_parity",
                 "M3_meeting_parity_x_speed", "M6_relative_rising_pairwise"]:
        fn = METHODS[name]
        opt, anyt, free = probe_map(n, maxspeed, fn)
        forb_by_lonely = anyt - opt    # map CAN reach these but not at a lonely tau
        forb_by_map = free - anyt      # map can NEVER reach (algebra, not loneliness)
        print("-" * 60)
        print(f"[{name}]")
        print(f"  free                : {len(free)}")
        print(f"  reachable any-tau   : {len(anyt)}")
        print(f"  reachable at LONELY : {len(opt)}")
        print(f"  forbidden by MAP    (free - anytau)   : {len(forb_by_map)}")
        if forb_by_map: print(f"    -> {sorted(forb_by_map)}")
        print(f"  forbidden by LONELINESS (anytau - opt): {len(forb_by_lonely)}")
        if forb_by_lonely: print(f"    -> {sorted(forb_by_lonely)}")

if __name__ == "__main__":
    print()
    for n, ms in [(3, 14), (4, 12), (5, 10)]:
        report_probe(n, ms)
        print()
