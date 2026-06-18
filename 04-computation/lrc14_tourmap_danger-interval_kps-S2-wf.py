# lrc14_tourmap_danger-interval_kps-S2-wf.py
# Creative tournament-generation from LRC(14) "danger-interval" data.
# Theme: danger arcs D_v = {tau : ||v tau|| < 1/14}, lonely intervals, nesting,
# blocking, overlap order. We seek NON-OBVIOUS arc rules that yield NON-TRANSITIVE
# tournaments and ask which iso classes are FORBIDDEN under loneliness.
#
# All decisions use exact Fractions. No floats.

from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product
import sys

# ---------------- EXACT M TOOL (verbatim from task) ----------------
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
        if v > b: b = v; at = t
    return b, at

# ---------------- tournament canonicalization ----------------
# A tournament on m vertices: adj[i][j] = True means arc i->j.
# We canonicalize by brute iso over all m! relabelings: take lexicographically
# minimal flattened upper-data string. Also compute #3-cycles and score seq.

def canon(adj, m):
    best = None
    verts = list(range(m))
    for p in permutations(verts):
        # relabel: new vertex a corresponds to old p[a]; arc a->b iff adj[p[a]][p[b]]
        bits = []
        for a in range(m):
            for b in range(m):
                if a != b:
                    bits.append('1' if adj[p[a]][p[b]] else '0')
        s = ''.join(bits)
        if best is None or s < best:
            best = s
    return best

def num_3cycles(adj, m):
    c = 0
    for i, j, k in combinations(range(m), 3):
        # count cyclic triangles among i,j,k
        # arcs
        def arc(a, b): return adj[a][b]
        # cyclic if i->j->k->i or i->k->j->i
        if arc(i, j) and arc(j, k) and arc(k, i): c += 1
        elif arc(j, i) and arc(i, k) and arc(k, j): c += 1
    return c

def score_seq(adj, m):
    return tuple(sorted(sum(1 for j in range(m) if i != j and adj[i][j]) for i in range(m)))

def ham_paths(adj, m):
    # count Hamiltonian paths (directed)
    cnt = 0
    for p in permutations(range(m)):
        ok = all(adj[p[i]][p[i+1]] for i in range(m-1))
        if ok: cnt += 1
    return cnt

# Free set sizes A000568: n=1..7 -> 1,1,2,4,12,56,456
FREE = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456}

# ---------------- danger-interval primitives ----------------
# Danger arc D_v on the circle [0,1): set of tau with ||v tau|| < 1/14.
# ||v tau|| < 1/14 around each integer/v center. Centers at k/v for k=0..v-1.
# Each danger sub-interval: (k/v - 1/(14v), k/v + 1/(14v)) mod 1.
# Half-width in tau is 1/(14 v). Total danger measure of D_v = v * 2/(14 v) = 1/7.
# So EVERY runner has the same total danger measure 1/7. Interesting: distinguishing
# them must be by POSITION/PHASE, not size.

# For a fixed lonely time tau0 (where observer is lonely for the WHOLE set, gap>=1/14),
# each runner v is at fractional position frac(v*tau0); its distance to 0 is ||v tau0||.
# The "danger zone center nearest to tau0" etc.

# We'll represent each runner's local danger picture at the optimum.

def frac(x):
    r = x - int(x)
    return r + 1 if r < 0 else r

# =====================================================================
# METHOD 1: NESTING / CONTAINMENT of safe arcs around the optimum.
# Vertex = runner. At optimum tau0, runner v has phase ph_v = frac(v tau0) in [0,1).
# Its nearest danger center is the integer translate; the SAFE half-interval the
# runner currently occupies is [c_v - 1/(14 v), c_v + 1/(14 v)] is the danger, the
# safe gap between consecutive danger centers has half-length depending on v.
# Define arc i->j by comparing the *signed clearance to the next danger center in the
# direction of increasing tau* scaled by speed: who hits danger first as tau advances.
# c_v(tau0) = distance (in tau) from tau0 forward to the next danger boundary of v.
# Arc i->j iff i reaches its danger boundary BEFORE j (smaller forward clearance) ...
# that's a total order => transitive. To break transitivity, use a CIRCULAR
# comparison: i->j iff (ph_j - ph_i mod 1) < 1/2  -- the overtaking snapshot (DEAD).
# Instead Method 1 uses forward-clearance-to-danger but on the v-SCALED circle with a
# WINDOW: i->j iff j's danger center lies inside i's CURRENT safe interval mapped...
# We implement a concrete non-snapshot rule below.

def forward_clearance(v, tau0):
    # tau distance forward until ||v tau|| hits 1/14 (entering danger), within current safe window.
    # phase p = frac(v*tau0). danger when ||.|| < 1/14, i.e. p in (k-1/14, k+1/14).
    # We're at a lonely time so ||v tau0|| >= 1/14, p's distance to nearest int >= 1/14.
    p = frac(v * tau0)
    # nearest integer boundary going forward (increasing tau => increasing p mod1 at rate v)
    # next danger entry: the next integer m with p approaching m-1/14 from below, or
    # m+1/14? entering danger when p crosses into (m-1/14, m+1/14).
    # forward in p: next boundary is the start of a danger band. Bands centered at integers.
    # current p in a safe band between (a+1/14) and (a+1 -1/14) for some integer a, OR could be
    # near boundary. distance forward in p to next danger start = ((a+1 -1/14) - p) where a=floor(p).
    a = int(p)  # p in [0,1)
    # safe band is (1/14, 1-1/14) within each unit; if p < 1/14 it's in danger (shouldn't at lonely)
    # forward danger entry start at 1 - 1/14 = 13/14 within the unit (then band around integer 1)
    dist_p = (F(13, 14) - p)
    if dist_p < 0:
        dist_p += 1
    # convert to tau distance: dp = v * dtau => dtau = dp / v
    return dist_p / v

def method1_adj(S, tau0):
    m = len(S)
    fc = [forward_clearance(S[i], tau0) for i in range(m)]
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            # i->j iff i has SMALLER forward clearance (i hits danger first)
            if fc[i] < fc[j]:
                adj[i][j] = True
            elif fc[i] > fc[j]:
                adj[i][j] = False
            else:
                adj[i][j] = (S[i] < S[j])  # tie-break by speed (deterministic)
    return adj

# =====================================================================
# METHOD 2: CIRCULAR-ARC OVERLAP at a fixed phase tau0.
# Vertex = runner. Each runner has a circular danger SET D_v subset of circle [0,1)
# consisting of v little arcs. Define a CHASE order using the nearest danger arc:
# let center_v = the danger center of v nearest to tau0 (an element of {k/v}), and let
# the SIGNED offset off_v = (tau0 - center_v) in (-1/(14v), ... ) clamped to nearest.
# Arc i->j iff the nearest-danger-center of i, advanced by speeds, "captures" j:
# Concretely use the CIRCULAR betweenness of three centers -> but that's for triples.
# Implement: i->j iff  frac( (center_i - center_j) * (S[i]+S[j]) ) < 1/2.
# This mixes phase difference with a pair-dependent winding (sum of speeds) => can be
# non-transitive because the multiplier depends on the pair.

def nearest_center(v, tau0):
    # nearest k/v to tau0
    k = round(float(v * tau0))  # approximate then correct exactly
    best = None; bestc = None
    for kk in (k-1, k, k+1):
        c = F(kk, v)
        d = abs(tau0 - c)
        if best is None or d < best:
            best = d; bestc = c
    return bestc

def method2_adj(S, tau0):
    m = len(S)
    cen = [nearest_center(S[i], tau0) for i in range(m)]
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            w = S[i] + S[j]
            val = frac((cen[i] - cen[j]) * w)
            if val == F(0):
                adj[i][j] = (S[i] < S[j])
            elif val < F(1, 2):
                adj[i][j] = True
            elif val > F(1, 2):
                adj[i][j] = False
            else:
                adj[i][j] = (S[i] < S[j])
    return adj

# =====================================================================
# METHOD 3: BLOCKING / WHO-FORCES-WHOM (binding structure of THM-524).
# At the optimum tau0, runner v is BINDING iff ||v tau0|| == M(S) (it sits exactly at
# the gap). Among all runners, define arc i->j by a "danger race": as tau moves AWAY
# from tau0 in the + direction, whose normalized distance ||v tau|| drops to <1/14
# first relative to the other -- but compared via the SIGN of the derivative.
# The signed velocity of ||v tau|| at tau0 is +/- v depending on which side of the
# nearest integer the phase sits. Define sgn_v = +1 if frac(v tau0) in (0,1/2) (moving
# away upward toward 1, distance decreasing means...) Let's define:
#   d_v(tau) = ||v tau||; near tau0 it's |frac(v tau) - round|.
#   left/right slope is +-v. The runner is "rising" (distance increasing with tau) or
#   "falling". Arc i->j iff i is FALLING faster toward danger than j, measured as the
#   tau-distance for i to reach 1/14 vs j, in the +direction.
# This is similar to forward_clearance but uses the ACTUAL signed motion of ||v tau||
# (both directions of danger), giving a DIFFERENT, finer order. To make it circular we
# compare on the pair circle of size (S[i]+S[j]).
# Concretely: time for runner v to reach danger (||v tau||=1/14) moving forward:
def time_to_danger_signed(v, tau0):
    p = frac(v * tau0)
    # distance to nearest integer
    d = nrm(v * tau0)  # in [0,1/2], current ||.||, >=1/14 at lonely
    # which integer is nearest? if p<=1/2 nearest is 0 (i.e. floor), distance = p; moving
    # forward (p increases) distance from 0 increases until p=1/2.
    # if p in (1/2,1) nearest is 1, distance = 1-p; moving forward distance decreases ->
    # reaches 1/14 when 1-p = 1/14 => p = 13/14.
    if p <= F(1, 2):
        # distance increasing; danger (at this integer) is behind us. The NEXT danger
        # entering is near p=13/14 (the next integer). dp forward = 13/14 - p
        dp = F(13, 14) - p
    else:
        # distance = 1-p decreasing toward 1/14 at p=13/14
        dp = F(13, 14) - p
        if dp < 0:
            dp = dp + 1  # passed it; wrap
    if dp < 0: dp += 1
    return dp / v  # tau time

def method3_adj(S, tau0):
    # arc i->j iff on the pair-circle, i's danger-arrival precedes j's, but winding by
    # (S[i]-S[j]) parity: flip orientation when (S[i]+S[j]) is even.
    m = len(S)
    ttd = [time_to_danger_signed(S[i], tau0) for i in range(m)]
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            base = ttd[i] < ttd[j]
            if ttd[i] == ttd[j]:
                base = (S[i] < S[j])
            # parity twist: if (S[i]+S[j]) even, reverse
            if (S[i] + S[j]) % 2 == 0:
                base = not base
            adj[i][j] = base
    return adj

# =====================================================================
# METHOD 4: SECTION TOURNAMENT (on-grid residues, the project's section lens).
# At grid time tau = a/14 (gcd(a,14)=1), runner i sits in SECTION r_i = (v_i*a) mod 14.
# Vertex = runner. Arc i->j by CIRCULAR ORDER of sections with a pair-dependent twist:
# i->j iff (r_i - r_j) mod 14 is in {1,2,...,6} (i.e. r_i is "clockwise-ahead" by <=6).
# Pure circular tournament on a 14-cycle of positions => non-transitive generically!
# Ties (same section, or diff exactly 7) broken by speed. This is the classic rotational
# tournament => has 3-cycles.
def method4_adj(S, a):
    m = len(S)
    r = [(S[i] * a) % 14 for i in range(m)]
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            d = (r[i] - r[j]) % 14
            if d == 0:
                adj[i][j] = (S[i] < S[j])
            elif 1 <= d <= 6:
                adj[i][j] = True
            elif d == 7:
                adj[i][j] = (S[i] < S[j])
            else:  # 8..13
                adj[i][j] = False
    return adj

# =====================================================================
# METHOD 5: PAIR-CROSSING / BINDING-SWITCH tournament (THM-524 sawtooth).
# Vertex = runner. For each ordered pair, THM-524 says M is attained at a binding-pair
# crossing tau = k/(v_a +- v_b). Define arc i->j by WHICH side of the optimum the binding
# crossing of the pair (i,j) sits, i.e. by the SIGN of (the crossing tau for pair (i,j))
# minus tau0, combined with the local order. Specifically the pair (i,j) has a crossing
# where ||v_i tau|| = ||v_j tau||; orient i->j if at that crossing v_i is the one whose
# danger band the pair sits inside (i is "inside" j's clearance). Implement via the
# binding crossing nearest tau0: c_ij = k/(v_i+v_j) closest to tau0; arc i->j iff
# frac(v_i * c_ij) < frac(v_j * c_ij) (i is closer to an integer => "more cornered").
def nearest_crossing(vi, vj, tau0):
    d = vi + vj
    k = round(float(tau0 * d))
    best = None; bc = None
    for kk in (k-1, k, k+1):
        c = F(kk, d)
        if c < 0: continue
        dist = abs(c - tau0)
        if best is None or dist < best:
            best = dist; bc = c
    return bc

def method5_adj(S, tau0):
    m = len(S)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            c = nearest_crossing(S[i], S[j], tau0)
            ni = nrm(S[i] * c); nj = nrm(S[j] * c)
            if ni < nj:
                adj[i][j] = True
            elif ni > nj:
                adj[i][j] = False
            else:
                adj[i][j] = (S[i] < S[j])
    return adj

# =====================================================================
# METHOD 6: DIFFERENCE-WINDING tournament (overtaking on the pair circle).
# Vertex = runner. The TWO runners i,j collide/overtake on the circle at rate |v_i-v_j|.
# At tau0 their relative phase is frac((v_i - v_j) tau0). Arc i->j iff this relative
# phase is in (0,1/2). This is the "relative overtaking" — DIFFERENT from the absolute
# overtaking snapshot because it uses the DIFFERENCE speed, and the winding number
# differs per pair, which can break transitivity.
def method6_adj(S, tau0):
    m = len(S)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            if i == j: continue
            rel = frac((S[i] - S[j]) * tau0)
            if rel == F(0):
                adj[i][j] = (S[i] < S[j])
            elif rel < F(1, 2):
                adj[i][j] = True
            elif rel > F(1, 2):
                adj[i][j] = False
            else:
                adj[i][j] = (S[i] < S[j])
    return adj

# ---------------- input families ----------------
# LRC-constrained inputs: primitive speed sets with a genuine lonely optimum.
# We enumerate small primitive speed sets of size m=3,4,5 and use tau0 = optimal tau.
# Primitive = gcd of all speeds = 1. We range speeds in a bounded window.

def gen_speed_sets(m, vmax):
    out = []
    for combo in combinations(range(1, vmax+1), m):
        if gcd_list(combo) == 1:
            out.append(combo)
    return out

def gcd_list(xs):
    g0 = 0
    for x in xs: g0 = gcd(g0, x)
    return g0

# ---------------- driver ----------------
def analyze(method_name, adj_fn, uses_grid, m_values, vmax_map, sample_cap=None):
    print("="*70)
    print(f"METHOD: {method_name}   (uses_grid={uses_grid})")
    print("="*70)
    for m in m_values:
        vmax = vmax_map[m]
        sets = gen_speed_sets(m, vmax)
        realized = {}  # canon -> (h, c3, score, example)
        nontrivial_seen = False
        count = 0
        for S in sets:
            S = list(S)
            if uses_grid:
                # range over grid times a in (Z/14)* and require loneliness on grid?
                # We use all a coprime to 14; vertex=runner; build section tournament.
                inputs = [a for a in range(1, 14) if gcd(a, 14) == 1]
                for a in inputs:
                    adj = adj_fn(S, a)
                    cn = canon(adj, m)
                    if cn not in realized:
                        h = ham_paths(adj, m)
                        c3 = num_3cycles(adj, m)
                        sc = score_seq(adj, m)
                        realized[cn] = (h, c3, sc, (tuple(S), a))
                        if h > 1: nontrivial_seen = True
                    count += 1
            else:
                gap, tau0 = M(S)
                if tau0 is None: continue
                adj = adj_fn(S, tau0)
                cn = canon(adj, m)
                if cn not in realized:
                    h = ham_paths(adj, m)
                    c3 = num_3cycles(adj, m)
                    sc = score_seq(adj, m)
                    realized[cn] = (h, c3, sc, (tuple(S), tau0))
                    if h > 1: nontrivial_seen = True
                count += 1
            if sample_cap and count > sample_cap:
                break
        free = FREE[m]
        print(f"  m={m}: #speed-sets-or-inputs scanned={count}, "
              f"realized {len(realized)}/{free} iso classes, nontrivial={nontrivial_seen}")
        # list realized classes
        hs = sorted(set(v[0] for v in realized.values()))
        c3s = sorted(set(v[1] for v in realized.values()))
        print(f"     H values realized: {hs}")
        print(f"     #3-cycle values realized: {c3s}")
        if len(realized) < free:
            print(f"     *** {free - len(realized)} classes NOT realized (candidate FORBIDDEN) ***")
        # show example of a nontrivial realized class
        for cn,(h,c3,sc,ex) in sorted(realized.items()):
            if h > 1:
                print(f"     e.g. nontrivial: H={h}, c3={c3}, score={sc}, from S={ex[0]}, key={ex[1]}")
                break
    print()

if __name__ == "__main__":
    vmax_map = {3: 14, 4: 12, 5: 10}
    methods = [
        ("M1 forward-clearance order", method1_adj, False),
        ("M2 circular-center winding", method2_adj, False),
        ("M3 signed danger-arrival parity-twist", method3_adj, False),
        ("M4 section rotational (on-grid)", method4_adj, True),
        ("M5 binding-crossing cornered", method5_adj, False),
        ("M6 difference-winding overtaking", method6_adj, False),
    ]
    for name, fn, grid in methods:
        analyze(name, fn, grid, [3, 4, 5], vmax_map)
