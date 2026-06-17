#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angleD_overtaking_tournament_mac-mini-2026-06-17-S3.py   (ANGLE D)

THE TOURNAMENT ANALOGY PROPER for LRC(14), via the REGIONS/SECTIONS reframe.

As tau: 0 -> 1 the runners overtake each other (swap cyclic order) at crossings
tau = k/(v_i - v_j); the wiring diagram is a periodic sequence of transpositions.

This script does three concrete things and reports the SHARPEST honest tournament/OCF
bridge:

  (1) CYCLIC ORDER + GAP SEQUENCE at the lonely time tau*. Observer loneliness
      = one gap (around section 0) >= 1/14.

  (2) THE OVERTAKING TOURNAMENT T*(tau): i -> j iff runner i is "ahead" of runner j
      at tau* (by fractional position on the circle, measured from the observer at 0).
      Speeds are DISTINCT so positions are distinct at a generic tau* => a genuine
      tournament. We compute:
        - H(T*) = # directed Hamiltonian paths (Redei: must be ODD; checks parity)
        - score sequence, whether T* is transitive
        - the OCF conflict graph Omega(T*) (vertices = directed odd 3-cycles, edges =
          share a vertex), and H_OCF = I(Omega,2) = # Ham paths via the OCF.
      KEY ARITHMETIC FACT to test: the "ahead at tau*" order is literally
      i -> j iff frac(v_i tau*) > frac(v_j tau*); for tau* = p/q this is a LINEAR
      ORDER, hence T* is TRANSITIVE => H(T*) = 1, Omega EMPTY. So the position
      tournament is trivial. We then test the NON-trivial alternative:

  (2b) THE OVERTAKING-PARITY TOURNAMENT. i -> j iff the number of times i has
      overtaken j on [0, tau*] has a fixed parity. This is the crossing count
      c_ij = floor((v_i - v_j) tau*) style count; the parity gives a relation
      that need NOT be transitive. Test H, Redei parity, Omega.

  (3) THE SECTION CONFLICT GRAPH. "each runner its own section" (perfect SDR) is a
      PROPER COLORING of a conflict graph whose vertices are runners and edges =
      "two runners are forced into the SAME 14-section at tau*" (i.e. floor(14*pos)
      collides). Compute independence number alpha and chromatic number chi for the
      sample configs, and compare to the project's Omega (where alpha / I(.,2)=H).
      Genuine bridge test: is "perfect SDR" <=> "this conflict graph is edgeless"
      <=> (analogue of) "Omega edgeless <=> H = 3^alpha (ideal gas)"?

All exact rationals via Fraction. stdlib only.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
from collections import Counter, defaultdict

# ----------------------------------------------------------------------
# exact M tool (validated, from prompt)
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

# ----------------------------------------------------------------------
# tournament tools
# ----------------------------------------------------------------------
def ham_path_count(adj, V):
    """# directed Hamiltonian paths. adj[(i,j)]=True iff i->j. V = list of vertices."""
    n = len(V); idx = {v: k for k, v in enumerate(V)}
    # DP over subsets ending at each vertex
    full = (1 << n) - 1
    # dp[mask][last] = # paths visiting 'mask', ending at 'last'
    dp = [[0] * n for _ in range(1 << n)]
    for k in range(n): dp[1 << k][k] = 1
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
    return {v: sum(1 for u in V if u != v and adj[(v, u)]) for v in V}

def is_transitive(adj, V):
    for a in V:
        for b in V:
            for c in V:
                if a != b and b != c and a != c:
                    if adj[(a, b)] and adj[(b, c)] and not adj[(a, c)]:
                        return False
    return True

def odd_cycles_up_to(adj, V, maxlen=5):
    """all directed cycles of odd length 3,5,... up to maxlen, as frozenset of vertices.
    (Omega vertices are the directed odd cycles; for the bridge we record vertex-sets.)
    We list distinct directed cycles; represent each by its vertex frozenset PLUS a
    canonical rotation tuple to distinguish cycles on the same vertex set."""
    cycles = []
    Vset = list(V)
    for L in range(3, maxlen + 1, 2):
        for combo in combinations(Vset, L):
            # all directed cycles on this vertex set
            for perm in permutations(combo[1:]):
                cyc = (combo[0],) + perm
                ok = all(adj[(cyc[k], cyc[(k + 1) % L])] for k in range(L))
                if ok:
                    # canonical: rotate to smallest start, keep direction
                    rots = [tuple(cyc[i:] + cyc[:i]) for i in range(L)]
                    canon = min(rots)
                    cycles.append(canon)
    return sorted(set(cycles))

def omega_graph(cycles):
    """conflict graph: vertices = cycles, edge iff share >=1 vertex."""
    n = len(cycles)
    supp = [frozenset(c) for c in cycles]
    adj = [set() for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if supp[a] & supp[b]:
                adj[a].add(b); adj[b].add(a)
    return adj

def indep_poly_at(adj, n, x):
    """I(G,x) = sum_{indep U} x^|U| (brute)."""
    val = 0
    for mask in range(1 << n):
        verts = [i for i in range(n) if (mask >> i) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a + 1, len(verts)):
                if verts[b] in adj[verts[a]]:
                    ok = False; break
            if not ok: break
        if ok:
            val += x ** len(verts)
    return val

# ----------------------------------------------------------------------
# graph independence / chromatic (brute, for the conflict graph)
# ----------------------------------------------------------------------
def indep_number(adj, n):
    best = 0
    for mask in range(1 << n):
        verts = [i for i in range(n) if (mask >> i) & 1]
        if len(verts) <= best: continue
        ok = True
        for a in range(len(verts)):
            for b in range(a + 1, len(verts)):
                if verts[b] in adj[verts[a]]:
                    ok = False; break
            if not ok: break
        if ok: best = len(verts)
    return best

def chromatic_number(adj, n):
    """greedy lower bound + exact via increasing k (small n)."""
    if n == 0: return 0
    edges = [(a, b) for a in range(n) for b in adj[a] if a < b]
    if not edges: return 1
    for k in range(1, n + 1):
        # try to k-color via backtracking
        color = [-1] * n
        def bt(v):
            if v == n: return True
            for c in range(k):
                if all(color[u] != c for u in adj[v]):
                    color[v] = c
                    if bt(v + 1): return True
                    color[v] = -1
            return False
        if bt(0): return k
    return n

# ======================================================================
print("=" * 80)
print("ANGLE D: the OVERTAKING TOURNAMENT and SECTION CONFLICT GRAPH for LRC(14)")
print("=" * 80)

CONFIGS = [
    ("TIGHT AP {1..13}", list(range(1, 14))),
    ("hardcore 84 {1..11,13,84}", [1,2,3,4,5,6,7,8,9,10,11,13,84]),
    ("drop6 u98 {1..5,7..13,98}", [1,2,3,4,5,7,8,9,10,11,12,13,98]),
]

# ----------------------------------------------------------------------
# (1) cyclic order + gap sequence at lonely time
# ----------------------------------------------------------------------
print("\n" + "#" * 80)
print("(1) CYCLIC ORDER + GAP SEQUENCE at the lonely time tau*")
print("#" * 80)
for name, S in CONFIGS:
    Mv, tau = M(S)
    pos = sorted(((v * tau) % 1, v) for v in S)
    order = [v for _, v in pos]
    fr = [p for p, _ in pos]
    # gaps around the circle, INCLUDING the wrap gap that contains the observer at 0
    gaps = []
    allp = [F(0)] + fr + [F(1)]  # observer at 0 and 1 (same point)
    # actually the observer is a STATIONARY runner at 0; loneliness = nearest runner dist
    # gap sequence between consecutive runner positions:
    aug = fr + [fr[0] + 1]
    seq = [aug[k+1] - aug[k] for k in range(len(fr))]
    # the observer sits at 0; distance to nearest runner:
    near = min(min(p, 1 - p) for p in fr)
    print(f"\n  {name}: M={Mv}, tau*={tau}")
    print(f"    cyclic order (by position): {order}")
    print(f"    gap sequence between runners: {[str(x) for x in seq]}")
    print(f"    observer (sec 0) nearest-runner dist = {near} = {float(near):.5f}  (lonely iff >= 1/14={float(F(1,14)):.5f})")
    print(f"    max runner-gap = {max(seq)} = {float(max(seq)):.5f}; observer's clear band straddles 0")

# ----------------------------------------------------------------------
# (2) the POSITION tournament (ahead at tau*) -- predicted TRANSITIVE
# ----------------------------------------------------------------------
print("\n" + "#" * 80)
print("(2) POSITION TOURNAMENT T*: i->j iff frac(v_i tau*) > frac(v_j tau*)")
print("#" * 80)
for name, S in CONFIGS:
    Mv, tau = M(S)
    V = list(S)
    fpos = {v: (v * tau) % 1 for v in V}
    adj = {}
    for a in V:
        for b in V:
            if a != b:
                adj[(a, b)] = fpos[a] > fpos[b]
    trans = is_transitive(adj, V)
    Hv = ham_path_count(adj, V) if len(V) <= 14 else None
    print(f"\n  {name}: transitive? {trans}   H(T*) = {Hv}")
    if trans:
        print("    -> position order is a TOTAL ORDER => transitive tournament => H=1, Omega EMPTY.")

# ----------------------------------------------------------------------
# (2b) the OVERTAKING-PARITY tournament (NON-trivial)
# i->j iff overtaking count parity. crossings of i,j on [0,tau*]:
#   relative position of i wrt j moves at speed (v_i - v_j); the # of full laps
#   of (v_i - v_j) by time tau* is floor(|v_i-v_j| tau*) (plus boundary). We use
#   c_ij = number of integers in (0, (v_i-v_j) tau*) ~ # times i laps j.
# Direction by sign of (frac of (v_i-v_j)tau* ) i.e. who is currently ahead in the
# relative frame; parity twist by c_ij.
# ----------------------------------------------------------------------
print("\n" + "#" * 80)
print("(2b) OVERTAKING-PARITY TOURNAMENT: i->j iff (current relative lead XOR crossing-parity)")
print("#" * 80)

def crossing_count(vi, vj, tau):
    """# of overtakings of j by i on (0, tau]: integer crossings of (vi-vj)*t.
    = floor((vi-vj)*tau) if vi>vj else floor((vj-vi)*tau) for the OTHER direction.
    We count overtakes = floor of |delta|*tau (number of times relative phase wrapped)."""
    d = abs(vi - vj)
    return int(d * tau)  # floor of nonneg rational

def parity_tournament(S, tau):
    V = list(S)
    adj = {}
    for a in V:
        for b in V:
            if a != b:
                # base direction: who is currently ahead in relative frame
                rel = nrm((a - b) * tau)  # in [0,1/2]; but we need a directed sign
                d = (a * tau - b * tau) % 1  # frac in [0,1)
                ahead = d < F(1, 2)  # a ahead of b in relative circular sense
                cc = crossing_count(a, b, tau)
                # twist by crossing parity
                flag = ahead ^ (cc % 2 == 1)
                adj[(a, b)] = flag
    # enforce tournament antisymmetry: ensure exactly one of (a,b),(b,a)
    for a in V:
        for b in V:
            if a < b:
                if adj[(a, b)] == adj[(b, a)]:
                    # tie-break by speed to keep it a tournament
                    adj[(a, b)] = a > b
                    adj[(b, a)] = not adj[(a, b)]
    return adj, V

for name, S in CONFIGS:
    Mv, tau = M(S)
    adj, V = parity_tournament(S, tau)
    trans = is_transitive(adj, V)
    Hv = ham_path_count(adj, V)
    sc = scores(adj, V)
    print(f"\n  {name}: tau*={tau}")
    print(f"    transitive? {trans}   H(parity tournament) = {Hv}  (Redei: must be ODD -> {Hv % 2 == 1})")
    print(f"    score sequence: {sorted(sc.values())}")
    if len(V) <= 13:
        cyc = odd_cycles_up_to(adj, V, maxlen=3)  # 3-cycles only for tractability
        om = omega_graph(cyc)
        nv = len(cyc)
        if nv <= 22:
            Hocf = indep_poly_at(om, nv, 2)
            edges = sum(len(a) for a in om) // 2
            print(f"    Omega(3-cycles): {nv} vertices, {edges} edges; I(Omega,2)={Hocf}")
        else:
            print(f"    Omega(3-cycles): {nv} vertices (too many for I(.,2) brute)")

# ----------------------------------------------------------------------
# (3) THE SECTION CONFLICT GRAPH
# vertices = runners; edge i~j iff at tau* they fall in the SAME 14-section
# (floor(14 * frac(v tau*)) collides). "perfect SDR" = edgeless = proper 1-coloring works.
# We compute alpha (max set of runners in distinct sections) and chi (min # of
# "section reuses"? no: chi = min colors so adjacent runners differ; with this
# adjacency, a proper coloring groups runners so same color => allowed... )
# Actually the natural object: SECTION ASSIGNMENT is a coloring of runners by
# section number (14 colors). "each runner its own section" = this coloring is a
# proper coloring of the COMPLETE graph restricted appropriately. Cleaner:
# Build conflict graph G_sec: i~j iff section(i)==section(j). Then # sections used
# = # connected components? No: i~j iff SAME section, so connected comps = sections.
# alpha(G_sec) = # distinct sections used (one rep per section). Perfect SDR iff
# G_sec is EDGELESS iff alpha = 13.
# ----------------------------------------------------------------------
print("\n" + "#" * 80)
print("(3) SECTION CONFLICT GRAPH G_sec: i~j iff same 14-section at tau*")
print("#" * 80)
for name, S in CONFIGS:
    Mv, tau = M(S)
    V = list(S)
    sec = {v: int((float((v * tau) % 1)) * 14) % 14 for v in V}
    # exact section via Fraction:
    sec = {v: int(((v * tau) % 1) * 14) for v in V}
    n = len(V)
    adj = [set() for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if sec[V[a]] == sec[V[b]]:
                adj[a].add(b); adj[b].add(a)
    edges = sum(len(x) for x in adj) // 2
    alpha = indep_number(adj, n)  # max # runners in DISTINCT sections
    chi = chromatic_number(adj, n)
    distinct_secs = len(set(sec.values()))
    sec_clear = 0 not in sec.values()
    occ = Counter(sec.values())
    print(f"\n  {name}: tau*={tau}")
    print(f"    section of each runner: {dict(sorted(sec.items()))}")
    print(f"    distinct sections used = {distinct_secs}/{n};  section 0 clear? {sec_clear}")
    print(f"    G_sec: {edges} edges; alpha (max distinct-section set) = {alpha}; chi = {chi}")
    print(f"    PERFECT SDR (edgeless)? {edges == 0}   (alpha={alpha}==n={n}? {alpha==n})")

# ----------------------------------------------------------------------
# (3b) ADJACENT-SECTION conflict graph (the band is 1/14, so being one section
# apart can still be "too close" if within 1/14). edge i~j iff |frac(v_i tau)-frac(v_j tau)| < 1/14.
# This is the REAL loneliness conflict: runners within the danger band of each other.
# It's the unit-distance-on-circle graph at threshold 1/14.
# ----------------------------------------------------------------------
print("\n" + "#" * 80)
print("(3b) PROXIMITY CONFLICT GRAPH G_prox: i~j iff circular dist(frac v_i tau, frac v_j tau) < 1/14")
print("#" * 80)
def cdist(x, y):
    d = abs(x - y) % 1
    return min(d, 1 - d)
for name, S in CONFIGS:
    Mv, tau = M(S)
    V = list(S)
    fp = {v: (v * tau) % 1 for v in V}
    n = len(V)
    adj = [set() for _ in range(n)]
    for a in range(n):
        for b in range(a + 1, n):
            if cdist(fp[V[a]], fp[V[b]]) < F(1, 14):
                adj[a].add(b); adj[b].add(a)
    edges = sum(len(x) for x in adj) // 2
    alpha = indep_number(adj, n)
    chi = chromatic_number(adj, n)
    print(f"\n  {name}: tau*={tau}  G_prox: {edges} edges; alpha={alpha}; chi={chi}")
    print(f"    edgeless (all runners >=1/14 apart)? {edges == 0}")

print("\n" + "=" * 80)
print("VERDICT printed in structured summary.")
print("=" * 80)
