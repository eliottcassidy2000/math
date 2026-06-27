#!/usr/bin/env python3
"""apex7_vs_forbiddenH7_bridge_audit.py

ADVERSARIAL audit of the claimed bridge:
    "apex-7 of LRC(14)"  <==>  "forbidden H=7 (Omega = K_3)".

The seductive coincidence: 14 = 2*7, and 7 = I(K_3,2) is the smallest forbidden
H-value (THM-029/THM-200). kps-S31y's slogan says a disproof of LRC(14) would
realize the forbidden K_3. This script DEMANDS EVIDENCE for the four candidate
bridges and reports honestly which (if any) survive.

TEST 1 (winding tournament can NEVER show H=7 -- so it is NOT a discriminator):
  For LRC(n) AP configs, sweep the winding tournament T(t) over generic t and
  confirm H(T(t)) avoids 7. (This is TRUE BUT VACUOUS: H=7 is forbidden for
  EVERY tournament, so "winding avoids 7" carries zero information about
  tightness. We quantify this.)

TEST 2 (apex tied-arc count -- the REAL apex-7 phenomenon at n=14):
  At t*=1/14 for AP {0..13}, count the TIED (antipodal, distance 1/2) arcs.
  Claim (S48): exactly 7 (the diameters (i,i+7)). Verify, and check the analog
  count is exactly floor(n/2) ... and is n/2 = 7 ONLY because n=14. Is "7" here
  the SAME 7 as I(K_3,2)? Test whether the tied-arc count ever equals 3 (which
  would be the K_3 vertex count) and whether the tied arcs form a K_3.

TEST 3 (THM-577 overlap graph thresholded at the apex -- the "3 mutually
  overlapping arcs" candidate for K_3):
  Build graph G on speeds: edge p~q iff the forbidden arcs Fb(p),Fb(q) OVERLAP
  with the apex offset ON, i.e. p+q > 14 (THM-577). Ask:
   (a) Do triangles (3 pairwise p+q>14) exist? (trivially yes for big speeds.)
   (b) Is there ANY map from such an overlap-triangle to a directed-odd-cycle
       conflict graph Omega with Omega = K_3? I.e. is the overlap graph the
       SAME object as Omega(winding tournament)? Test by direct comparison.

TEST 4 (the honest control): does the winding tournament's Omega(T(t)) EVER
  equal K_3 for ANY t and ANY config? (It cannot, since that needs H=7.)
  And does the per-x conflict structure of the LRC cover ever LITERALLY produce
  3 pairwise-conflicting directed odd cycles in a tournament? Search small n.

Author: research audit, 2026-06-27.
"""
import sys
from fractions import Fraction as F
from math import comb, gcd
from itertools import combinations, permutations
from collections import Counter, defaultdict

sys.stdout.reconfigure(line_buffering=True)


# --------------------------------------------------------------------------
# Winding tournament + invariants
# --------------------------------------------------------------------------
def winding_adj(speeds, t):
    """adj[i][j]=1 iff i->j: frac((s_i - s_j) t) in (0,1/2). t a Fraction.
    Returns (adj, n_tied) where n_tied = #pairs at exact distance 1/2 (ties)."""
    n = len(speeds)
    adj = [[0] * n for _ in range(n)]
    tied = 0
    for i in range(n):
        for j in range(i + 1, n):
            f = ((speeds[i] - speeds[j]) * t) % 1
            if f == F(1, 2):
                tied += 1
                # tie-break lower index -> higher index for a genuine tournament
                adj[i][j] = 1
            elif 0 < f < F(1, 2):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj, tied


def is_tournament(adj):
    n = len(adj)
    for i in range(n):
        if adj[i][i]:
            return False
        for j in range(i + 1, n):
            if adj[i][j] + adj[j][i] != 1:
                return False
    return True


def H_count(adj):
    """#Hamiltonian paths via Held-Karp."""
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]
        for v in range(n):
            val = row[v]
            if not val or not (S & (1 << v)):
                continue
            av = adj[v]
            for u in range(n):
                if (S & (1 << u)) or not av[u]:
                    continue
                dp[S | (1 << u)][u] += val
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def directed_3cycles(adj):
    """List of frozensets {a,b,c} forming a directed 3-cycle."""
    n = len(adj)
    out = []
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            out.append(frozenset((a, b, c)))
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            out.append(frozenset((a, b, c)))
    return out


def conflict_graph_3cycles(cycles):
    """Conflict graph Omega restricted to 3-cycles: vertices=cycles, edge iff
    they SHARE a vertex (pairwise conflict). Returns (V, edges as set of frozenset)."""
    V = list(range(len(cycles)))
    E = set()
    for a, b in combinations(V, 2):
        if cycles[a] & cycles[b]:
            E.add(frozenset((a, b)))
    return V, E


def is_K3(V, E):
    return len(V) == 3 and len(E) == 3


# --------------------------------------------------------------------------
# TEST 1 + TEST 4: winding tournament sweep -- H spectrum + does Omega=K_3?
# --------------------------------------------------------------------------
def winding_sweep(speeds, denom):
    """Sweep t = a/denom for a=1..denom-1. Report H spectrum, tie events,
    and whether the 3-cycle conflict graph is EVER exactly K_3."""
    Hs = Counter()
    seven = 0
    ever_K3 = False
    ever_3cyc_allshare = 0  # configs with exactly 3 three-cycles all pairwise sharing
    tie_events = 0
    for a in range(1, denom):
        t = F(a, denom)
        adj, tied = winding_adj(speeds, t)
        if tied:
            tie_events += 1
        assert is_tournament(adj)
        H = H_count(adj)
        Hs[H] += 1
        if H == 7:
            seven += 1
        # check the 3-cycle conflict structure
        cyc = directed_3cycles(adj)
        # Omega uses ALL odd cycles; but K_3 (3 vtxs) requires exactly 3 odd
        # cycles total all pairwise sharing. 3-cycles are a subset; if there are
        # exactly 3 directed 3-cycles AND no longer odd cycles, check K_3.
        if len(cyc) == 3:
            V, E = conflict_graph_3cycles(cyc)
            if is_K3(V, E):
                # would need to also confirm no 5-cycles etc., but H tells all:
                # H=7 <=> Omega=K_3 exactly. Since H!=7 always, this never fully holds.
                ever_3cyc_allshare += 1
                if H == 7:
                    ever_K3 = True
    return Hs, seven, ever_K3, ever_3cyc_allshare, tie_events


# --------------------------------------------------------------------------
# TEST 2: apex tied-arc count for LRC(2m) AP {0..2m-1} at t=1/(2m)
# --------------------------------------------------------------------------
def apex_tied_arcs(m_apex):
    """For AP {0,...,m_apex-1} at t = 1/m_apex, count tied (distance-1/2) arcs.
    Returns (count, list of tied pairs)."""
    speeds = list(range(m_apex))
    t = F(1, m_apex)
    tied = []
    for i in range(m_apex):
        for j in range(i + 1, m_apex):
            f = ((speeds[i] - speeds[j]) * t) % 1
            if f == F(1, 2):
                tied.append((i, j))
    return len(tied), tied


# --------------------------------------------------------------------------
# TEST 3: THM-577 apex-thresholded overlap graph
# --------------------------------------------------------------------------
def overlap_ON(p, q, apex=14):
    """THM-577: forbidden arcs Fb(p),Fb(q) have the offset sum ON iff p+q>apex.
    (For coprime p<q; we just test the threshold predicate.)"""
    return p + q > apex


def overlap_graph(speeds, apex=14):
    V = list(range(len(speeds)))
    E = set()
    for a, b in combinations(V, 2):
        if overlap_ON(speeds[a], speeds[b], apex):
            E.add(frozenset((a, b)))
    return V, E


def triangles(V, E):
    out = []
    Eset = E
    for a, b, c in combinations(V, 3):
        if (frozenset((a, b)) in Eset and frozenset((b, c)) in Eset
                and frozenset((a, c)) in Eset):
            out.append((a, b, c))
    return out


# ==========================================================================
def main():
    print("#" * 74)
    print("# APEX-7  vs  FORBIDDEN-H=7  BRIDGE AUDIT")
    print("#" * 74)

    # ---------------- TEST 1 + 4 ----------------
    print("\n" + "=" * 74)
    print("TEST 1 & 4: winding tournament H-spectrum & does Omega EVER = K_3?")
    print("=" * 74)
    print("Config: observer 0 + AP runners. Sweep all generic t = a/denom.")
    configs = {
        "LRC(5) AP {0,1,2,3,4}": (list(range(5)), 2*3*5*7),
        "LRC(6) AP {0..5}":      (list(range(6)), 2*3*5*7),
        "LRC(7) AP {0..6}":      (list(range(7)), 2*2*3*3*5*7),  # 1260
        "GW-ish {0,1,3,4,7}":    ([0,1,3,4,7],   2*3*5*7),
        "non-tight {0,1,2,4,8}": ([0,1,2,4,8],   2*3*5*7),
    }
    for name, (sp, dn) in configs.items():
        Hs, seven, everK3, n3share, ties = winding_sweep(sp, dn)
        spectrum = sorted(Hs)
        allodd = all(h % 2 == 1 for h in Hs)
        print(f"\n  {name}  (denom={dn})")
        print(f"    H spectrum            : {spectrum}")
        print(f"    all H odd             : {allodd}")
        print(f"    H=7 occurrences       : {seven}   (forbidden => MUST be 0)")
        print(f"    Omega == K_3 ever     : {everK3}   (needs H=7 => MUST be False)")
        print(f"    grid pts w/ exactly 3  three-cycles all-sharing: {n3share}")
        print(f"    tie events (degenerate t): {ties}")

    print("""
  VERDICT TEST 1/4: "winding tournament avoids H=7" is TRUE but VACUOUS.
  H=7 is forbidden for EVERY tournament (THM-200), so the winding tournament of
  ANY speed set -- tight, loose, or random -- avoids 7. It cannot discriminate
  tight from loose. Omega(T(t)) is NEVER K_3 for any t, any config (would need
  H=7). So the apex obstruction is NOT literally a forbidden-H=7 event in T(t).""")

    # ---------------- TEST 2 ----------------
    print("\n" + "=" * 74)
    print("TEST 2: apex tied-arc count for AP {0..m-1} at t=1/m  (the REAL apex-7)")
    print("=" * 74)
    print("  m   #tied-arcs   = floor(m/2)?   tied pairs (sample)")
    for m in [4, 6, 8, 10, 12, 14, 16, 18]:
        cnt, pairs = apex_tied_arcs(m)
        pred = m // 2 if m % 2 == 0 else 0
        flag = "<== m=14 gives 7" if m == 14 else ""
        print(f"  {m:2d}   {cnt:3d}          {pred:3d} ({cnt==pred})        "
              f"{pairs[:3]}{'...' if len(pairs)>3 else ''}  {flag}")
    print("""
  The apex tied-arc count is m/2 for even m. At m=14 it is 7. This "7" is
  n/2 = 14/2, the count of ANTIPODAL DIAMETER PAIRS -- it is the SAME 7 only
  because the LRC denominator is 14. It is NOT I(K_3,2)=7 and it is NOT 3 (the
  vertex count of K_3). The tied arcs (i,i+7) form a PERFECT MATCHING on 14
  vertices (7 disjoint edges), which is the OPPOSITE of a K_3 (a matching has
  NO triangles). So "apex tied-arc count = 7" and "K_3 has I=7" are different
  sevens that happen to coincide numerically at this one denominator.""")

    # Is the tied-arc set ever a triangle / K_3?
    cnt14, pairs14 = apex_tied_arcs(14)
    # build graph of tied arcs and look for triangles
    Vt = list(range(14))
    Et = set(frozenset(p) for p in pairs14)
    tri = triangles(Vt, Et)
    print(f"  Apex tied-arc graph at m=14: {len(Et)} edges (a perfect matching), "
          f"triangles = {len(tri)} (a matching is triangle-free).")

    # ---------------- TEST 3 ----------------
    print("\n" + "=" * 74)
    print("TEST 3: THM-577 apex-thresholded overlap graph (p+q>14) -- is it K_3-related?")
    print("=" * 74)
    test_sets = {
        "AP {1..13} (tight)":            list(range(1, 14)),
        "GW {1..11,13,24}":              list(range(1, 12)) + [13, 24],
        "small {1,2,3}":                 [1, 2, 3],
        "three big {12,13,14}":          [12, 13, 14],
        "three mid {7,8,9}":             [7, 8, 9],
    }
    for name, sp in test_sets.items():
        V, E = overlap_graph(sp, apex=14)
        tri = triangles(V, E)
        print(f"\n  {name}: speeds={sp}")
        print(f"    overlap edges (p+q>14): {len(E)}")
        print(f"    overlap triangles      : {len(tri)}  "
              f"{'(e.g. ' + str([tuple(sp[i] for i in t) for t in tri[:2]]) + ')' if tri else ''}")
    print("""
  VERDICT TEST 3: the apex-thresholded overlap graph (edge iff p+q>14) has
  triangles whenever 3 speeds pairwise sum > 14 (e.g. {12,13,14}). But this is
  an UNDIRECTED OVERLAP graph on SPEEDS, not a directed-odd-cycle conflict graph
  Omega on a tournament. Its triangles are NOT directed-odd-cycle K_3 atoms.
  Crucially: a triangle here = "3 forbidden arcs that pairwise overlap a lot" =
  speeds whose forbidden regions overlap = a COVERING-DENSE cluster, which is
  the GOOD (lonely-helping) direction, not a conflict obstruction. There is no
  morphism shown from "overlap triangle" to "Omega=K_3". The K_3=7 link is
  numerical, not structural, on this object too.""")

    print("\n" + "#" * 74)
    print("# OVERALL: the three candidate bridges all FAIL to realize a literal")
    print("# H=7 / Omega=K_3 event. The genuine apex-7 phenomenon is the 7 TIED")
    print("# DIAMETER ARCS (a perfect matching, order-2 antipodal symmetry), NOT")
    print("# a 3-clique. '7 = n/2' and '7 = I(K_3,2)' are distinct sevens.")
    print("#" * 74)


if __name__ == "__main__":
    main()
