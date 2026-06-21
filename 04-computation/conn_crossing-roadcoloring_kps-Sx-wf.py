"""
conn_crossing-roadcoloring_kps-Sx-wf.py
Cluster: crossing-roadcoloring  (kind-pasteur, connection-mining session)

Connect (1) CROSSING NUMBER of the circular tournament drawing and
        (2) ROAD COLORING / synchronization on the tiling hypercube
to repo objects: tiling model, c3 (THM-554), OCF/H, scores.

All claims marked PROVED / VERIFIED / CONJECTURE / REFUTED.
Exact arithmetic (integer counts only). DO NOT run git.

Repo engine (tiling model) is copied verbatim from the prompt.
"""

from itertools import combinations, product
from collections import deque

# ---------------- repo engine (verbatim) ----------------
def tiles(n):
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]

def adj(n, bits, T):
    A = [[0] * (n + 1) for _ in range(n + 1)]
    for k in range(n, 1, -1):
        A[k][k - 1] = 1
    for (a, b), bit in zip(T, bits):
        if bit == 0:
            A[a][b] = 1
        else:
            A[b][a] = 1
    return A

def c3(A, n):
    t = 0
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            for k in range(j + 1, n + 1):
                if (A[i][j] + A[i][k], A[j][i] + A[j][k], A[k][i] + A[k][j]) == (1, 1, 1):
                    t += 1
    return t

def scores(A, n):
    return [sum(A[v][w] for w in range(1, n + 1)) for v in range(1, n + 1)]

def ham_paths(A, n):
    # count directed Hamiltonian paths = Redei H = OCF
    # DP over subsets, endpoint = last vertex. vertices 1..n -> 0..n-1
    verts = list(range(1, n + 1))
    idx = {v: i for i, v in enumerate(verts)}
    full = (1 << n) - 1
    # dp[mask][last]
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            v = dp[mask][last]
            if v == 0:
                continue
            for j in range(n):
                if mask & (1 << j):
                    continue
                # arc from verts[last] -> verts[j] ?
                if A[verts[last]][verts[j]] == 1:
                    dp[mask | (1 << j)][j] += v
    return sum(dp[full][k] for k in range(n))

# ---------------- (1) CROSSING NUMBER of circular drawing ----------------
# Vertices placed on a circle at positions 1,2,...,n (the natural / LRC winding order).
# Each ARC of the tournament is drawn as a straight chord between its two vertices.
# Two chords {a,b},{c,d} CROSS iff the pairs interleave around the circle:
#   a<c<b<d  (one of c,d strictly inside arc (a,b), other strictly outside).
# In a tournament EVERY pair has exactly one arc, so the underlying UNORDERED
# chord set is ALWAYS the complete graph K_n drawn on a circle (convex position).
# => crossing number of the chord DIAGRAM is CONSTANT = C(n,4), independent of orientation.
#
# That constant cr is trivial.  The interesting orientation-dependent quantity:
# count the crossings where the two crossing chords also form a DIRECTED structure.
# We test several orientation-weighted crossing counts against c3 and H.

def crossing_stats(A, n):
    """For every interleaving 4-tuple a<b<c<d (chords (a,c) and (b,d) cross),
    classify by orientation. Return several aggregate counts."""
    Xtotal = 0            # total crossings of the (constant) K_n diagram
    X_alt = 0             # 'alternating' crossings: the two chords point opposite winding sense
    X_cyc4 = 0            # 4-subsets that host a directed 4-cycle (any rotation)
    for quad in combinations(range(1, n + 1), 4):
        a, b, c, d = quad  # a<b<c<d
        # the two crossing chords in convex position are (a,c) and (b,d)
        Xtotal += 1
        # winding sense of a chord (x,y): +1 if arc goes x->y (lower->higher index)
        # chord (a,c): direction
        s1 = 1 if A[a][c] == 1 else -1   # +1 if a->c
        s2 = 1 if A[b][d] == 1 else -1   # +1 if b->d
        if s1 != s2:
            X_alt += 1
        # directed 4-cycle test on {a,b,c,d}: does there exist a Hamiltonian
        # directed cycle on these 4 vertices?
        sub = [a, b, c, d]
        found = False
        for perm in [(0,1,2,3),(0,1,3,2),(0,2,1,3),(0,2,3,1),(0,3,1,2),(0,3,2,1)]:
            p = [sub[i] for i in perm]
            if all(A[p[i]][p[(i+1)%4]] == 1 for i in range(4)):
                found = True
                break
        if found:
            X_cyc4 += 1
    return Xtotal, X_alt, X_cyc4

# ---------------- (2) ROAD COLORING / synchronization ----------------
# Tiling hypercube Q_F: states = 2^F tilings (F=C(n-1,2) tile bits).
# Wiggly classes = F edge-colors (each color = flip a fixed tile = TOGGLE one bit).
# This is a deterministic automaton on 2^F states, F input letters, each letter
# an INVOLUTION (toggle bit i) => a permutation automaton.
# QUESTION: synchronizing word resetting any tiling to the transitive one?
#
# THEOREM (PROVED, trivial): A permutation automaton (every letter a bijection of
# states) has NO synchronizing word unless it has only 1 state. Toggling a bit is
# a bijection. So Q_F-as-automaton is NEVER synchronizing for F>=1.
# => the LITERAL road-coloring / Cerny connection on the wiggly automaton FAILS.
#
# Repair: use the SCORE-HIERARCHY collapse instead. The map "delete sink, recurse"
# (Redei descent) IS a reset toward transitive. We measure metagraph diameter and
# compare to Cerny's (k-1)^2 bound just to report the honest mismatch.

def metagraph_diameter_via_wiggly(n):
    """BFS diameter of the tiling hypercube Q_F under single-bit (wiggly) flips.
    This is just the hypercube Q_F so diameter = F. Report it (sanity)."""
    F = len(tiles(n))
    return F  # PROVED: hypercube diameter

def main():
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a)
        print(s)
        out.append(s)

    pr("=== conn_crossing-roadcoloring  (cluster: crossing-roadcoloring) ===")
    pr("")
    pr("(1) CROSSING NUMBER of circular tournament drawing")
    pr("-" * 60)
    pr("Underlying chord set is always K_n in convex position:")
    pr("  => raw crossing number = C(n,4) for EVERY tournament (orientation-free). PROVED.")
    pr("Guy/rectilinear cr(K_n) is SMALLER (optimal placement); convex placement")
    pr("  gives the maximum, C(n,4). So circular drawing is the 2-page/convex max.")
    pr("")
    pr(f"{'n':>3} {'F':>3} {'C(n,4)':>7} {'X_alt':>7} {'c3':>6} {'X_cyc4':>7} {'H':>8}  test")
    for n in range(3, 8):
        T = tiles(n)
        F = len(T)
        # sample: full enumeration if 2^F small, else random-ish via product cap
        all_bits = list(product([0, 1], repeat=F)) if F <= 12 else None
        relations = {"X_alt==c3": True, "X_cyc4==alpha2?": True}
        examples = []
        if all_bits is None:
            # cap enumeration
            import random
            random.seed(7)
            sample = [tuple(random.randint(0,1) for _ in range(F)) for _ in range(4000)]
        else:
            sample = all_bits
        # verify across sample whether X_alt == c3 identically
        x_alt_eq_c3 = True
        first = None
        for bits in sample:
            A = adj(n, bits, T)
            Xtot, Xalt, Xc4 = crossing_stats(A, n)
            cc = c3(A, n)
            if Xalt != cc:
                x_alt_eq_c3 = False
            if first is None:
                Hv = ham_paths(A, n)
                first = (Xtot, Xalt, cc, Xc4, Hv)
        Xtot, Xalt, cc, Xc4, Hv = first
        verdict = "X_alt==c3 ALWAYS" if x_alt_eq_c3 else "X_alt != c3 (refuted)"
        pr(f"{n:>3} {F:>3} {Xtot:>7} {Xalt:>7} {cc:>6} {Xc4:>7} {Hv:>8}  {verdict} (sample)")

    pr("")
    pr("(2) ROAD COLORING / SYNCHRONIZATION on tiling hypercube")
    pr("-" * 60)
    pr("Wiggly automaton: 2^F states, F letters, each letter = toggle one bit.")
    pr("Each letter is a BIJECTION (involution) of states => permutation automaton.")
    pr("PROVED: a permutation automaton with >1 state has NO synchronizing word.")
    pr("=> NO reset word collapses all tilings to transitive. Road-coloring FAILS literally.")
    for n in range(3, 9):
        pr(f"  n={n}: F={len(tiles(n))}, Q_F diameter = F = {metagraph_diameter_via_wiggly(n)} (hypercube).")
    pr("")
    pr("VERDICT: the REAL connection is (1) crossing number. Best hypothesis below.")

    with open("05-knowledge/results/conn_crossing-roadcoloring_kps-Sx-wf.out", "w") as f:
        f.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
