"""
lrc_personal_space_tournament_opus_S138.py

THE PERSONAL-SPACE TOURNAMENT (owner: invent runner-tournaments from pair statistics with
creative binary cutoffs; relativity hint: from each runner's perspective all others move
relatively, and pairwise rates are exactly antisymmetric).

Object: for the 14 runners {0, v_1..v_13} (observer included), runner i's PERSONAL SPACE at
time x is PS_i(x) = (distance to nearest phase below) + (distance to nearest above) — the
sum of the two gaps flanking phase_i.  Facts:
  * CONSERVATION: sum_i PS_i(x) = 2 exactly (each circular gap is counted by both flanking
    runners).  The mean personal space is a zero-sum resource: sum_i E[PS_i] = 2.
  * RELATIVITY: PS_i is frame-invariant (translating all speeds by -v_i puts i at rest;
    pairwise approach rates are antisymmetric, distances are shared) — the tournament below
    is therefore well-defined on the SPEED SET, no frame choice needed.
  * LONELINESS: the observer is lonely at x iff both its flanking distances are >= 1/14 —
    in particular PS_0(x) > 1/7 is necessary; the anchored-tail program (gap∋0 = PS of the
    observer) is personal-space statistics.  The 2-anchor tail adds the antipode's space.

TOURNAMENT: orient i -> j  iff  E_x[PS_i] > E_x[PS_j]  (exact rationals; ties -> by speed).
This is the "who lives more spaciously on average" dominance order.  Questions probed:
  (1) is it transitive (it is a total preorder by construction — the INTERESTING data is
      the exact E[PS_i] profile and WHERE the observer ranks);
  (2) does the observer's rank / the profile shape separate tight vs loose vs record
      families;
  (3) profile symmetry: reversal E -> max+min-E permutes runners i <-> mirror(i) — is the
      profile mirror-symmetric on palindromic families (prediction: yes) and asymmetric on
      mirror pairs (measures the "spontaneous symmetry breaking" of mac-mini-S43)?

Note E[PS_i] = E[gap-left_i] + E[gap-right_i]; per subcell of the order-cell decomposition
each flanking gap of a fixed runner is affine, so E[PS_i] is exact per runner.
"""
from fractions import Fraction as F
import sys, time

sys.path.insert(0, ".")
from lrc_exact_mu_ordercells_opus_S136 import order_cells

def EPS_exact(V):
    """Exact E_x[PS_v] for each v in V (speeds, may include 0). Returns dict v -> Fraction."""
    V = sorted(V)
    tot = {v: F(0) for v in V}
    for a, b in order_cells([v for v in V if v != 0] or V):
        mid = (a + b) / 2
        # sorted phases with owners at midpoint; per cell the order is constant
        pts = sorted((((v * mid) % 1), v) for v in V)
        n = len(pts)
        # per cell, each adjacent pair (i, i+1) gap = affine c*x + d with
        # c = v_{i+1} - v_i adjusted by wrap; integrate and credit both owners
        for s in range(n):
            p1, v1 = pts[s]
            p2, v2 = pts[(s + 1) % n]
            if s < n - 1:
                c = F(v2 - v1)
                d = ((v2 * mid) % 1) - ((v1 * mid) % 1) - c * mid
            else:
                c = F(v2 - v1)
                val = ((v2 * mid) % 1) + 1 - ((v1 * mid) % 1)
                d = val - c * mid
            # integral of gap over (a,b)
            I = c * (b * b - a * a) / 2 + d * (b - a)
            tot[v1] += I
            tot[v2] += I
    return tot

def profile(V, name):
    t0 = time.time()
    eps = EPS_exact(V)
    S = sum(eps.values())
    order = sorted(eps.items(), key=lambda kv: (-kv[1], kv[0]))
    ranks = {v: r for r, (v, _) in enumerate(order)}
    obs_rank = ranks.get(0, None)
    print(f"\n {name}  (sum E[PS] = {S} {'== 2 OK' if S == 2 else '*** SUM != 2 ***'})")
    row = "   " + "  ".join(f"{v}:{float(e):.4f}" for v, e in sorted(eps.items()))
    print(row)
    print(f"   dominance order (desc): {[v for v, _ in order]}   observer rank: {obs_rank}")
    # mirror symmetry check: reversal v -> (max+min) - v on the MOVING part
    mv = [v for v in V if v != 0]
    c = max(mv) + min(mv)
    mirrored_ok = all(abs(float(eps[v]) - float(eps.get(c - v, F(-1)))) < 1e-12
                      for v in mv if (c - v) in eps)
    print(f"   moving-part mirror profile symmetric: {mirrored_ok}   [{time.time()-t0:.0f}s]")

def main():
    print("=" * 100)
    print("PERSONAL-SPACE TOURNAMENT: exact E[PS] profiles (observer 0 included; sum = 2 law)")
    print("=" * 100)
    fams = {
        "AP {0,1..13}": [0] + list(range(1, 14)),
        "record (2^4 1^4 2^4) {0,2,..}": [0, 2, 4, 6, 8, 10, 11, 12, 13, 14, 16, 18, 20, 22],
        "prim-sat 2AP+13": [0, 2, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 22, 24],
        "GW": [0] + list(range(1, 12)) + [13, 24],
        "deep well {1..12,182}": [0] + list(range(1, 13)) + [182],
        "far-elt D18 {1..12,19}": [0] + list(range(1, 13)) + [19],
        "mirror A {1..12,20}": [0] + list(range(1, 13)) + [20],
        "mirror B (reversal of A)": [0] + [1] + list(range(9, 21)),
    }
    for name, V in fams.items():
        if len(set(V)) != 14:
            print(f"\n {name}: SKIP (needs 14 distinct)"); continue
        profile(V, name)

if __name__ == "__main__":
    main()
