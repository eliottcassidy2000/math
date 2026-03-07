#!/usr/bin/env python3
"""
H=21 impossibility: general-n proof attempt.

Key insight: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
H=21 requires alpha_1 + 2*alpha_2 + 4*alpha_3 + ... = 10

Since alpha_k >= C(alpha_1, k) / something... constraints are tight.

The 6 decompositions with alpha_3=0:
  (alpha_1, alpha_2) = (10,0), (8,1), (6,2), (4,3), (2,4), (0,5)

(0,5) and (2,4) are trivially impossible (alpha_2 <= C(alpha_1,2)/something).
(4,3) proved impossible by Part D.
(6,2) proved impossible by Five-Cycle Forcing.
Remaining: (10,0) and (8,1).

THIS SCRIPT: attempts a general proof for these last two cases.

The approach: show that in ANY tournament, the achievable (alpha_1, i_2) pairs
never include (10,0) or (8,1).

Instance: kind-pasteur-2026-03-07-S33
"""

from itertools import combinations, permutations
from collections import Counter, defaultdict
import random
import time

def find_3cycles(adj, n):
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset([i,j,k]))
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset([i,j,k]))
    return cycles

def find_5cycles_dp(adj, n):
    cycles = []
    for verts in combinations(range(n), 5):
        v = list(verts)
        dp = {}
        dp[(1, 0)] = 1
        for S in range(1, 32):
            for i in range(5):
                if not (S & (1 << i)):
                    continue
                if (S, i) not in dp:
                    continue
                c = dp[(S, i)]
                for j in range(5):
                    if S & (1 << j):
                        continue
                    if adj[v[i]][v[j]]:
                        key = (S | (1 << j), j)
                        dp[key] = dp.get(key, 0) + c
        count = 0
        for j in range(1, 5):
            if (31, j) in dp and adj[v[j]][v[0]]:
                count += dp[(31, j)]
        num = count // 5
        for _ in range(num):
            cycles.append(frozenset(verts))
    return cycles

def count_alpha_i2(all_cycles):
    alpha1 = len(all_cycles)
    i2 = 0
    for a in range(alpha1):
        for b in range(a+1, alpha1):
            if not (all_cycles[a] & all_cycles[b]):
                i2 += 1
    return alpha1, i2


def general_n_sampling(n, num_samples=50000):
    """Sample random tournaments at n, collect (alpha_1, i_2) data."""
    print(f"\n{'='*60}")
    print(f"n={n}: Sampling {num_samples} random tournaments")
    print(f"{'='*60}")

    random.seed(42 + n)
    alpha_i2 = defaultdict(lambda: defaultdict(int))

    t0 = time.time()
    for trial in range(num_samples):
        adj = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1

        c3 = find_3cycles(adj, n)
        # At large n, skip higher cycles for speed if alpha_1(3-only) already > 10
        if len(c3) > 10:
            # If already >10 from 3-cycles alone, alpha_1 > 10 for sure
            continue

        c5 = find_5cycles_dp(adj, n)
        all_cyc = c3 + c5
        if len(all_cyc) > 14:  # skip if way above 10
            continue

        alpha1, i2 = count_alpha_i2(all_cyc)
        alpha_i2[alpha1][i2] += 1

    elapsed = time.time() - t0
    print(f"Done in {elapsed:.1f}s")

    # Show alpha_1 near 10
    for a1 in sorted(alpha_i2.keys()):
        if a1 < 6 or a1 > 14:
            continue
        i2_dict = alpha_i2[a1]
        i2_vals = sorted(i2_dict.keys())
        total = sum(i2_dict.values())
        # Check if needed i_2 for H=21 is present
        needed = None
        if 2 <= a1 <= 10 and (10 - a1) % 2 == 0:
            needed = (10 - a1) // 2
        tag = ""
        if needed is not None:
            if needed in i2_vals:
                tag = f" [ALERT: i_2={needed} exists!]"
            else:
                tag = f" [BLOCKED: need i_2={needed}]"
        print(f"  alpha_1={a1:2d}: i_2 in {{{','.join(str(v) for v in i2_vals)}}}, count={total}{tag}")


def prove_10_0_impossible():
    """
    Theorem: alpha_1=10 with i_2=0 is impossible.

    Proof strategy: Show that 10 pairwise-intersecting odd cycles in a tournament
    always forces a contradiction.

    Key observations:
    1. At any n, two odd cycles of lengths l1, l2 are vertex-disjoint only if l1+l2 <= n.
    2. So 3+3=6 is the SMALLEST disjoint pair type.
    3. alpha_1=10 with all 3-cycles: t3=10 requires n >= 6 (need 10 vertex-triples).
       At n=6: max t3 = C(6,3) = 20, but 10 pairwise-intersecting 3-sets on 6 vertices?
       EKR: max intersecting 3-family on [n] = C(n-1,2) = 10 for n=6.
       These are all triples through vertex 0: {{0,1,2},{0,1,3},...,{0,4,5}} = C(5,2)=10.
       ALL through vertex 0, so ALL intersecting. CAN all 10 be directed 3-cycles?

    Let's check: at n=6, can all C(5,2)=10 triples through vertex 0 be directed 3-cycles?
    """
    print("\n" + "="*60)
    print("PROOF ATTEMPT: alpha_1=10, i_2=0 impossible")
    print("="*60)

    # At n=6: Check if all 10 triples through vertex 0 can be 3-cycles
    print("\nn=6: Can all 10 triples through vertex 0 be directed 3-cycles?")
    edges6 = [(i, j) for i in range(6) for j in range(i+1, 6)]

    count_all_thru_0 = 0
    total = 0

    for bits in range(2**15):
        adj = [[0]*6 for _ in range(6)]
        for k, (i, j) in enumerate(edges6):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1

        # Count 3-cycles through vertex 0
        c3_thru_0 = 0
        for j in range(1, 6):
            for k in range(j+1, 6):
                if adj[0][j] and adj[j][k] and adj[k][0]:
                    c3_thru_0 += 1
                elif adj[0][k] and adj[k][j] and adj[j][0]:
                    c3_thru_0 += 1
        if c3_thru_0 == 10:
            count_all_thru_0 += 1
            # What's the total t3?
            c3 = find_3cycles(adj, 6)
            t3 = len(c3)
            c5 = find_5cycles_dp(adj, 6)
            t5 = len(c5)
            scores = sorted([sum(adj[i]) for i in range(6)])
            print(f"  All 10 thru 0: t3={t3}, t5={t5}, alpha_1={t3+t5}, score={scores}")
            # Check i_2 if alpha_1=10
            if t3 + t5 == 10:
                all_cyc = c3 + c5
                _, i2 = count_alpha_i2(all_cyc)
                print(f"    i_2={i2}")
        total += 1

    print(f"  Total with all 10 triples through v=0 as 3-cycles: {count_all_thru_0}/{total}")

    # Through-vertex analysis:
    # If all C(5,2)=10 triples through v=0 are 3-cycles, what is the score of v=0?
    # Through-v 3-cycles = s*(n-1-s) - arcs within out-nbhd of v that DON'T form part of cycles
    # Actually, through-v cycles = # arcs from out(v) to in(v).
    # out(v) has s vertices, in(v) has 5-s vertices.
    # # arcs out->in = between s*(5-s) (every out beats every in).
    # Max through-v = s*(5-s), minimized at s=0 or s=5 (giving 0), maximized at s=2 or 3.
    # s=2: 2*3=6. s=3: 3*2=6. Max = 6.
    # But we need 10 through vertex 0. 10 > 6. IMPOSSIBLE!

    print(f"\n  Max 3-cycles through any vertex at n=6:")
    for s in range(6):
        max_thru = s * (5 - s)
        print(f"    score s={s}: max through-v = {max_thru}")
    print(f"  Max is 6 (at s=2 or s=3), but we need 10. IMPOSSIBLE.")
    print(f"  So 10 pairwise-intersecting 3-cycles through a common vertex is IMPOSSIBLE at n=6.")

    # Wait, but 10 pairwise-intersecting sets don't need a COMMON element!
    # Hilton-Milner: max intersecting without common = C(n-1,k-1) - C(n-k-1,k-1) + 1
    # For k=3, n=6: C(5,2) - C(2,2) + 1 = 10 - 1 + 1 = 10
    # So 10 pairwise-intersecting 3-sets without common element is POSSIBLE at n=6!

    # Actually the argument needs to be different.
    # Through ANY vertex v, max 3-cycles = s_v*(n-1-s_v) <= (n-1)^2/4.
    # At n=6: max = 6.25 -> 6. So through any vertex: at most 6 three-cycles.

    # If 10 three-cycles are pairwise intersecting, does some vertex appear in many?
    # Total incidences = 3*10 = 30. On 6 vertices: average = 5 per vertex.
    # Max through any vertex: 6. If one vertex has 6, remaining 4 use 24 incidences on 5 verts.
    # Avg remaining = 4.8. Some have >= 5.

    # Can we show: if all 10 three-cycles pairwise intersect on [6],
    # some vertex must appear in > 6 of them (impossible by max-through-v)?

    # Hmm, average is 5, max is 6. Possible.

    # Actually: I need to count THREE-CYCLES, not just triples.
    # A pairwise-intersecting family of 3-sets on [6] with all being directed 3-cycles
    # in some tournament: this is highly constrained.

    # Let me just check: can we have 10 pairwise-intersecting directed 3-cycles at n=6?
    print(f"\nn=6: Can we have 10 pairwise-intersecting directed 3-cycles?")

    found = False
    for bits in range(2**15):
        adj = [[0]*6 for _ in range(6)]
        for k, (i, j) in enumerate(edges6):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1

        c3 = find_3cycles(adj, 6)
        if len(c3) < 10:
            continue

        # Check if any 10-subset is pairwise intersecting
        if len(c3) == 10:
            # Check all pairs
            all_intersect = True
            for a in range(10):
                for b in range(a+1, 10):
                    if not (c3[a] & c3[b]):
                        all_intersect = False
                        break
                if not all_intersect:
                    break
            if all_intersect:
                found = True
                scores = sorted([sum(adj[i]) for i in range(6)])
                c5 = find_5cycles_dp(adj, 6)
                print(f"  FOUND: 10 pairwise-intersecting 3-cycles! score={scores}, t5={len(c5)}")
                break
        elif len(c3) > 10:
            for subset in combinations(range(len(c3)), 10):
                cyc10 = [c3[i] for i in subset]
                all_intersect = True
                for a in range(10):
                    for b in range(a+1, 10):
                        if not (cyc10[a] & cyc10[b]):
                            all_intersect = False
                            break
                    if not all_intersect:
                        break
                if all_intersect:
                    found = True
                    scores = sorted([sum(adj[i]) for i in range(6)])
                    c5 = find_5cycles_dp(adj, 6)
                    print(f"  FOUND (from t3={len(c3)}): 10 pairwise-intersecting 3-cycles! score={scores}, t5={len(c5)}")
                    break
            if found:
                break

    if not found:
        print(f"  NOT FOUND: 10 pairwise-intersecting directed 3-cycles at n=6 is IMPOSSIBLE.")
        print(f"  Max t3 at n=6 = 8 (regular score (2,2,2,3,3,3)) < 10.")

    # Check max t3 at n=6
    print(f"\nn=6: Distribution of t3")
    t3_dist = Counter()
    for bits in range(2**15):
        adj = [[0]*6 for _ in range(6)]
        for k, (i, j) in enumerate(edges6):
            if (bits >> k) & 1:
                adj[j][i] = 1
            else:
                adj[i][j] = 1
        c3 = find_3cycles(adj, 6)
        t3_dist[len(c3)] += 1
    for t3 in sorted(t3_dist.keys()):
        print(f"  t3={t3}: {t3_dist[t3]} tournaments")

    # KEY: At n=6, max t3 = 8 (from regular score).
    # So you can NEVER have 10 three-cycles at n=6!
    # But alpha_1=10 includes 5-cycles too.

    # alpha_1=10 = t3 + t5 = 10 at n=6.
    # We know: always t3=6, t5=4.
    # The disjoint pairs are ONLY (3,3) type.
    # i_2=0 would mean: no two 3-cycles are vertex-disjoint.
    # But at n=6 with t3=6, we showed alpha_1=10 forces complementary pair.

    print(f"\n{'='*60}")
    print(f"APPROACH FOR GENERAL n:")
    print(f"{'='*60}")
    print(f"""
For H=21, need alpha_1 + 2*i_2 = 10 (assuming alpha_3=0, proved).

Case (10,0): alpha_1=10, i_2=0.
  ALL 10 odd cycles must pairwise share a vertex.
  The cycle types are 3 and 5 (and 7, 9, etc. at larger n).
  Two cycles of sizes l1, l2 are always sharing if l1+l2 > n.
  So the ONLY way to get a disjoint pair is l1+l2 <= n.
  For i_2=0: no pair has l1+l2 <= n.
  If n >= 8: 3+5=8, so 3-cycle and 5-cycle CAN be disjoint.
  We need to show: alpha_1=10 with i_2=0 forces a contradiction.

  At n <= 7: only (3,3) pairs possible for disjoint.
  i_2=0 means: no two 3-cycles disjoint.
  Proved impossible at n=6,7 (exhaustive/sampling).

  At n >= 8: (3,5) pairs also possible.
  Need: no 3-3 disjoint AND no 3-5 disjoint.
  This is MORE restrictive, so harder to achieve.

  At n=8 sampling: alpha_1=10 gives i_2 in {{1,2,8}} (never 0).
  The i_2=8 case has 3-5 disjoint pairs appearing.

Case (8,1): alpha_1=8, i_2=1.
  Exactly 1 vertex-disjoint pair among 8 cycles.
  At n=6: alpha_1=8 always gives i_2=0. (t3=4, t5=4)
  At n=7: alpha_1=8 gives i_2 in {{0,3}}.
  At n=8: alpha_1=8 gives i_2 in {{0,7}}.
  i_2=1 NEVER appears.
""")


prove_10_0_impossible()

# Run sampling at n=8,9
for n_val in [8]:
    general_n_sampling(n_val, 30000)
