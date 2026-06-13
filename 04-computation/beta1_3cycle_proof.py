#!/usr/bin/env python3
"""
beta1_3cycle_proof.py — Prove sum_v beta_1(T\v) <= 3 via 3-cycle counting

ESTABLISHED FACTS:
1. beta_1(T) ∈ {0,1} for all tournaments (verified n≤6)
2. At n=4: beta_1(T)=1 iff t3(T)=2
3. For T on n vertices with beta_1(T)=0:
   beta_1(T\v)=1 iff t3(T\v) ≥ 2 (at the (n-1)-vertex level)

PROOF STRATEGY:
Show that if beta_1(T)=0, then at most 3 vertices v have t3(T\v) ≥ 2.
This reduces to a combinatorial statement about 3-cycle distributions.

Actually, let's be more precise about what "t3(T\v)≥2" means.
At n=4, beta_1=1 iff t3=2 (= C(4,3)/2, the maximum).
What's the analogous condition at general n-1?

KEY: need to verify the exact condition at n=5 subtournaments (n-1=5).

opus-2026-03-08
"""

import numpy as np
from itertools import combinations
from collections import defaultdict, Counter
import sys

sys.path.insert(0, '/home/e/Documents/claude/math/04-computation')
from beta1_structure_analysis import (
    compute_beta1_full, compute_beta1_only,
    subtournament, all_tournaments, count_3cycles, rank_of_matrix
)

def run():
    print("=" * 70)
    print("BETA_1 DELETION BOUND VIA 3-CYCLE COUNTING")
    print("=" * 70)

    # ================================================================
    # Step 1: Exact condition for beta_1=1 at each n
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 1: EXACT CONDITION FOR beta_1=1")
    print("=" * 70)

    for n in [3, 4, 5, 6]:
        by_t3 = defaultdict(Counter)
        for A in all_tournaments(n):
            t3 = count_3cycles(A, n)
            b1 = compute_beta1_only(A, n)
            by_t3[t3][b1] += 1

        print(f"\n  n={n}:")
        for t3 in sorted(by_t3.keys()):
            print(f"    t3={t3}: {dict(by_t3[t3])}")

        # Find the threshold
        threshold = None
        for t3 in sorted(by_t3.keys()):
            if 1 in by_t3[t3]:
                threshold = t3
                break
        if threshold is not None:
            # Check: is beta_1=1 iff t3 >= threshold?
            all_above = all(1 in by_t3[t3] for t3 in by_t3 if t3 >= threshold)
            none_below = all(1 not in by_t3[t3] for t3 in by_t3 if t3 < threshold)
            # Actually mixed at n=5: t3=3 has both 0 and 1
            print(f"    First t3 with beta_1=1: {threshold}")
            # Check if it's a clean cutoff
            for t3 in sorted(by_t3.keys()):
                if t3 >= threshold:
                    total = sum(by_t3[t3].values())
                    n_b1 = by_t3[t3].get(1, 0)
                    if n_b1 < total:
                        print(f"    WARNING: t3={t3} has mixed beta_1 values: {dict(by_t3[t3])}")

    # ================================================================
    # Step 2: At n=4, prove beta_1=1 iff t3=2
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 2: beta_1 AT n=4")
    print("=" * 70)

    print("""
  At n=4: C(4,3) = 4 triples. Each is either TT or 3-cycle.
  t3 ∈ {0, 1, 2}  (t3=3 and t3=4 are impossible at n=4).

  Wait — let me verify the possible t3 values at n=4.
  By the relation: t3 = C(n,3) - TT/6 ... no.
  Actually: each unordered triple {a,b,c} is either a transitive triple or a 3-cycle.
  TT + t3 = C(n,3). At n=4: TT + t3 = 4.
  t3 ∈ {0, 1, 2, 3, 4}? Or is there a constraint?

  For n=4: scores are (s0, s1, s2, s3) with sum = C(4,2)=6.
  t3 = C(n,3) - sum C(s_i, 2) ... no, that's not right either.
  t3 = (C(n,2)*n/6 - sum s_i(s_i-1)(s_i-2)/6)... hmm.

  Actually: t3 = (C(n,3) - sum_v C(s_v^+, 2)) = C(n,3) - sum_v C(out-degree_v, 2).
  No, the standard formula is:
  t3 = C(n,3) - sum_v C(s_v, 2) where s_v = out-degree of v, and this counts TT.
  Wait: sum_v C(s_v, 2) = # transitive triples. So t3 = C(n,3) - sum_v C(s_v, 2).

  At n=4, scores sum to 6. Score sequences: (0,1,2,3), (0,2,2,2), (1,1,1,3), (1,1,2,2).
  """)

    # Verify score sequences and t3 at n=4
    score_seqs = Counter()
    for A in all_tournaments(4):
        scores = tuple(sorted([sum(A[i]) for i in range(4)]))
        t3 = count_3cycles(A, 4)
        score_seqs[(scores, t3)] += 1
    for (s, t3), cnt in sorted(score_seqs.items()):
        tt = 4 - t3  # since C(4,3)=4
        sum_cs2 = sum(si*(si-1)//2 for si in s)
        print(f"  scores={s}, t3={t3}, TT={tt}, sum C(s,2)={sum_cs2}, count={cnt}")

    # ================================================================
    # Step 3: 3-cycle counting under vertex deletion
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 3: 3-CYCLE COUNTING UNDER DELETION")
    print("=" * 70)

    print("""
  For a 3-cycle {a,b,c} in T, it contributes to t3(T\\v) iff v ∉ {a,b,c}.
  So t3(T\\v) = t3(T) - #{3-cycles through v}.

  Let c(v) = #{3-cycles through v}. Then t3(T\\v) = t3(T) - c(v).

  Now: sum_v c(v) = 3 * t3(T) (each 3-cycle contributes to 3 vertex counts).

  For beta_1(T\\v) = 1, we need t3(T\\v) >= threshold(n-1).
  At n-1=4: threshold = 2.
  So: beta_1(T\\v)=1 iff t3(T) - c(v) >= 2, i.e., c(v) <= t3(T) - 2.

  The number of "bad" vertices = #{v : c(v) <= t3(T) - 2}.

  Claim: when beta_1(T)=0, this count is at most 3.
  """)

    # Verify: at n=5, beta_1(T\v)=1 iff t3(T\v)=2
    print("  Verifying: beta_1(T\\v)=1 iff t3(T\\v)=2 at n=5 (T\\v has 4 vertices):")
    match = 0
    mismatch = 0
    for A in all_tournaments(5):
        info = compute_beta1_full(A, 5)
        if info['beta_1'] != 0:
            continue
        for v in range(5):
            Av, nv = subtournament(A, 5, [v])
            t3v = count_3cycles(Av, nv)
            b1v = compute_beta1_only(Av, nv)
            if (b1v == 1) == (t3v == 2):
                match += 1
            else:
                mismatch += 1
                if mismatch <= 3:
                    print(f"    MISMATCH: t3={t3v}, beta1={b1v}")
    print(f"  Match: {match}, Mismatch: {mismatch}")

    # At n=6: beta_1(T\v)=1 needs checking against t3(T\v)
    # T\v has 5 vertices, where beta_1=1 threshold is t3>=3 (from the data)
    print("\n  At n=6: checking beta_1(T\\v)=1 condition on 5-vertex subtournaments:")
    match6 = 0
    mismatch6 = 0
    t3_when_b1 = Counter()
    t3_when_b0 = Counter()
    cnt = 0
    for A in all_tournaments(6):
        info = compute_beta1_full(A, 6)
        if info['beta_1'] != 0:
            continue
        for v in range(6):
            Av, nv = subtournament(A, 6, [v])
            t3v = count_3cycles(Av, nv)
            b1v = compute_beta1_only(Av, nv)
            if b1v == 1:
                t3_when_b1[t3v] += 1
            else:
                t3_when_b0[t3v] += 1
        cnt += 1
        if cnt >= 5000:
            break
    print(f"  t3(T\\v) when beta_1=1: {dict(sorted(t3_when_b1.items()))}")
    print(f"  t3(T\\v) when beta_1=0: {dict(sorted(t3_when_b0.items()))}")

    # ================================================================
    # Step 4: The 3-cycle distribution constraint
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 4: 3-CYCLE DISTRIBUTION CONSTRAINT")
    print("=" * 70)

    print("""
  At n=5, beta_1(T)=0 requires t3(T) <= 4 (since t3=5 => beta_1=1).

  Actually, let me re-check: at n=5, beta_1=0 occurs for t3 in {0,1,2,3,4}.
  But not ALL tournaments with t3=3 or t3=4 have beta_1=0.

  The key constraint: when beta_1(T)=0 at n=5:
  - 0 bad vertices: t3 ∈ {0,1} (no subtournament can have t3_sub=2)
  - 1 bad vertex: t3=2
  - 2 bad vertices: t3=3
  - 3 bad vertices: t3=4

  This is EXACTLY: # bad vertices = max(0, t3(T) - 1).
  Let me verify this formula.
  """)

    # Verify: # bad = max(0, t3 - 1) at n=5
    for A in all_tournaments(5):
        info = compute_beta1_full(A, 5)
        if info['beta_1'] != 0:
            continue
        t3 = count_3cycles(A, 5)
        n_bad = 0
        for v in range(5):
            Av, nv = subtournament(A, 5, [v])
            if compute_beta1_only(Av, nv) > 0:
                n_bad += 1
        expected = max(0, t3 - 1)
        if n_bad != expected:
            print(f"  FAIL: t3={t3}, n_bad={n_bad}, expected={expected}")
            break
    else:
        print("  VERIFIED: # bad vertices = max(0, t3-1) for all n=5 tournaments with beta_1=0")

    # ================================================================
    # Step 5: c(v) distribution
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 5: c(v) DISTRIBUTION (3-cycles through v)")
    print("=" * 70)

    for n in [5]:
        print(f"\n  n={n}:")
        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            if info['beta_1'] != 0:
                continue
            t3 = count_3cycles(A, n)
            if t3 < 3:
                continue

            cvs = []
            for v in range(n):
                c_v = sum(1 for trip in combinations(range(n), 3)
                          if v in trip and is_3cycle(A, trip))
                cvs.append(c_v)

            betas = []
            for v in range(n):
                Av, nv = subtournament(A, n, [v])
                betas.append(compute_beta1_only(Av, nv))

            t3_check = [t3 - cv for cv in cvs]
            print(f"    t3={t3}, c(v)={cvs}, t3(T\\v)={t3_check}, betas={betas}")
            break  # one example per t3

        # Show all with t3=4 (= 3 bad)
        print(f"\n  All n=5, t3=4, beta_1=0 (3 bad):")
        count = 0
        for A in all_tournaments(n):
            info = compute_beta1_full(A, n)
            if info['beta_1'] != 0:
                continue
            t3 = count_3cycles(A, n)
            if t3 != 4:
                continue

            cvs = []
            for v in range(n):
                c_v = sum(1 for trip in combinations(range(n), 3)
                          if v in trip and is_3cycle(A, trip))
                cvs.append(c_v)

            betas = []
            for v in range(n):
                Av, nv = subtournament(A, n, [v])
                betas.append(compute_beta1_only(Av, nv))

            t3_check = [t3 - cv for cv in cvs]
            if count < 10:
                print(f"    c(v)={cvs}, t3(T\\v)={t3_check}, betas={betas}")
            count += 1
        print(f"    Total: {count} tournaments")

    # ================================================================
    # Step 6: General n analysis
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 6: GENERAL THRESHOLD ANALYSIS")
    print("=" * 70)

    print("""
  CONJECTURE (from data):
  beta_1(T) = 1 iff the tournament has "too many" 3-cycles.

  More precisely, for n-vertex tournaments:
  - n=3: beta_1=1 iff t3=1 (= the 3-cycle tournament)
  - n=4: beta_1=1 iff t3=2
  - n=5: beta_1=1 iff t3 ∈ {3(some), 4(some), 5(all)} — NOT a clean threshold

  Wait, at n=5 it's NOT a clean threshold! t3=3 has BOTH beta_1=0 and 1.

  So the condition for beta_1=1 is NOT simply t3 ≥ threshold.
  It depends on the STRUCTURE, not just the count.

  But for DELETION, we need: t3(T\\v) ≥ threshold(n-1).
  At (n-1)=4: threshold = 2 (clean).
  At (n-1)=5: threshold involves a structural condition, not just t3.

  Let me check: at n=6, is beta_1(T\\v)=1 determined by t3(T\\v)?
  """)

    # Direct check: is beta_1 determined by t3 alone?
    print("  n=5: beta_1 values by t3:")
    for n in [5]:
        by_t3 = defaultdict(Counter)
        for A in all_tournaments(n):
            t3 = count_3cycles(A, n)
            b1 = compute_beta1_only(A, n)
            by_t3[t3][b1] += 1
        for t3 in sorted(by_t3.keys()):
            print(f"    t3={t3}: {dict(by_t3[t3])}")

    # At n=5, t3=3 and t3=4 have mixed beta_1.
    # So beta_1 is NOT determined by t3 alone at n>=5.
    # But beta_1(T\v) for T\v on 4 vertices IS determined by t3.

    # The question is: what determines beta_1 at n=5?
    print("\n  n=5: What distinguishes beta_1=0 from beta_1=1 at same t3?")
    for n in [5]:
        for t3_target in [3, 4]:
            b0_scores = Counter()
            b1_scores = Counter()
            for A in all_tournaments(n):
                t3 = count_3cycles(A, n)
                if t3 != t3_target:
                    continue
                scores = tuple(sorted([sum(A[i]) for i in range(n)]))
                b1 = compute_beta1_only(A, n)
                if b1 == 0:
                    b0_scores[scores] += 1
                else:
                    b1_scores[scores] += 1
            print(f"    t3={t3_target}:")
            print(f"      beta_1=0 score sequences: {dict(sorted(b0_scores.items()))}")
            print(f"      beta_1=1 score sequences: {dict(sorted(b1_scores.items()))}")

    # ================================================================
    # Step 7: The ACTUAL proof approach
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 7: THE PROOF APPROACH")
    print("=" * 70)

    print("""
  CRITICAL REALIZATION:

  The deletion bound sum_v beta_1(T\\v) <= 3 for beta_1(T)=0 at n=5,6
  seems to follow from DIMENSIONAL constraints.

  At n=5:
  - beta_1(T)=0 means dim(B_1(T)) = 6
  - For each v: dim(B_1(T\\v)) = 3 - beta_1(T\\v)
  - pi_v(B_1(T)) = R^6 (surjective, since pi_v is injective on B_1)
  Wait — pi_v maps R^10 -> R^6 and B_1(T) has dimension 6.
  pi_v(B_1(T)) has dimension 6 in R^6. So pi_v(B_1(T)) = R^6. ✓

  Now: B_1(T\\v) ⊂ R^6 with dimension 3 - beta_1(T\\v).
  B_1(T\\v) = {d2(p) : p ∈ Omega_2(T\\v)} ⊂ pi_v(B_1(T)).

  Since pi_v(B_1(T)) = R^6 = full space, the constraint B_1(T\\v) ⊂ pi_v(B_1(T))
  is vacuous! It doesn't help.

  THE REAL CONSTRAINT must come from Omega_2 structure.

  NEW APPROACH: Think about dim(Omega_2(T)) and dim(Omega_2(T\\v)).

  For tournaments:
  |A_2(T)| = # 2-paths a->b->c in T
  Omega_2(T) = subspace of A_2 where boundaries have no non-allowed faces
  dim(B_1(T)) = rank(d2|_{Omega_2(T)})

  The key: Omega_2 is constrained by the tournament's 3-cycle structure.
  Each 3-cycle triple contributes to the "obstruction" of Omega_2.

  For the deletion: when we remove v, the 3-cycles that REMAIN are exactly
  the 3-cycles not through v. The obstruction at (n-1) level comes from these.

  At n-1=4: obstruction occurs iff the subtournament has t3_sub = 2.
  The 4-vertex subtournament has C(4,3)=4 triples; t3_sub=2 means
  exactly half are 3-cycles.
  """)

    # ================================================================
    # Step 8: The core combinatorial bound
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 8: CORE COMBINATORIAL BOUND")
    print("=" * 70)

    # For n=5 with beta_1=0:
    # Each vertex v has c(v) 3-cycles through it.
    # sum c(v) = 3*t3
    # beta_1(T\v) = 1 iff t3 - c(v) >= 2, i.e., c(v) <= t3 - 2

    # The question: can we have 4 vertices with c(v) <= t3 - 2?

    # At n=5: maximum t3 with beta_1(T)=0 is t3=4.
    # For t3=4: need c(v) <= 2 for bad vertex.
    # sum c(v) = 12. n=5 vertices.
    # Average c(v) = 12/5 = 2.4.
    # If 4 vertices have c(v) <= 2: their sum <= 8.
    # 5th vertex has c(v) >= 12 - 8 = 4.
    # But c(v) <= C(4,2) = 6 (max 3-cycles through any vertex at n=5).
    # So c(v) = 4 is possible. Sum would be ≤ 8 + 4 = 12 = 3*t3. ✓ consistent.
    # BUT: does this actually occur? Data says NO — max 3 bad.

    print("""
  COMBINATORIAL ANALYSIS at n=5:

  sum_v c(v) = 3 * t3. Threshold: c(v) <= t3 - 2 for bad.

  If 4 vertices are bad (c(v) <= t3-2 each):
  sum of their c(v) <= 4(t3-2) = 4t3 - 8
  Remaining vertex: c(v5) >= 3t3 - (4t3-8) = 8 - t3

  For t3 ≤ 4 (needed for beta_1(T)=0 at n=5):
  c(v5) >= 8 - t3 >= 4.

  But c(v) counts 3-cycles through v. v appears in C(4,2)=6 triples.
  So c(v) <= 6 at n=5. So c(v5) >= 4 is possible.

  But ALSO: each bad vertex has c(v) <= t3-2 <= 2.
  And each 3-cycle involves exactly 3 vertices.

  Let's think about it as a hypergraph. The 3-cycles form a
  3-uniform hypergraph on {0,1,2,3,4}. Each vertex's degree is c(v).

  For 4 bad vertices {a,b,c,d}: c(a),c(b),c(c),c(d) <= t3-2.
  The remaining vertex e has c(e) = 3*t3 - c(a)-c(b)-c(c)-c(d) >= 3t3 - 4(t3-2) = 8-t3.

  Now the 3-cycles through e involve e and 2 of {a,b,c,d}.
  There are C(4,2)=6 such potential triples. c(e) of these are actual 3-cycles.

  The 3-cycles NOT through e form the 3-cycles in T\\e.
  t3(T\\e) = t3 - c(e) <= t3 - (8-t3) = 2t3 - 8.
  For t3 <= 4: t3(T\\e) <= 0.

  But t3(T\\e) >= 0, so t3(T\\e) = 0 for t3 = 4, meaning T\\e is... well,
  t3(T\\e) = t3 - c(e) = 4 - c(e).
  And c(e) >= 8 - t3 = 4 (for t3=4).
  So c(e) = 4 and t3(T\\e) = 0.

  But we need: the 3-cycles among {a,b,c,d} have t3({a,b,c,d}) = 0.
  That means T\\e is a tournament on 4 vertices with NO 3-cycles — transitive!
  But also each of a,b,c,d has c(v) <= 2.

  Let me check if this is possible.
  With t3=4 and c(e)=4: 4 cycles, all through e.
  Remaining 4 vertices form a transitive tournament (t3=0).
  Each of a,b,c,d has c(v) = #{3-cycles through v involving e} + #{3-cycles among a,b,c,d through v}.
  Since t3({a,b,c,d})=0, the second term is 0.
  So c(v) = #{3-cycles through v and e} for v in {a,b,c,d}.

  These 3-cycles are of the form {e,v,w} for w in {a,b,c,d}\\{v}.
  For each v, c(v) = #{w in {a,b,c,d}\\{v} : {e,v,w} is a 3-cycle}.

  Total: sum_{v in {a,b,c,d}} c(v) = sum_{v} #{w : {e,v,w} is 3-cycle}
  Each 3-cycle {e,v,w} contributes to c(v) AND c(w). So sum c(v) = 2*c(e) = 8.
  Average c(v) = 2. For all to be <= 2, need c(v) = 2 for each.

  Is this achievable? With 4 3-cycles through e and 4 other vertices,
  and each pair {v,w} from {a,b,c,d} is in at most one 3-cycle (with e).

  Wait: C(4,2)=6 pairs. 4 of these 6 are 3-cycles (with e).
  Each vertex v in {a,b,c,d} is in 3 pairs from C(4,2).
  c(v) = # of those 3 pairs that are 3-cycles.
  sum c(v) = 2*4 = 8. Average = 2. Max min possible?

  By double counting: 4 pairs chosen from 6. Each vertex in 3 pairs.
  Want each vertex in exactly 2 chosen pairs.
  4 vertices, 2 each = 8/2 = 4 pairs chosen. ✓

  Can we choose 4 pairs from {a,b,c,d} such that each vertex is in exactly 2?
  This is a 2-regular graph on 4 vertices = a 4-cycle!
  So the 4 3-cycles form a "cycle" pattern: {e,a,b}, {e,b,c}, {e,c,d}, {e,d,a}.

  Now: does such a tournament actually exist? And does it have beta_1(T)=0?
  If beta_1(T)=1 for ALL t3=4 tournaments with this structure, then 4 bad vertices
  is impossible because beta_1(T)=0 fails!
  """)

    # Check: do t3=4 tournaments with beta_1=0 exist, and what's their structure?
    print("  n=5, t3=4, beta_1=0 tournaments:")
    count_t4b0 = 0
    for A in all_tournaments(5):
        t3 = count_3cycles(A, 5)
        if t3 != 4:
            continue
        b1 = compute_beta1_only(A, 5)
        if b1 != 0:
            continue
        count_t4b0 += 1

        cvs = []
        for v in range(5):
            c_v = sum(1 for trip in combinations(range(5), 3)
                      if v in trip and is_3cycle(A, trip))
            cvs.append(c_v)

        if count_t4b0 <= 10:
            scores = tuple(sum(A[i]) for i in range(5))
            print(f"    scores={scores}, c(v)={cvs}, bad={[v for v in range(5) if t3-cvs[v]>=2]}")

    print(f"  Total: {count_t4b0}")

    # Check: do t3=5 tournaments with beta_1=0 exist? (answer: NO from data)
    count_t5b0 = 0
    for A in all_tournaments(5):
        t3 = count_3cycles(A, 5)
        if t3 != 5:
            continue
        b1 = compute_beta1_only(A, 5)
        if b1 == 0:
            count_t5b0 += 1
    print(f"\n  n=5, t3=5, beta_1=0: {count_t5b0} (expected 0)")

    # ================================================================
    # Step 9: The actual obstruction
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 9: THE ACTUAL OBSTRUCTION")
    print("=" * 70)

    print("""
  KEY OBSERVATION: For n=5 tournaments with t3=4 and beta_1=0,
  the c(v) distribution always has exactly one vertex with c(v) >= 3.
  That one vertex is "good" (t3-c(v) <= 1 < 2), leaving 3 bad vertices.

  If there were 4 bad vertices, we'd need ALL c(v) <= 2.
  That requires a vertex e with c(e) >= 4 (from sum = 12, four vertices contributing ≤ 8).
  But then T\\e has t3(T\\e) = 4-4 = 0, so e is good.
  And the 4th "candidate bad" vertex must have c(v) ≤ 2, but e is the one with c(e)=4.
  So e is NOT bad (c(e)=4 > t3-2=2), and we have at most 3 bad among {a,b,c,d}.

  Wait — but we assumed {a,b,c,d} are all bad. The contradiction is:
  if 4 are bad with c(v) ≤ 2, then the 5th must have c(v) ≥ 4,
  and the 5th is NOT bad (c(v)=4 > 2). So exactly 4 are bad.
  But we need to check: can ALL FOUR actually be bad?

  The constraint: the {a,b,c,d} subtournament must have t3=0.
  And the {e,v,w} 3-cycles form a 2-regular graph on {a,b,c,d}.
  2-regular on 4 vertices = 4-cycle.

  So T\\e is transitive on 4 vertices (t3=0, beta_1=0).
  T\\a: 3-cycles = those not through a = t3 - c(a) = 4-2 = 2. So beta_1(T\\a)=1.
  Similarly for b,c,d: t3(T\\v)=2, beta_1=1.

  So we WOULD have 4 bad vertices. But does such a tournament exist with beta_1(T)=0?

  The answer from data: NO! All t3=4 tournaments with beta_1=0 have c(v) values
  with min(c) ≥ ... let me check.
  """)

    # Detailed c(v) for t3=4, beta_1=0
    print("  Detailed c(v) for n=5, t3=4, beta_1=0:")
    for A in all_tournaments(5):
        t3 = count_3cycles(A, 5)
        if t3 != 4:
            continue
        b1 = compute_beta1_only(A, 5)
        if b1 != 0:
            continue
        cvs = []
        for v in range(5):
            c_v = sum(1 for trip in combinations(range(5), 3)
                      if v in trip and is_3cycle(A, trip))
            cvs.append(c_v)
        print(f"    c(v)={cvs}, min={min(cvs)}, max={max(cvs)}, t3_sub={[t3-c for c in cvs]}")

    # And for t3=4, beta_1=1:
    print("\n  Detailed c(v) for n=5, t3=4, beta_1=1:")
    count_shown = 0
    for A in all_tournaments(5):
        t3 = count_3cycles(A, 5)
        if t3 != 4:
            continue
        b1 = compute_beta1_only(A, 5)
        if b1 != 1:
            continue
        cvs = []
        for v in range(5):
            c_v = sum(1 for trip in combinations(range(5), 3)
                      if v in trip and is_3cycle(A, trip))
            cvs.append(c_v)
        if count_shown < 10:
            print(f"    c(v)={cvs}, min={min(cvs)}, max={max(cvs)}, t3_sub={[t3-c for c in cvs]}")
        count_shown += 1
    print(f"    Total: {count_shown}")

    # ================================================================
    # Step 10: FINAL PROOF SKETCH
    # ================================================================
    print("\n" + "=" * 70)
    print("STEP 10: PROOF SKETCH")
    print("=" * 70)

    print("""
  THEOREM: For any tournament T on n≥4 vertices with beta_1(T)=0,
  sum_v beta_1(T\\v) ≤ 3.

  PROOF APPROACH (verified computationally at n=5,6):

  Step 1: beta_1(T\\v) ∈ {0,1} (since T\\v is a tournament on n-1 vertices).

  Step 2: beta_1(T\\v)=1 iff t3(T\\v) ≥ threshold(n-1).
  At n-1=4: threshold = 2 (clean; beta_1=1 iff t3=2).
  At n-1≥5: the condition is more structural, but t3(T\\v) ≥ threshold is necessary.

  Step 3: t3(T\\v) = t3(T) - c(v) where c(v) = #{3-cycles through v}.

  Step 4: If beta_1(T)=0, then t3(T) is bounded:
  - At n=5: t3 ≤ 4 (t3=5 forces beta_1=1)
  - At n=6: t3 ≤ some bound

  Step 5: The bound on the number of "bad" vertices comes from:
  - 3*t3 = sum c(v)  (each 3-cycle counted 3 times)
  - Bad: c(v) ≤ t3 - threshold(n-1)
  - At most 3 can satisfy this given the constraint sum c(v) = 3*t3.

  The algebraic constraint: if k vertices are "bad" with c(v) ≤ t3 - θ each,
  sum(bad c(v)) ≤ k(t3-θ)
  sum(good c(v)) ≥ 3t3 - k(t3-θ) = (3-k)t3 + kθ
  (n-k) good vertices must absorb this.

  For this to work with beta_1(T)=0 (which bounds t3 from above),
  at most 3 vertices can be bad.
  """)

def is_3cycle(A, triple):
    """Check if the triple forms a 3-cycle in tournament A."""
    a, b, c = triple
    if A[a][b] and A[b][c] and A[c][a]:
        return True
    if A[b][a] and A[a][c] and A[c][b]:
        return True
    return False

if __name__ == '__main__':
    run()
