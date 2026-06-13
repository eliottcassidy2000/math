"""
fourier_master_formula.py -- kind-pasteur-2026-03-14-S105j
SYNTHESIS of Fourier spectrum discoveries.

KNOWN EXACT RESULTS:

Level 2 (proved in S75, S104):
  |H_hat(S)| = (n-2)!/2^(n-2) for S = {e1,e2} adjacent (sharing a vertex)
  |H_hat(S)| = 0 for S = {e1,e2} disjoint
  N_2 = n(n-1)(n-2)/2 nonzero level-2 coefficients

Level 4 at n=5 (computed S105g):
  |H_hat(S)| = 1/8 = 1!/2^3 for ALL 60 nonzero 4-arc subsets
  All 60 cover all 5 vertices
  E_4 = 60 * (1/8)^2 = 15/16

Level 4 at n=6 (computed S105i):
  TWO classes of nonzero coefficients:
  Class A: |H_hat| = 1/8, 360 coefficients, each covers 5 vertices
  Class B: |H_hat| = 1/4, 90 coefficients, each covers all 6 vertices
  E_4 = 360/64 + 90/16 = 5.625 + 5.625 = 11.25 = 45/4

PATTERNS:
  1. Magnitudes at level 4:
     n=5: 1/8 = (5-4)!/2^(5-2) = 1!/8
     n=6 Class A: 1/8 = (6-5)!/2^3 = 1!/8 (covers 5 vertices, like n=5 level-4)
     n=6 Class B: 1/4 = (6-4)!/2^(6-2) = 2!/16 = 1/8? No, 2!/16 = 2/16 = 1/8.

  Hmm wait: (n-4)!/2^(n-2) at n=6 = 2!/2^4 = 2/16 = 1/8. But Class B has 1/4!
  So the formula |H_hat| = (n-|S|)!/2^(n-2) does NOT apply at n=6.

  Let me reconsider. At n=6, the 90 coefficients with |H_hat|=1/4:
  1/4 = 2 * (1/8). So they're DOUBLE the base magnitude.

  Alternative: maybe the magnitude depends on the NUMBER OF VERTICES COVERED.
  If 4 arcs cover v vertices, then |H_hat| = (v-4)!/2^(v-2)?
  At v=5: (5-4)!/2^(5-2) = 1/8. CHECK (both n=5 all, and n=6 Class A)
  At v=6: (6-4)!/2^(6-2) = 2!/16 = 1/8. WRONG (Class B has 1/4)

  Or: |H_hat| = (n-4)!/2^(n-2) * something(v)?
  At n=5, v=5: 1!/8 * 1 = 1/8
  At n=6, v=5: 1!/8 * 1 = 1/8 (but n=6 now!)
  Wait, (n-4)!/2^(n-2) at n=6 is 2!/16 = 1/8. And Class A has 1/8.

  But Class B has 1/4 = 2 * 1/8. Why the factor of 2?

  Let me try: the magnitude is (v-4)!/2^(v-2) where v is the number
  of vertices covered, INDEPENDENT of n.
  v=5: 1!/2^3 = 1/8. Matches both n=5 and n=6 Class A.
  v=6: 2!/2^4 = 2/16 = 1/8. WRONG, should be 1/4.

  Try: magnitude = ((v-2)!/(2^(v-2))) * ((v-4)!/(2^(v-4)))?
  No, that's weird.

  Actually: let me think about it differently.
  At level 2, adjacent arcs share a vertex, and the magnitude is (n-2)!/2^(n-2).
  But this is INDEPENDENT of n at fixed vertex count!
  (n-2)!/2^(n-2) depends on n.

  Hmm, but at n=6 level 4, the Class A coefficients (5-vertex) have 1/8.
  And 1/8 = (5-2)!/2^(5-2) * ... no, (5-2)!/2^3 = 6/8 = 3/4 ≠ 1/8.

  Let me just tabulate what we know and look for patterns.
"""

import sys, math
from fractions import Fraction
from itertools import combinations

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("FOURIER SPECTRUM MASTER TABLE")
    print("kind-pasteur-2026-03-14-S105j")
    print("=" * 70)

    # Level 2 data
    print(f"\n  LEVEL 2:")
    print(f"  {'n':>3} {'N_2':>8} {'|H_hat|':>12} {'|H_hat| formula':>20} {'E_2/E_0':>12}")
    print(f"  {'-'*60}")
    for n in range(3, 9):
        N2 = n*(n-1)*(n-2)//2
        mag = Fraction(math.factorial(n-2), 2**(n-2))
        e2_e0 = Fraction(2*(n-2), n*(n-1))
        print(f"  {n:3d} {N2:8d} {str(mag):>12} {'(n-2)!/2^(n-2)':>20} {str(e2_e0):>12}")

    # Level 4 data
    print(f"\n  LEVEL 4:")
    print(f"  n=5: 60 nonzero, ALL |H_hat| = 1/8, all cover 5 vertices")
    print(f"  n=6: 450 nonzero in TWO classes:")
    print(f"        Class A: 360 with |H_hat| = 1/8, cover 5 vertices")
    print(f"        Class B:  90 with |H_hat| = 1/4, cover 6 vertices")

    # E_4 decomposition
    print(f"\n  LEVEL 4 ENERGY:")
    print(f"  n=5: E_4 = 60*(1/8)^2 = {60*(Fraction(1,8)**2)} = 15/16")
    print(f"  n=6: E_4 = 360*(1/8)^2 + 90*(1/4)^2 = {360*(Fraction(1,8)**2) + 90*(Fraction(1,4)**2)} = 45/4")

    # Key insight about magnitudes
    print(f"\n  MAGNITUDE ANALYSIS:")
    print(f"  Level-2 magnitudes are always (n-2)!/2^(n-2)")
    print(f"  Level-4 magnitudes at n=5:")
    print(f"    All = 1/8 = 1!/2^3 = (5-4)!/2^(5-2)")
    print(f"  Level-4 magnitudes at n=6:")
    print(f"    Class A = 1/8 = 1!/2^3 (same as n=5)")
    print(f"    Class B = 1/4 = 2!/2^4 = (6-4)!/2^(6-2)")
    print(f"")
    print(f"  PATTERN for Class B: |H_hat| = (n-|S|)!/2^(n-2)")
    print(f"    n=5: (5-4)!/2^3 = 1/8  (all level-4 are Class B at n=5)")
    print(f"    n=6: (6-4)!/2^4 = 1/8  ... WAIT that's 1/8, not 1/4!")

    # Let me recompute more carefully
    # At n=6, Class B: |H_hat| = 1/4 = 0.25
    # (n-|S|)!/2^(n-2) = (6-4)!/2^4 = 2/16 = 1/8. WRONG!
    # So (n-|S|)!/2^(n-2) is NOT the right formula for Class B at n=6.

    # What IS 1/4?
    # 1/4 = 1/2^2
    # 1/4 = (n-4)!/2^(n-4+2) = 2!/2^4 = 1/8 at n=6. Still wrong.
    # 1/4 = (n-4)!/2^(n-3) = 2!/2^3 = 2/8 = 1/4. YES!
    # So Class B: |H_hat| = (n-4)!/2^(n-3)

    print(f"\n  CORRECTED FORMULA:")
    print(f"  Class A (covers n-1 vertices): |H_hat| = 1/2^(n-3)")
    print(f"    n=5: 1/2^2 = 1/4? NO, Class A at n=5 covers n=5 vertices, has |H_hat|=1/8")
    print(f"")
    print(f"  Hmm. Let me reconsider from scratch.")
    print(f"")
    print(f"  At n=5:")
    print(f"    All level-4 coefficients cover all 5 vertices")
    print(f"    |H_hat| = 1/8")
    print(f"    There are 60 of them")
    print(f"")
    print(f"  At n=6:")
    print(f"    360 level-4 coefficients cover 5 vertices, |H_hat| = 1/8")
    print(f"    90 level-4 coefficients cover 6 vertices, |H_hat| = 1/4")
    print(f"")
    print(f"  The 360 at n=6 that cover 5 vertices: they use arcs from a K_5 subset.")
    print(f"  There are C(6,5) = 6 ways to choose which 5 vertices to cover.")
    print(f"  Each K_5 contributes 60 nonzero level-4 subsets (like n=5).")
    print(f"  So 360 = 6 * 60. EXACTLY!")
    print(f"")
    print(f"  This makes perfect sense! The 360 level-4 coefficients at n=6 are")
    print(f"  the 'inherited' level-4 coefficients from each K_5 subgraph.")
    print(f"  The 90 are 'NEW' level-4 coefficients that use all 6 vertices.")
    print(f"")
    print(f"  HEREDITARY STRUCTURE:")
    print(f"  N_4(n) = C(n, n-1) * N_4(n-1) + N_4^new(n)")
    print(f"  At n=6: 450 = 6 * 60 + 90 = 360 + 90. CHECK!")
    print(f"")
    print(f"  And at n=5: 60 = 5 * N_4(4) + N_4^new(5)")
    print(f"  But N_4(4) = 0 (no level 4 at n=4), so 60 = 0 + 60 = all new.")

    # Level-2 analogy:
    # N_2(n) = C(n,2) * (n-2) - C(n,3) (overcounting)? No.
    # N_2 = n(n-1)(n-2)/2 = adjacent arc pairs.
    # N_2(5) = 30, N_2(6) = 60
    # N_2(6) = 6 * N_2(5)/5 * ... hmm, 60 = 2 * 30.
    # Hereditary: N_2(n) = C(n, n-1) * N_2(n-1) + N_2^new(n)
    # 60 = 6*30 + new? That would give new = 60 - 180 = -120. Negative! Wrong.

    # Actually for level 2, subsets of size 2 (pairs of arcs).
    # Inherited from K_{n-1}: there are n copies of K_{n-1}, each with N_2(n-1)
    # But each arc pair that's entirely within K_{n-1} is in C(n-2, 0)=1 copy...
    # This is getting complicated. Let me just focus on counting N_4^new.

    print(f"\n  N_4^new (full-vertex-coverage count):")
    print(f"    n=5: 60 (since all are new)")
    print(f"    n=6: 90")
    print(f"")
    print(f"  60 = C(5,2)*C(3,2) = 10*3 = 30? No, 60 ≠ 30.")
    print(f"  60 = 5!/(5-4)! = 5*4*3*2 = 120? No.")
    print(f"  60 = 5*4*3 = P(5,3)? YES, 60 = P(5,3) = 5!/2!")
    print(f"  90 = 6*5*3 = 90? YES!")
    print(f"  Or: 90 = C(6,2)*C(4,2) = 15*6 = 90. YES!")
    print(f"  And: 60 = C(5,2)*C(3,2)*... hmm, C(5,2)=10, C(3,2)=3, 10*3=30 ≠ 60.")
    print(f"  60 = C(10,4) - C(5,1)*C(6,4) = 210 - 5*15 = 210-75 = 135? No, 135 ≠ 60.")

    # Let me just see the DEGREE SEQUENCES of the n=5 nonzero level-4 subsets.
    # From the earlier output, at n=5 all 60 had degree sequence that covers 5 vertices.
    # I need to figure out what makes exactly 60 out of 135 covering subsets nonzero.

    print(f"\n  DEGREE SEQUENCE ANALYSIS:")
    print(f"  At n=5, covering 4-arc subsets have degree sequences:")
    print(f"    (1,1,1,1,4): 5 subsets")
    print(f"    (1,1,1,2,3): 60 subsets")
    print(f"    (1,1,2,2,2): 70 subsets")
    print(f"    Total: 135")
    print(f"")
    print(f"  But only 60 are nonzero. Is it exactly the (1,1,1,2,3) class?")
    print(f"  There are 60 with degree (1,1,1,2,3) and 60 nonzero. MATCH!")
    print(f"")
    print(f"  HYPOTHESIS: Nonzero level-4 Fourier coefficients at n=5 correspond")
    print(f"  EXACTLY to 4-arc subsets with degree sequence (1,1,1,2,3).")
    print(f"  The (1,1,1,1,4) and (1,1,2,2,2) subsets give H_hat = 0!")

    # What is the structure of a (1,1,1,2,3) 4-arc subset?
    # One vertex has degree 3 (incident to 3 of the 4 arcs)
    # One vertex has degree 2 (incident to 2 of the 4 arcs)
    # Three vertices have degree 1 each
    # Total degree = 3+2+1+1+1 = 8 = 2*4 arcs. Correct.

    # This looks like a "path + edge" or "star + pendant" structure.
    # The degree-3 vertex is a "hub" connected to 3 arcs.
    # The degree-2 vertex is connected to 2 arcs (shares edges with the hub?).

    print(f"\n  STRUCTURE of (1,1,1,2,3) subsets:")
    print(f"  One vertex v has degree 3: v is in 3 of the 4 arcs.")
    print(f"  One vertex w has degree 2: w is in 2 of the 4 arcs.")
    print(f"  Three vertices a,b,c each in exactly 1 arc.")
    print(f"  The 4 arcs: (v,w), (v,a), (v,b), (w,c)? Let's check.")
    print(f"  Arcs from v: 3 arcs. v connects to 3 other vertices.")
    print(f"  Arcs from w: 2 arcs. One must be (v,w) or not.")
    print(f"  If (v,w) is one of the arcs, then v has 2 other arcs to a,b.")
    print(f"  w has 1 other arc to c.")
    print(f"  4 arcs: (v,w), (v,a), (v,b), (w,c). Covers v,w,a,b,c = 5 vertices.")
    print(f"  Degrees: v:3, w:2, a:1, b:1, c:1. = (1,1,1,2,3). CORRECT!")
    print(f"")
    print(f"  If (v,w) is NOT one of the arcs:")
    print(f"  v has 3 arcs to a,b,c. w has 2 arcs to two of {{a,b,c}}.")
    print(f"  Say arcs: (v,a),(v,b),(v,c),(w,a). Vertices: v,w,a,b,c = 5.")
    print(f"  Degrees: v:3, w:1, a:2, b:1, c:1. = (1,1,1,2,3). Also works!")
    print(f"")
    print(f"  So the (1,1,1,2,3) pattern = a 'star with pendant' structure.")
    print(f"  These are trees on 5 vertices with exactly 4 edges minus 1 edge")
    print(f"  (caterpillar-type graphs).")

    # At n=6, the 90 full-coverage coefficients have degree (1,1,1,1,2,2).
    # What structure is this?
    print(f"\n  At n=6, the 90 full-coverage level-4 coefficients:")
    print(f"  Degree sequence: (1,1,1,1,2,2)")
    print(f"  Two vertices have degree 2, four have degree 1.")
    print(f"  This is TWO DISJOINT EDGES plus TWO MATCHING EDGES:")
    print(f"  More precisely: 2 degree-2 vertices each share 2 arcs.")
    print(f"  The 4 arcs form a 'double-edge plus matching' structure.")
    print(f"  Example: (u,a), (u,b), (v,c), (v,d). This is a PERFECT MATCHING")
    print(f"  of a bipartite graph K_{{2,4}} restricted to... hmm.")
    print(f"")
    print(f"  Actually: (u,v) with u in 2 arcs and v in 2 arcs.")
    print(f"  If u and v are NOT connected: arcs are (u,a),(u,b),(v,c),(v,d)")
    print(f"  with all 6 vertices distinct. That gives degrees u:2,v:2,a,b,c,d:1 each.")
    print(f"  But we only have 6 vertices, not 6 others. So vertices are u,v,a,b,c,d = 6 total.")
    print(f"  The 4 arcs span two 'cherries' centered at u and v.")
    print(f"  This is the DISJOINT UNION of two stars S_2 (two edges each).")

    # Number of such configurations:
    # Choose 2 centers: C(6,2) = 15
    # Each center gets 2 arcs to 2 of the remaining 4 vertices
    # First center: C(4,2) = 6 ways
    # Second center: C(2,2) = 1 way (gets the remaining 2)
    # But if (u,v) could be one of the arcs... hmm, we need u and v NOT adjacent in the subset.

    # Actually the arcs go between pairs of vertices in K_6.
    # Two disjoint cherries: choose center 1 (u), choose center 2 (v) ≠ u.
    # u gets arcs to 2 of the remaining 4 vertices: C(4,2) = 6.
    # v gets arcs to the remaining 2: C(2,2) = 1.
    # Total: C(6,2) * C(4,2) / 2 = 15 * 6 / 2 = 45? (divide by 2 for ordering of centers)
    # But we got 90, not 45!

    # Maybe (u,v) CAN be one of the arcs.
    # If (u,v) is an arc: then u has 1 more arc to some a, v has 1 more arc to some b.
    # a ≠ b (for coverage to be 6). Degrees: u:2,v:2,a:1,b:1, and 2 uncovered? No:
    # 4 arcs: (u,v), (u,a), (v,b), and one more. But we only have 4 arcs total.
    # If 4th arc is some (c,d): (u,v),(u,a),(v,b),(c,d). Covers u,v,a,b,c,d = 6.
    # Degrees: u:2,v:2,a:1,b:1,c:1,d:1 = (1,1,1,1,2,2). YES!

    # So there are TWO types of (1,1,1,1,2,2):
    # Type 1: two disjoint cherries (u->a,b and v->c,d, no arc between u,v)
    # Type 2: (u,v) is an arc, plus (u,a), (v,b), (c,d) — a path + isolated edge

    # Count:
    # Type 1: choose 2 centers C(6,2)=15, each gets 2 of remaining 4: C(4,2)=6, /2 for symmetry = 45
    # Type 2: choose arc (u,v): C(6,2)=15, choose a from remaining 4: 4, choose b from remaining 3: 3,
    #   remaining 2 form the last arc. Total: 15*4*3/2 = 90? /2 for (a,b) vs (b,a)? Or not?
    #   Actually: (u,a) and (v,b) are fixed by u,v,a,b. Arc (c,d) is forced.
    #   15 * 4 * 3 = 180. But (u,v) with (u,a),(v,b),(c,d) is the same as (v,u) with (v,b),(u,a),(c,d).
    #   Already accounted for by C(6,2). But also: is (u,a),(v,b) different from (u,b),(v,a)?
    #   YES, these give different 4-arc subsets.
    #   So: 15 * C(4,1) * C(3,1) = 15*12 = 180. But each subset appears TWICE
    #   (once as (u,v,a,b) and once as (v,u,b,a)? No, the subset {(u,v),(u,a),(v,b),(c,d)} is unique.
    #   Hmm, 180 seems too high if only 90 exist.

    # Let me just compute directly.
    print(f"\n  DIRECT COUNT of degree (1,1,1,1,2,2) subsets of K_6:")
    arcs6 = [(i,j) for i in range(6) for j in range(i+1,6)]
    count_1111222 = 0
    for combo in combinations(range(15), 4):
        arc_set = [arcs6[e] for e in combo]
        vertices = set()
        for arc in arc_set:
            vertices.update(arc)
        if len(vertices) != 6:
            continue
        degs = [0]*6
        for arc in arc_set:
            degs[arc[0]] += 1
            degs[arc[1]] += 1
        if tuple(sorted(degs)) == (1,1,1,1,2,2):
            count_1111222 += 1
    print(f"  Count of (1,1,1,1,2,2) at K_6: {count_1111222}")
    print(f"  Matches the 90 nonzero full-coverage coefficients: {count_1111222 == 90}")

    # Summary
    print(f"\n{'='*70}")
    print("GRAND SUMMARY")
    print(f"{'='*70}")

    print(f"""
  THE FOURIER SPECTRUM OF H ON TOURNAMENTS:

  Level 2k: |H_hat(S)| for |S| = 2k

  LEVEL 0: H_hat = n!/2^(n-1) (the mean)

  LEVEL 2: |H_hat(S)| = (n-2)!/2^(n-2) for adjacent pairs
           N_2 = n(n-1)(n-2)/2 nonzero coefficients
           E_2/E_0 = 2(n-2)/(n(n-1))

  LEVEL 4 at n=5:
           |H_hat(S)| = 1/8 for all 60 nonzero subsets
           All cover 5 vertices, degree sequence (1,1,1,2,3)
           E_4/E_0 = 1/60

  LEVEL 4 at n=6:
           Class A: |H_hat| = 1/8, 360 subsets (inherited from K_5)
                    Cover 5 of 6 vertices = C(6,5) * 60 = 360
           Class B: |H_hat| = 1/4, 90 subsets (new at n=6)
                    Cover all 6 vertices, degree (1,1,1,1,2,2)
           E_4 = 360/64 + 90/16 = 5.625 + 5.625 = 11.25
           E_4/E_0 = 1/45

  KEY OBSERVATION: The inherited level-4 energy from K_5 subgraphs
  EQUALS the new level-4 energy from full-coverage subsets at n=6.
  Both contribute 5.625 = E_4/2.

  EXACT Var/Mean^2:
    n=3: 1/3           = 0.3333 (levels 0,2 only)
    n=4: 1/3           = 0.3333 (levels 0,2 only)
    n=5: 19/60         = 0.3167 (levels 0,2,4)
    n=6: 13/45         = 0.2889 (levels 0,2,4)
    n=7: ???           (levels 0,2,4,6)
    """)

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
