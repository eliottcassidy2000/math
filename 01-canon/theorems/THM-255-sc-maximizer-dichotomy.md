# THM-255: SC Maximizer Dichotomy at n=6 Regular

**Status:** PROVED (exhaustive computation, n=6)
**Filed by:** kind-pasteur-2026-03-20-S1

## Statement

For the regular score class (2,2,2,3,3,3) on n=6 vertices, the 2640 labeled tournaments
decompose into exactly 4 isomorphism types by their independence polynomial I(Omega(T), x):

| Type | IP | H | Count | Anti-aut? | alpha_2 | c3 | c5_dir |
|------|-----------|-----|-------|-----------|---------|-----|--------|
| A (SC-BIBD) | (1,14,4) | 45 | 240 | YES | 4 | 8 | 6 |
| B (SC-rich) | (1,20,1) | 45 | 240 | YES | 1 | 8 | 12 |
| C (SC-weak) | (1,16,2) | 41 | 720 | YES | 2 | 8 | 8 |
| D (NSC) | (1,19,1) | 43 | 1440 | NO | 1 | 8 | 11 |

Maximum H=45 is achieved ONLY by SC tournaments (Types A and B).

## Key Observations

1. **All regular n=6 tournaments have c3=8.** The 3-cycle count is a score-class invariant.

2. **Two routes to max H=45:** Both satisfy alpha_1 + 2*alpha_2 = 22.
   - Route A: High alpha_2 (4 disjoint 3-cycle pairs from orbit pairing), low alpha_1 (14)
   - Route B: High alpha_1 (20 directed odd cycles), low alpha_2 (1)

3. **NSC achieves alpha_1 + 2*alpha_2 = 21.** One unit short of the maximum, giving H=43.

4. **SC-weak achieves alpha_1 + 2*alpha_2 = 20.** Even worse, giving H=41 < H(NSC)=43.

5. **The anti-automorphism enables both extremes:** SC tournaments can achieve EITHER
   maximum disjoint pairs (Route A) OR maximum total cycles (Route B), but not both.
   NSC tournaments occupy a constrained middle ground.

## Route A Mechanism (Orbit Pairing)

For Type A, the involutory anti-aut sigma = (0,1)(2,3)(4,5) pairs vertices into 3 orbits.
The 4 vertex-disjoint 3-cycle pairs are:
- {0,2,4} + {1,3,5}
- {0,2,5} + {1,3,4}
- {0,3,4} + {1,2,5}
- {0,3,5} + {1,2,4}

Each pair picks one vertex from each orbit, and sigma maps each selection to its complement.
ALL 4 complementary pairs have both halves cyclic (comp_both=4).

## Route B Mechanism (Cycle Maximization)

For Type B, sigma distributes arcs uniformly, creating the maximum number of directed
5-cycles (12, vs 11 for NSC, 8 for SC-weak, 6 for SC-BIBD). With alpha_1=20 total
directed odd cycles (8 three-cycles + 12 five-cycles), Route B achieves H = 1+40+4 = 45.

## Proof

Exhaustive computation over all 2^15 = 32768 tournaments on 6 vertices. All 2640 regular
tournaments classified by H, IP, anti-automorphism status. Zero exceptions.

Script: `04-computation/sc_maximizer_orbit_fixed.py`
Results: `05-knowledge/results/sc_maximizer_orbit_fixed.out`

## n=7 Contrast

At n=7 regular (3,3,3,3,3,3,3), the mechanism FLIPS:
- H=189 maximizer: c3=14, disjoint_3_pairs=7 (LOWEST)
- H=175: c3=14, disjoint_3_pairs=14 (highest)
- H=171: c3=14, disjoint_3_pairs=10

The n=7 maximizer wins via alpha_1 (more total directed odd cycles), not alpha_2.
This means ANY algebraic proof of the SC Maximizer conjecture must handle both mechanisms.

## Related

- OPEN-Q-016: SC Maximizer conjecture (this theorem is the n=6 case)
- THM-024: Involution anti-aut existence
- THM-029: H=7 impossibility
- Kind-pasteur-2026-03-06-S18e: SC maximizer mechanism draft
