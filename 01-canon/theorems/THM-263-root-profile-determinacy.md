# THM-263: Root Cycle Profile Determines H at n=5

**Status:** PROVED (exhaustive verification)
**Filed by:** kind-pasteur-2026-03-21-S18

## Statement

For tournaments on n=5 vertices, the sorted root cycle profile uniquely determines H(T).

**Definition:** The root cycle profile of a tournament T on [n] is the sorted tuple (c(alpha_1), ..., c(alpha_m)) where m = C(n,2) and c(alpha_{ij}) = number of directed 3-cycles passing through the arc corresponding to the positive root alpha_{ij} = e_i - e_j of A_{n-1}.

**Theorem:** At n=5, the sorted root cycle profile is a complete invariant for H. There are exactly 9 distinct profiles mapping bijectively to H values {1, 3, 5, 9, 9, 11, 13, 15, 15} (note: two profiles give H=9 and two give H=15, reflecting score sequence differences within the same H).

## Proof

Exhaustive verification over all 2^10 = 1024 tournaments on 5 vertices. The 9 profiles are:

| Profile | H |
|---------|---|
| [0,0,0,0,0,0,0,0,0,0] | 1 |
| [0,0,0,0,0,0,0,1,1,1] | 3 |
| [0,0,0,0,0,1,1,1,1,2] | 5 |
| [0,0,0,1,1,1,1,1,1,3] | 9 |
| [0,0,0,1,1,1,1,1,2,2] | 9 |
| [0,0,1,1,1,1,2,2,2,2] | 11 |
| [0,1,1,1,1,1,1,2,2,2] | 13 |
| [1,1,1,1,1,1,1,1,1,3] | 15 |
| [1,1,1,1,1,2,2,2,2,2] | 15 |

## Boundary

This property is SHARP at n=5 and FAILS at n=6:

- **n=6:** 3 out of 29 profiles map to 2 H values each. All failures have H-gap exactly 4, corresponding to one alpha_2 contribution in the OCF. The 3-cycle profile captures alpha_1 completely but cannot see alpha_2 (pairs of vertex-disjoint odd cycles).

- **n=7:** 51 out of 149 sampled profiles have multiple H values. Failure is endemic.

## Lie-theoretic interpretation

The root cycle profile encodes how the signed A_{n-1} roots participate in triangular cycles. At n=5, every odd cycle is either a 3-cycle or a 5-cycle, and the 5-cycles are determined by the 3-cycle structure (since the complement of a 3-subset in [5] is a 2-subset, which is too small for an independent cycle). So the 3-cycle root profile contains ALL cycle information at n=5.

At n >= 6, independent cycle pairs become possible (alpha_2 > 0), and 5-cycles carry information not captured by 3-cycle profiles. The H-gap of 4 is exactly 2^2 = the OCF weight of one additional independent cycle pair.

## Source
root_profile_determinacy_s18.py, root_profile_5cycle_s18.py
