# THM-265: Anatomy of SRCP Failures — Sorted Profiles Lose Topological Information

**Status:** PROVED (exhaustive analysis at n=6)
**Filed by:** kind-pasteur-2026-03-21-S19

## Statement

At n=6, the third failure of the c3-only SRCP (profile [0,0,1,1,1,1,1,1,2,2,2,2,2,2,3] mapping to H in {33, 37}) persists even when c5 per-arc counts are added. The (c3, c5) SORTED profile is identical for both H values:

c3 profile: [0,0,1,1,1,1,1,1,2,2,2,2,2,2,3] (same)
c5 profile: [3,3,3,3,3,3,3,3,3,3,4,4,4,4,4] (same)

Yet the two tournaments have different:
- 3-cycle vertex sets: {0,1,3} and {0,2,3} in H=33 replaced by {0,3,4} and {0,3,5} in H=37
- Omega edge count: 65 (H=33) vs 64 (H=37)
- alpha_2: 1 (H=33) vs 2 (H=37)
- H gap: 37 - 33 = 4 = exactly one additional independent cycle pair

## The Mechanism

SORTING the per-arc counts discards the information about WHICH arcs have which count. Two tournaments can have identical sorted (c3, c5) profiles but different UNSORTED profiles — different topological arrangements of cycles on the arcs. This topological difference changes which cycle vertex sets overlap, changing Omega's edges, changing alpha_2, and thus changing H by 4.

## Key Facts

1. The SRCP does NOT determine Omega(T). At n=6, one (c3,c5) profile maps to 2 distinct Omegas.
2. The SRCP does NOT determine H at n=6 with (c3, c5) alone (1 failure with H gap = 4).
3. The failure is because sorting loses topological information.
4. The UNSORTED profile (which is NOT isomorphism-invariant) WOULD distinguish the two tournaments.
5. At n=7, the sorted (c3, c5) profile has 0 collisions in 20K samples — so the topological information lost by sorting is somehow recovered at n=7 (perhaps because more arcs provide more constraints).

## Implications for the Conjecture

The SRCP conjecture ("sorted full profile determines H") requires:
- Either enough cycle lengths that the sorted profile becomes sufficient despite the topological loss
- Or a deeper mechanism that compensates for the sorting loss at larger n
- The conjecture FAILS in the strict sense at n=6 with (c3,c5) but may hold with (c3,c5,c7) or at all n with enough cycle lengths

## Source
srcp_proof_strategy_s19.py
