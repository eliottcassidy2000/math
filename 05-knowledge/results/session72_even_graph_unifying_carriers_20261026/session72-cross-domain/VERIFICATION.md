# Exact verification receipt

Input pin: `4950589d056a64f8807774cf7da893dcafcd63a0`. All arithmetic tests are local, deterministic and standard-library only. Snapshot verification is mandatory, checking both SHA256 and the original Git blob identity for all13 inputs. No whole-repository validation is claimed by this isolated lane.

## Commands and outcomes

- `python3 -B check.py > check.normal.json` — PASS.
- `python3 -O -B check.py > check.optimized.json` — PASS.
- `cmp check.normal.json check.optimized.json` — byte-identical.
- `python3 -B check.py --inject-failure` — expected nonzero exit1.
- `python3 -O -B check.py --inject-failure` — expected nonzero exit1; optimized execution does not disable the guard.

The injected-failure receipt is `hostile.normal-optimized.json`. The paired arithmetic receipts contain all exact counts and selected native automorphism results.

## Closed finite universes

1. **Cyclic alignment.** Every period profile of arity1..4 with each period1..4 (340 profiles), plus `(4,6,8)`, `(3,5,15)`, `(4,6,10,15)`, `(6,10,15)`, `(3,9,15,45)`, `(15,21,35)`, `(3,9,9)`:347 profiles,45,520 states,2,250 orbit signatures. Direct orbit generation is compared with pairwise phase signatures. Every realized signature is reconstructed prime by prime.
2. **Image compatibility.**21,534 candidate edge-label assignments, exhaustively enumerated for each selected profile whose edge-label space has size at most50,000. Triangle congruences are compared with actual orbit signatures, not assumed valid by construction.
3. **Sparse observation.** Every simple edge subgraph for all347 profiles:17,100 graphs. Actual signature projection is injective precisely when each prime-power-induced graph is connected. The ordinary-star and all-three-pairs hostiles are retained.
4. **Finite fields.** All8 binary normal-basis words in F8. The map `F8 tensor F8 tensor F8 -> F8^9` is represented by a27-by-27 binary matrix of rank27. All729 tensor basis-pair products preserve multiplication in every factor; a rank-only coincidence would not pass.
5. **Native tournaments.** All62 r-bit colour blocks for r1..5 have the expected uniquely decoding strong-component sequences. All104 constructed tournaments for `(r,m)=(1,3),(2,3),(1,5)` are checked.24,272 vertex-pair minimal-module closures recover the outer carousel blocks from adjacency alone. Direct degree-pruned permutation search counts all automorphisms of all8 twelve-vertex cases and4 twenty-four-vertex controls:12 direct censuses. No generic full automorphism census at larger sizes is claimed.
6. **Arithmetic cocycle.** All961 pairs of binary words of lengths0..4; exact carry composition. Actual valuation-word merge483/469->17. The110 legal macros from source7 and credit1, including the macro86 `5 -> 2 -> 9` reference-shape refill. Corrected full-match sources5 and21 at word1000.
7. **Geometry.** All5,461 split words of lengths0..6 over totals1,2,3,4, tested on13 half-integral inputs from-1 to5:70,993 exact legality/output comparisons. Equal formal identity/unequal domain hostile, plus the independent auditor's same-formal-centre/unequal-partial-commutation hostile (`2x` on `(0,1)` and `3x` on `(0,2)`). All14,641 ordered pairs of integer vectors in `[-5,5]^2` for the Gaussian/Lorentz identity; the norm65 hostile is checked separately. Two-involution equal-product/different-connectivity hostile. Five-square local-seam/overlap hostile, with area sum5 and union area4. This is not a dissection search or complete Penrose geometry replay.
8. **Older moment mechanism.** Exact Bareiss factorial-Hankel determinants for degrees0..10; p-by-p and (p+1)-by-(p+1) ranks for p2,3,5,7,11,13; cyclic13 alias and kernel/nonideal hostile. All proper subset guards for even/odd parity laws on cubes k2..7, with the actual full conjunction probability gap. Native tournament/Penrose realization is inherited from the pinned source reports, not silently reconstructed by this abstract probability check.

## Reproduction boundaries

All universal conclusions use the written proofs, not extrapolation of these finite ranges. `seal_sources.py` is a one-time local preparation helper and is **not part of the curated publication manifest**; replay needs only the immutable source snapshots, source ledger and check.py. External PDF/HTML bytes are not packaged. EXTERNAL_SOURCES.json pins only primary source URLs, checked locations, exact imported statements and non-implications.

The initial local checker caught an implementation issue in SCC sorting: Python empties the list during key calculation, so iterating that same list in the sort key erased the order statistic. The key now closes over a saved copy. The mathematical colour-decoding claim was unchanged; the repaired checker was rerun in both modes before this receipt was sealed.

Independent audit results and their separate exact universes are owned by the response/alignment auditor. Parent publication should bind that audit to the final hashes in MANIFEST.json rather than conflating its larger/different test sets with this producer receipt.

The independent auditor also corrected a prose scope issue in GA's explanatory corollary: common centres and commutation concern formal total affine extensions. Actual partial maps need fixed-point domain membership and equality of composite domains. The report now says this explicitly and the hostile is executed in both modes. No pre-correction theorem was published by this lane.
