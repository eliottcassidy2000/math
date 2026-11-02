# Review and exact reproduction

Status: same-agent mathematical audit plus independent consequence
implementations inside one script. No separate reviewer or formal build.

From the repository root:

```sh
python3 05-knowledge/results/clocks_carriers_genus_20261102/verify.py
python3 -O 05-knowledge/results/clocks_carriers_genus_20261102/verify.py
python3 agents/check_docs.py
git diff --check
```

Normal and optimized JSON must agree byte for byte; no Python assertions
are used as checks. These finite checks supplement the written universal
arguments, not the other way around.

## Declared universes and consequences

1. All fixed-path tilings of orders 3..7, at every valid h, have exactly the
   claimed corner-state cardinalities. Actual endpoint deletions decode
   every corner coordinate. There are 268,560 deletion updates and 200,832
   attachment updates (all input rows, orders through 6). Reflection and
   complement are checked. All fixed-path labelled HP languages through 6
   are distinct, supplementing the inherited all-order reconstruction theorem.
2. The inverse Collatz grammar is expanded through depth 16, with 65,536
   distinct sources at that depth and 131,070 verified macro edges. Every
   native step is replayed, not just the affine endpoint identity. The
   residue table proves two choices for every odd unit; source labels never
   collide. The unbounded growth proof is separate from this finite expansion.
3. Genus identities are checked at 17,257,65537, with the composite F5
   hostile. An independent implementation builds actual permutations S,R
   on P1(F_p), checks their product is translation, and computes cycle
   Euler characteristic at 12 prime levels through 65537. Four-channel
   arithmetic is checked at k1..8. No unknown Fermat primality is assumed.
4. All connection sets at odd orders 3,5,7,9,11,13,15,17 are enumerated
   under their full unit group. The 15/17 orbit counts agree at 16, with
   free multiplier actions. Explicit permutations check the 375- and 1215-
   element block automorphism subgroups of two order 15 tournaments.
   No full composite-order isomorphism classification is inferred.
   The auxiliary Fermat-product exponent is checked at 4,096 increments,
   including the exact equality boundary; 129 Pascal parity rows give an
   independent product encoding. All 31 nonempty subsets of the known
   Fermat primes have precisely the four predicted count collisions.
   The proof covers unknown Fermat primes through monotonicity on ALL
   auxiliary binary indices; these finite controls alone would not do so.
5. There are 694 exhaustive single/two-leaf attachment controls on all
   graphs of orders 2..4. Independent subset DP gives Q15's endpoint counts
   exactly at 8 and 9; a direct edge check gives Q18's three leaves.
   The bridge theorem is proved for arbitrary graphs, not inferred from this.

## Proof audit and repaired scope

- Corner counts exclude the shared ordered labels and clock metadata;
  six orientation bits are not six total bits for a self-contained labelled
  object. Every transition explicitly decrements the budget.
- The final re-marking hostile uses n5,h1 and a legal cost 1 re-marking,
  so it tests an enlargement at the SAME clock. An earlier n4,h0 example
  illustrated a fresh query only and was replaced before publication.
- Maximal macro length 7 does not imply disjoint native paths: even tails
  can merge. The growth proof uses distinct odd source labels, not disjointness.
- Contracting the trivial Collatz cycle does not remove the growth obstacle.
  The hyperbolic layout is for macro edges, with native clocks retained.
- Genus is computed for X0, not X0+. The prime-level formula is not applied
  to composite F5; its actual genus is calculated separately.
- The group action at composite15 is a multiplier action. Free multipliers
  do not imply full automorphism group C15; the block examples refute that.
- A000055/A000568 relate through marked tree and orientation data. The
  inherited Burnside-mark criterion and corrected pair/triple boundary are
  stronger than an unmarked numerical analogy and remain the dependencies.
- The complete free-action collision classification excludes n=1. Its
  auxiliary exponent has a three-way plateau at indices 0,1,2; removing
  index 0 leaves the first pair, 3/5. Later plateaus are isolated pairs.
- The existing square-sum all-n construction is not rediscovered or weakened
  into a finite observation. The new question is marked-port response.

No false mathematical theorem from this packet was promoted. The input's
modular-curve conflation and the prior navigation's stale n56 sentence are
recorded in the correction ledger, with their strongest survivors.

## Canonical dependencies (ID plus exact path)

- [THM-796-three-sorted-recursive-tiling-line-node-incidence](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md)
- [THM-871-fermat-rung-rigidity](../../../01-canon/theorems/THM-871-fermat-rung-rigidity.md)
- [THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation](../../../01-canon/theorems/THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md)
- [THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md)
- [THM-4836-coherent-carrier-sets-of-tournaments-and-the-pair-or-triple-law](../../../01-canon/theorems/THM-4836-coherent-carrier-sets-of-tournaments-and-the-pair-or-triple-law.md)
- [THM-4801-collatz-meets-the-constructible-polygons-trunk-five-initial-segment-principle-collatz-closed-family](../../../01-canon/theorems/THM-4801-collatz-meets-the-constructible-polygons-trunk-five-initial-segment-principle-collatz-closed-family.md)
- [THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion](../../../01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md)
- [THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet](../../../01-canon/theorems/THM-4505-small-graphs-encoding-arithmetic-square-sum-zigzag-fences-collatz-alphabet.md)
- [THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags](../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md)
- [THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections](../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md)
- [FRACTRAN_PRECISION](../../../05-knowledge/results/stability_priority_moduli_20261102/FRACTRAN_PRECISION.md)
- [FRACTRAN_CLOCKS](../../../05-knowledge/results/stability_priority_spectrum_20261102/FRACTRAN_CLOCKS.md)
- [GUARDED_PATHS](../../../05-knowledge/results/stability_748_guarded_encodings_20261102/GUARDED_PATHS.md)
