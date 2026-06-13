# HYP-1834: LRC has an operation-shadow recursion

**Status:** EXPLORATORY; supported by S378 exact feature extraction and prior S356/S362/S365/S366 machinery.

## Claim

The Lonely Runner recursion in the number of moving runners is not governed by a
single scalar monotone.  It should be tracked by an operation-shadow state vector:

- additive/Dirichlet completion: the initial-segment equality spine;
- multiplicative/divisibility skeleton: speeds divisible by `n=k+1` and nonunit
  residue channels;
- boundary state: unit endpoint skeleton of size `phi(n)`;
- endpoint dynamics: protection pressure, peel depth, and terminal core size;
- speed-set operation closure: additive gates, multiplicative gates, and divisor
  edges inside the speed set;
- micro-staircase state: scalar-ramp distance, missed-cell packages, and repair
  deficit near the frontier.

In this framing, addition is the complete transitive shadow and multiplication
is the sparse surviving residue.  LRC moves between them: Dirichlet gives the
additive equality spine, while endpoint protection and descent are forced
through divisibility.

## Evidence

`natural_lrc_recursive_modes_s378.py` connects three exact ledgers.

First, the natural operation shadow repeats S365/S366 at a metric level:
addition on `[1000]` has all `499500` order edges, while nonunit multiplication
has only `5070`, density `0.010150`.  Product-sum minima aligned with LRC
runner count change by factor-packing type rather than smoothly in `n`.

Second, initial-segment LRC features through `n=22` show an arithmetical
recursion.  Prime `n` maximizes the unit skeleton (`phi(n)=n-1`), while
composite `n` opens nonunit quotient channels.  Peel depth and product-sum seed
type jump at different prime/composite transitions.

Third, the hard-family feature vectors separate tight examples from near
disproofs.  For example, the `n=14` seven-ladder has tiny positive gap ratio
`0.005411`, `84` boundary witnesses, five speeds divisible by `14`, and
quotient layer `higher_denoms_5_types`; the initial `n=14` spine is
boundary-only with unit quotient layer `unit_mod_14`.

## Predictions

1. LRC proof searches should rank speed sets by a vector, not only by max gap:
   `(max_gap_ratio, phi(n), divisor gates, nonunit residues, peel depth,
   protection pressure, scalar distance, missed-cell repair deficit)`.
2. Composite frontier cases such as `n=14` and `n=15` should be attacked by
   joint additive/multiplicative state transitions: unit-skeleton peeling plus
   divisibility-channel descent.
3. Counterexample searches should generate endpoint-protection patterns and
   operation-closure signatures first, then solve for integer speed sets.
4. A repo-native LRC TDA module should treat endpoint-protection incidence
   graphs and micro-staircase missed-cell hypergraphs the way `tournament_tda.py`
   treats SCC defect, deletion residue, and odd-cycle support.

## See

- `04-computation/natural_lrc_recursive_modes_s378.py`
- `05-knowledge/results/natural_lrc_recursive_modes_s378.out`
- `07-reflections/natural-lrc-recursive-modes-s378.md`
- HYP-1820, HYP-1826, HYP-1829, HYP-1830, HYP-1831, HYP-1833
