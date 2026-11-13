# Session 315: Dyadic pairing, exact self-converse gaps, and joint recurrences

This session joins the [previous metagraph incidence and score recurrences](../session314_metagraph_constraint_algebra_20261113/README.md)
to a two-vertex growth law. The common principle is **retain the marked
interface until the exact gluing is performed**. A symmetry may collapse
a full product to diagonal and paired terms, but forgetting the terms
that pair is not the same as counting them.

## Enumeration and metagraph

* [Paired sites and actual blue moves](PAIRED_SITES.md): each self-converse
  class is a simultaneous orbit of a half-size tournament, an independent
  undirected graph, and (in odd order) one bit per site. A two-state
  transfer gives the exact joint blue-line grade polynomial, and a
  signed strong self-converse inventory gives $1/\mathcal T(-z)$
  for the ordinary tournament-class series $\mathcal T$.
* [Half-cut Hamilton transfer](HALF_CUTS.md): full Hamilton counts
  include nontransversal cuts. After reflection, these give doubled
  terms; fixed cuts give quadratic forms or squares. The same transfer
  reads out the number of fixed blue paths.
* [Pointed Burnside recurrence and prime sieve](SPECTRUM_AND_SIEVE.md):
  one cycle-insertion recurrence gives A000568, the self-converse
  count, and the merged vertex count under different evaluations.
  Joint Hamilton/fixed-path fibers filter admissible site-cycle
  orders by the gcd of the two counts. The blue-line grading is not
  sufficient for the full black edge wiring.

## Exact counts versus residue counts

* [Four rooted twin charts and a mirrored wrap](DYADIC_SPECTRUM.md):
  every odd residue modulo every power of two occurs in a strong
  self-converse tournament. A fixed bank with arbitrary transitive
  toothpick lengths is still polylogarithmically sparse in absolute
  count. The result cannot be turned into an exact target claim by
  increasing the precision.
* [Global low-count self-converse gaps](SPECTRUM_AND_SIEVE.md): the
  missing positive odd values below 75 are
  $7,19,21,23,35,39,43,47,53,55,59,63$.
  Busch's strong bound and mirrored SCC decomposition make the
  result global for this finite target range. In particular the
  ordinary 23 and 63 must be sought outside this symmetry sector.
* [Exact bridge flip laws](MIRROR_BRIDGES.md): one mirrored return
  gives a difference of squares or an $x^2+2y^2$ norm. Two switches
  expose a real deletion-port interaction, not independent digits.

## Older structures brought into the construction

[Square-sum and prime filters](HISTORY_AND_CLOCKS.md) gives the
23–26 ordered two-cover aperture, simultaneous unit targets on
separate 23- and 25-path clocks at the same first-25-prime-sieve
surviving length, and different six-prime filters inside the Ogg set.
The Monster/Baby Monster valuation pattern is a literal five-cycle
with a separate 13 feeder under composite-base divisor valuation;
the actual 2A class index instead gives a two-cycle.

## Reproduction

Run `python replay_pairs.py`, `python replay_half_cuts.py`,
`python replay_cycle_pointing.py`, `python scsmall.py`,
`python bridge.py`, `python replay_dyadic.py`, and
`python replay_history.py` from this folder. These standalone controls
check the formulas and finite support used above. They do not
replace the all-order arguments.

**Open:** Chan–Pak coverage, a compressed recursion for the full
black incidence, and a changing-core construction that turns
a modular address into the specified integer with the necessary
endpoint data.
