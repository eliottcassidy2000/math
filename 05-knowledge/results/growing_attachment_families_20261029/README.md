# Whole families of arithmetic attachments

2026-10-29. **PROVED deductions from CITED classical inputs / FINITE-EXACT
controls / OPEN dynamic closure. Same-agent proof and implementation review.**
This packet grows families by their generating operations and gives uniform
laws for all compatible finite configurations within explicit size budgets.
It makes no literature-priority claim and creates no new canon theorem ID.

The main gain is a pair of charged interfaces: **rank, factor degree and
deleted points** over finite fields; **rank, conductor and source height**
for prime-indexed attachments. Displayed query count is often the wrong
cost. Joint dependencies and realization losses remain part of the object.

## Results

All proofs and boundaries are in [PROOFS](PROOFS.md).

| Result | Scope and mechanism |
|---|---|
| P1: every finite polynomial character bank | Factor-parity code of rank r; Weil gives an explicit count for every compatible answer. The same relative bound covers every Boolean union, with no union penalty. |
| P2: actual rooted Paley composition | Sequentially adjoining a legal root adds exactly one independent linear factor. Uniform joint colored-copy bounds hold for every root placement outside the deleted set. |
| P3: every finite binary code/matroid | Products of r actual attachment signs realize any affine code of dimension r, once the explicit field-size inequality holds. All products and supported powers are covered together. |
| P4: growing prime banks | Siegel--Walfisz permits conductor Q<=log^A X; the inherited quadratic large sieve then handles all compatible colored tournament patterns in a growing sufficient window. |
| P5: one infinite Apollonius discriminant family | Centers (0,0),(t,0),(0,1), radii 1,2,3; all four-character even-parity words have their predicted density outside explicit bad characteristics and zero parameters. This counts discriminant signs, not yet distinct geometric circles. |
| P7: every finite legal 3n+1 or 3n-1 valuation word | An explicit determinant product identifies all exceptional primes for that chart. Outside them, its intermediate quadratic-sign patterns obey the same character bound, with the original integer source and clock retained. |

For the general polynomial-bank Paley theorem, a convenient sufficient
relative-error budget is

    (D+z+h+k+1) 2^(r+h+k) = o(sqrt(q)),

where D is factor-support degree, z the number of deleted points, h named
roots and k ordered new vertices. In the undecorated case this includes
h+k <= (1/2-eta)log_2 q for every fixed eta>0. The one-root estimate is
classical (Graham--Spencer; Bollobas--Thomason; Jones); its sequential use
here retains the bank and the actual root identities.

For prime banks Q<=log^A X, the sufficient growing-pattern window is
k <= (1-eta)sqrt(log_2 X). Neither window is claimed optimal. The prime
asymptotic constants are ineffective; finite checks do not establish a
usable starting height.

## Inheritance and portfolio

**Anchor:** lift fixed arithmetic attachments to whole growing families.
**Niche:** the older Apollonius discriminant specialization problem.
**Wildcard:** apply the same theorem to every finite valuation word while
retaining source and clock, rather than hypothesizing independent clocks.

Concept board: factor-parity codes; quantitative character cancellation;
actual ordered roots; excluded sources; common Frobenius transport;
valuation-word affine lines. Comparing these lanes produced P7: a fully
dependent affine value pair can have nearly independent quadratic signs.
That is a change of observer with measurable information loss.

Closest proved mechanism: Session107 [CHARACTER_STRUCTURE, A/B](../session107_quasirandom_arithmetic_20261029/characters/CHARACTER_STRUCTURE.md),
and [CONDITIONAL_CANCELLATION](../session107_quasirandom_arithmetic_20261029/CONDITIONAL_CANCELLATION.md).
These already give fixed banks, arbitrary large-fiber spectral counting,
and exact affine source clocks. They are inherited, not rediscovered.

Canonical hostile: the twelve-vertex 81/27 pair, inherited from
[DEFECT_TRANSFER D6](../overnight_rank_return_20261029/DEFECT_TRANSFER.md).
Equal cyclic-triple data still do not identify classes. P6 keeps the
factor-three orbit normalization; no response matrix was recomputed.

Corrected near miss: square multiplication preserves a sign only away
from its new zeros; root reuse can make a supposedly new bit dependent.
An entire family of square queries can have rank zero and empty source
domain. Also, fixed-modulus PNT does not by itself justify a growing bank.

Least-used relevant sidecar: [Apollonius squareclass code, sections 3--6](../apollonius_squareclass_code_20261019/APOLLONIUS_SQUARECLASS_CODE.md).
Its generic graph-incidence rank did not prove numerical distribution.
P5 supplies an explicit one-parameter application, preserving the
distinction between discriminant data and geometric realization.

Recovered context: [LEGENDRE_SIDECAR](../root27_machine_bridges_20261026/anchor/LEGENDRE_SIDECAR.md)
already records the known countable random tournament construction and
why extension existence alone does not imply ordered-prefix densities.
The [prior arithmetic packet](../quasirandom_arithmetic_structure_20261029/README.md)
supplies the simplex/Hamming duality and Sylvester hostile. Incoming
f8bc2099e7 was read before derivation; P6 connects its colored recurrence
without claiming a new A000568 term.

## Exact evidence

[check.py](check.py) uses only Python's standard library. Reproduce:

```text
python -X utf8 -B check.py --output checks.json
python -X utf8 -B -O check.py --output optimized.json
```

Paths may be absolute; the script works from any current directory.
[checks.json](checks.json) records:

* 38,192 complete root configurations and 400,376 sign-pattern queries:
  q=3,7,11,19,23 through four roots (capped by q), and q=27,31,43 through
  three. All combinations and all signs in those universes are included.
* Every one of the 16 binary subspaces of F2^3 on seven field sizes,
  including all Boolean unions of its basis-answer fibers.
* 2,186 literal ordered tuples, independently counted by permutation/cut
  filtering and successive bitset intersections, in four rooted/colored
  examples including F27.
* A literal F27 Frobenius transport check, plus a counterexample when
  vertices move but the named root stays fixed.
* The Apollonius family at seven primes through 3511, all eight even words;
  the Hamming/simplex dual; exact CRT fibers; square-zero and reused-root
  failures.
* All 240 words of lengths 1--4 with valuations in {1,2,3}, for both fixed
  signs, and 816 legal prime charts among 7,11,19,43,139. Another 384 charts
  are rejected by the determinant certificate. Direct integer trajectories
  independently verify every retained valuation and intermediate value.

Normal and optimized receipts agree; see [AUDIT](AUDIT.md). These are finite
checks of the stated envelopes and guards, not proofs of convergence.

## Transfer boundaries and next targets

| Source -> target | Preserved | Lost / required sidecar |
|---|---|---|
| Polynomials -> affine code | Every multiplicative sign relation | Z, factor support and field size determine actual realization. |
| Root incidences -> product queries | Joint sign vector and its code | A product query is not itself a new tournament vertex. |
| Uniform source cylinder parameter -> sign patterns | Actual integer valuation word and residue sampling law | Intercepts, source density and word order cannot be discarded. |
| Discriminants -> parity code | Square/nonsquare data | Geometric degeneracy and distinct-solution guards remain. |
| Colored arithmetic samples -> class frequencies | Labelled counts and symmetry weights | No new exact enumeration or metagraph complement law. |

The next sharp questions are whether actual complement compilers preserve
a controlled rank/degree budget, whether the growing-pattern windows can
be enlarged, and how many extra divisors complete the Apollonius geometry
guards. General smallest-carrier and all-size selected-path recovery
remain open. Source files and primary-source locators are pinned in
[sources.json](sources.json).
