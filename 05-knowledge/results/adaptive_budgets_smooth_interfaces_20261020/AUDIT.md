# Proof review and reproduction

**Same-agent adversarial review and separately implemented controls.** No
independent-agent review, formal proof assistant check or novelty claim is
made. The analytic theorem is a deduction using a cited classical theorem;
finite tests do not establish its asymptotics. Universal Collatz entry and
an improved plane-chromatic bound remain open.

## Proof gate

1. **Degree budget, not merely an abstract matrix.** All past lengths are
   qc with c odd and y-smooth; q>y makes their largest prime exactly q.
   Future lengths c remain unprocessed. The common inequality `(q+1)sum S<=N`
   makes every future subset affordable from every past subset. It also
   keeps every individual probe in every state's future vector.
2. **Injectivity is on the actual retained response.** The divisor-closed
   smooth bank has `G=Z diag(phi) Z^T`; the triangular incidence matrix is
   invertible. A radical-only bank would fail. Ordinary weights are positive;
   the lower bound is not exported to signed cancellation.
3. **Rank has fixed degree on both shores.** The pigeonhole cell fixes
   sum A, hence past degree q sum A and future degree sum A. The discarded
   scalars are nonzero row/column factors. The exponential dot-product
   kernel is strictly positive definite on distinct embedded points by
   its full polynomial-moment expansion. No rank claim is made after all
   future columns have been summed into one number.
4. **The cited smooth estimate is used in range.** With
   `y=exp(B_N/2)` and `z=exp((log N)/2-eta B_N)`, u tends to infinity and
   is less than `sqrt(y)`. Odd-part fibres cost only O(log log z) in log D.
   The positive fixed eta supplies room in the actual degree budget before
   it is allowed to tend to zero in the liminf argument. The old upper
   bound is uniform inside batches, not only at their endpoints.
5. **Clocks and denominator alignment.** All prefix denominators in the DP
   are exactly `2^K_j`. A shifted positive emission of +1 is multiplied by
   `2^K_j` after the step, not the previous denominator. The draft error
   is recorded in MISTAKES.md and corrected before any push. The separate
   suffix formula verifies positive charge minus reset dividends against
   every actual endpoint in the charge universe.
6. **Two-block padding is conditional on certified arrival.** The appended
   exponent-two word is physically valid at 1. The original source remains
   in the full appended cylinder. The actual integer and unaccelerated
   clocks include every repetition. This gives a compression theorem and
   a reformulation of convergence, not a proof of universal arrival.
7. **The positive-charge helper is OPEN.** Every tested first exit exists
   by direct finite traversal; its existence is not assumed universally.
   The unproved step drops a nonnegative dividend from an already paid
   expression. A random parity or independent-digit model is not used.
8. **Geometry is induced and labelled.** All 45 pairwise squared distances
   are tested; the exact 18-edge graph has nine named ports. The two-wing
   projection is proved by a constructive extension, and a same-palettes
   hostile proves that named rim positions matter. The gadget is itself
   three-colourable and cannot be iterated along unit edges to escape its
   common triangular lattice. It is not a counterexample to THM-4614's
   different field and fixed eight-port theorem.
9. **Half-plane scope.** Translation plus graph compactness gives equality
   only for unrestricted chromatic number. Moving an entire configuration
   is not adding a context at fixed absolute ports. Piecewise congruences
   may lose cross-piece edges, so Banach–Tarski supplies no colour transfer.
10. **Farey content.** Both fractions are reduced before the determinant-one
    test. The gcds divide adjacent integers and are coprime; this makes
    `gh=d` equivalent to the full prime-power condition `d|a(a+1)`. At d=9,
    the radical admits false events. These are graph edges, not a forced
    tournament orientation or a Collatz transition.

## Finite universes and controls

`verify.py --limit 200001` performs 301,194 explicit checks. It tests four
complete smooth banks, their Smith determinants, all subset responses and
equal-degree kernel submatrices. It enumerates the entire normalized
six-colour relation of the three-wing graph and directly checks its distances.
For all 100,000 odd sources 3..200001 it computes first-exit optimal block
budgets with exact integer DP and replays every selected envelope using
rationals. It also tests terminal padding and exact appended cylinders at
twelve named sources. The step caps are declared and none is hit.

The separate implementation imports no producer or inherited module. It
performs 1,306,782 explicit checks. It checks the same complete 100,000-source histogram, enumerating **all 80,909
partitions** needed to verify each hard source's minimum budget and its
reported optimal value. It rebuilds the geometry with independent triangular
lattice norm coordinates and enumerates all six-colour assignments directly.
It checks Farey content on a=1..500, d=1..200. For merged counts it compares
partition Burnside formulas with literal signed edge-orbit actions for every
permutation through n=6, reproduces the incoming n<=9 table, and tests the
twice-odd transition formula on all odd partitions of mass at most ten.

`charge_probe.py` uses a forward integer recurrence on all 500,000 odd
sources 3..1000001, with a 2000-step cap. The separate implementation instead
sums suffix contributions explicitly and checks both the exact signed
identity and the paid positive bound on the entire same universe. Its cap
is 1000 steps, which is still above the observed maximum 111. The largest
discarded-dividend/descent ratio is
`215402584767891/281474976710656`, at source63.

The full receipts retain the exact universes and results. Normal and
optimized Python runs are compared byte for byte; all load-bearing checks
use explicit exceptions rather than assert. File and dependency hashes are
in manifest.json. These are independently implemented calculations by the
same agent, not an independent referee's mathematical audit.

## Primary-source pins

* Hildebrand–Tenenbaum, *Integers without large prime factors*, Journal de
  Theorie des Nombres de Bordeaux 5 (1993), 411--484. The author's PDF at
  https://tenenb.perso.math.cnrs.fr/PPP/Survey_friables.pdf was inspected at
  Corollary 1.3 and its stated range on PDF page 7; Theorem 1.1 and (1.7)
  give a second route to the required estimate. No RH input is used.
* De Bruijn–Erdos, *A colour problem for infinite graphs and a problem in
  the theory of relations* (1951), Theorem 1, checked in the author archive:
  https://www.renyi.hu/~p_erdos/1951-01.pdf . The institutional copy returned
  HTTP 403, so no claim of inspecting that copy is made. The theorem adopts
  the usual choice framework; no measurable variant is imported.

## Corrected draft and stopping boundaries

The charge recurrence initially undercounted exponent-one emissions by a
factor two; source3's `(1,4)` word exposed it. The corrected full replay still
passes, but its worst ratio and conclusion use only the corrected receipts.
No erroneous result was pushed. The source27 anchor correction is inherited
from the previous packet and is not reintroduced here.

The new sharp width law closes only the designated prime-batch interface's
leading scale. The geometric relation is nontrivial but its common lattice
is three-colourable. The terminal-padding construction does not create a
missing first arrival. These limitations generate the next questions listed
in CONNECTIONS.md rather than claims of solving the original conjectures.
