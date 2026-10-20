# Independent publication audit: finite-state debt barriers

Date: 2026-10-20. Status: **PASS with no further source changes required**.

The audited addendum is `FINITE_STATE_BARRIERS_ADDENDUM.md`. Its frozen
version already includes the integer `L>=1` boundary condition recorded in
`FINITE_STATE_INDEPENDENT_REVIEW.md`. This audit independently checks the
proofs and supplies a second self-contained exact computational path. It
does not rely on the earlier review's additional unaccompanied test counts.

The base guarded reset-debt packet, its publication manifest, and the
repository were not changed. This addendum audit is a separate add-only
packet; no theorem text overlay is needed.

## 1. Positive barriers: necessity and sufficiency

Theorem 1 is correct for a nonempty finite vertex set, finitely many edges,
positive real multipliers, and nonnegative real debts. Parallel edges,
loops, isolated vertices, acyclic components, and the empty edge set are
included.

The source/target orientation in the edge certificate is correct:

    lambda_e b_target - b_source >= delta_e.

For chronological paths, multiplying each slack by the product of preceding
multipliers makes the intermediate barrier terms cancel. The exact identity
is the displayed equation (2), with debt weighted in the same chronological
order as in the base note. On a cycle this yields

    (Lambda-1)b_start >= Delta >= 0.

Strict positivity of the barrier excludes Lambda below one. When Lambda is
one, the positive coefficients in Delta force each nonnegative edge debt
to be zero. A product-one cycle with positive debt is therefore genuinely
infeasible, not a missing strictness case.

Every closed walk decomposes into simple directed cycles, including loops.
Its product is the product of their multipliers. If each simple cycle has
product at least one, a product-one closed walk can contain only neutral
simple cycles; each of their debts vanishes. This justifies checking simple
cycles, including in a multigraph.

For sufficiency, the finite simple-cycle minimum gives a positive epsilon
on paid cycles. Subtracting epsilon from each positive-debt edge's logarithm
leaves all cycle sums nonnegative. Shortest-path distances satisfy
`d_v<=d_u+w_e`; their negations consequently satisfy
`r_u-r_v<=w_e`, with the required orientation. Exponentiation gives
nonnegative original slack, strictly positive on each positive-debt edge.
This remains true for paid edges not belonging to any cycle. Finite common
scaling then pays all debts. If there are no paid edges, scaling by one is
already valid.

## 2. Exact rational feasibility and constructive certificates

The equivalence with the rational system `b_v>=1` is correct. A positive
feasible vector may be multiplied by a constant at least one until its
minimum coordinate reaches one; nonnegative debts ensure that this does
not invalidate any inequality. Conversely, a solution of that system is
already a positive barrier. Rational coefficients and nonempty rational
polyhedra give rational certificates.

The fully rational construction is valid independently of that general
polyhedron fact. For N at least one,

    (1+1/N)^p <= 1+(2^p-1)/N,

by bounding each `N^(-j)` by `N^(-1)` in the binomial expansion. The stated
choice of N therefore leaves every modified cycle product at least one.
The added source reaches every vertex with multiplier one. Removing a
cycle of product at least one cannot increase a path product, so a minimum
product occurs on a simple path and is a strictly positive rational.

Writing that minimum as d_v gives `d_v<=d_u mu_e`. Thus
`b_v^0=1/d_v` has original edge slack at least
`(kappa-1)b_u^0` on paid edges, and at least zero elsewhere. Rational common
scaling produces a rational barrier. No floating-point equality test is
used or needed.

The independent checker compares the cycle criterion against rational
basic-feasible-vertex enumeration, not Fourier–Motzkin elimination. This is
a complete feasibility test here: any nonempty feasible polyhedron lies in
the coordinate orthant b>=1; minimizing sum(b) has compact nonempty
sublevel sets and an extreme-point minimizer. Enumerating nonsingular
sets of active constraints therefore finds a feasible vertex whenever one
exists. Each candidate is solved and checked with exact fractions.

For constructive barriers, the independent implementation uses
multiplicative Floyd–Warshall, rather than the producer's synchronous
Bellman–Ford relaxation. Both algorithms are applied only when exact cycle
products satisfy the theorem's hypothesis.

## 3. Strict cycles and the finite-head conclusion

The exponential multiplier argument is sound. Repeatedly removing simple
cycles from an m-edge walk leaves at most V-1 edges, and each removed cycle
has at most V edges. The finitely many simple cycles have a positive
minimum logarithmic product; the finitely many residual simple paths have
a finite lower bound on their logarithmic products. Equation (5) follows.
If there are no cycles, no infinite path exists and the claim is vacuous.

The finite-head conclusion concerns graph/block endpoints `x_m`. With
`B=max_v b_v` and Lambda tending to infinity, equation (3) gives
`x_m<=B+o(1)`. For integer x_m, the gap to `floor(B)+1` is positive, so
eventually `x_m<=floor(B)`. This is not a claim that every intermediate
literal or odd-step height within a block is in that head.

The finite head still needs separate analysis; nothing forces ROOT or even
excludes a non-target cycle in an arbitrary affine system. The abstract
fixed point `3 -> (3+3)/2 = 3` with positive multiplier gain is a valid
hostile. It is correctly not represented as a Collatz counterexample.

## 4. Finite reset banks: exact guard failure

Theorem 2 is correct with a finite nonempty bank of nonempty exact odd
valuation words, each containing a reset. Its maximum length L is at least
one. For every integer t>=1, the specified source
`n=2^(L+1)t-1` is odd and at least three. At each of the first L odd steps,

    3 U^j(n)+1 = 2[3^(j+1) 2^(L-j)t-1],

and the bracket is odd. Thus the valuation is exactly one, not merely at
least one. All these endpoints strictly increase and remain above ROOT.
Every proposed bank word therefore fails its valuation guard at its first
required reset. Additional guard restrictions cannot restore legality.

For any nonempty all-ones word, lambda is below one and physical height
debt is positive. It cannot satisfy the unit budget. It also increases
every positive source. Consequently both the unit-budgeted-bank corollary
and the genuinely descending-bank version have the claimed scope.

## 5. Finite-state universal coverage obstruction

Theorem 3 is correct with all its stated hypotheses: V>=1 finite states,
finitely many fixed nonempty valuation words, an integer length bound
L>=1, the actual physical height parameters, and positive state barriers.

For every t>=1 and source `2^(VL+1)t-1`, the first VL valuations are one.
Any lawful path of V graph edges would consume at most VL valuations.
Every edge word along it would therefore be all ones. A repeated state
gives a nonempty closed subwalk, whose multiplier is `(2/3)^K<1` for a
positive integer K, contradicting Theorem 1.

There is also a direct independent proof of the same contradiction:
an all-ones edge has lambda below one and positive debt, so its barrier
inequality forces the target barrier to exceed the source barrier strictly.
The all-ones subgraph is acyclic and admits at most V-1 consecutive edges.
This proof confirms the orientation and the exact V-edge cutoff.

The conclusion holds from every initial graph state, including a state
chosen as a function of the source. The witnesses do not reach ROOT during
the relevant prefix. The explicit L>=1 convention covers the empty-edge
case too, where taking an otherwise vacuous maximum L=0 would allow source
1 and invalidate the separate non-ROOT sentence. The frozen theorem
already has this boundary repair.

The result forbids universal continued coverage in this specific model.
It does not forbid source entry after an unbounded preliminary search,
unbounded word lengths through parameterized schemes, infinitely many
states or unbounded auxiliary state, a different signed/matrix resource,
or useful certification of restricted source classes. These are essential
restrictions, not optional interpretations.

## 6. Reproduced and independent finite controls

The producer checker reproduces **83,452 explicit checks**, byte for byte,
in normal and optimized Python. Both outputs match the supplied receipts.

The independent checker reproduces **146,537 explicit checks**, also
identically in normal and optimized Python. Its declared universe is:

- All 2,401 labelled two-state graphs with an absent edge or one of
  lambda in {2/3,1,3/2} and debt in {0,2} at each ordered vertex pair;
  448 are feasible. This parameter grid differs from the producer's.
- 400 seeded three-state and 50 seeded four-state multigraphs, allowing
  loops and parallel edges, with at most six edges each.
- Empty graphs on one through five vertices, neutral zero-debt loops and
  two-cycles, positive-debt neutral-cycle failures in both positions,
  positive-debt acyclic edges, and cycles containing an edge below one.
- 1,200 seeded rational path-slack identities, with arbitrary changing
  barriers, independently checking the exact excess identity.
- 640 actual all-ones source prefixes, lengths 1–128 and t=1–5,
  including exact valuations and strict growth above ROOT.
- All 1,359 reset-containing valuation words of lengths 1–5 with letters
  in {1,2,3,4}, checked against six uncovered sources with common L=5.
- All 6,561 two-state labelled physical graphs with an absent edge or one
  word in {(1),(2),(3),(1,1),(1,2),(2,1),(2,2),(1,1,3)} at each ordered
  pair. Exactly 1,375 satisfy the barrier criterion. Each feasible graph
  is given an exact rational barrier and checked for failure of a V-edge
  path at three all-ones witnesses, from every initial state.
- The explicit physical changing-barrier example with barriers (5,8),
  source 31, lawful first edge to 47, and failed next-word guard.

Random instances use the fixed seed 904712 and are finite-exact samples,
not class censuses. All computations use exact integers/fractions and
exception checks that survive optimization. Neither checker imports the
other. The all-graph and arbitrary-length claims rest on the proofs.

Run from the published package directory:

    python3 -B verify_finite_state_barriers.py > /tmp/finite-state-producer.json
    python3 -B -O verify_finite_state_barriers.py > /tmp/finite-state-producer-O.json
    cmp /tmp/finite-state-producer.json /tmp/finite-state-producer-O.json
    cmp /tmp/finite-state-producer.json finite-state-receipt.json
    python3 -B independent_finite_state_audit.py > /tmp/finite-state-independent.json
    python3 -B -O independent_finite_state_audit.py > /tmp/finite-state-independent-O.json
    cmp /tmp/finite-state-independent.json /tmp/finite-state-independent-O.json
    cmp /tmp/finite-state-independent.json finite-state-independent.normal.json
    sha256sum -c FINITE_STATE_SHA256SUMS

## 7. Publication verdict and limits

The frozen addendum, its producer checker and receipts, the earlier review,
and this independent checker/receipt package may be published unchanged as
an add-only extension of the guarded reset-debt packet. Source identity and
the base-note dependency are recorded in `FINITE_STATE_SOURCE_INVENTORY.json`.

No additional mathematical correction or source overlay is required.
There is no universal Collatz convergence, ROOT-entry, or infinite-bank
coverage conclusion, and no native signed-row 7X estimate is repaired.
The results are elementary finite-graph barrier theorems and an explicit,
guard-aware impossibility theorem for the stated finite physical-height
model. No external peer review or literature novelty is claimed.
