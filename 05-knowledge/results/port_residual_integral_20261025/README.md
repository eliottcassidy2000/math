# Closed endpoint-four repairs and integral tournament recurrences

2026-10-25. **PROVED, independently reviewed:** the two remaining local
port branches; their exhaustion of all closed endpoint-four lines from
order15; an integral normalization of the adaptive clone recurrence;
pair-free and positive-interpolation parity laws; and a bounded-cycle
algorithm for fixed binary precision of A000568. **FINITE-EXACT** controls
have their own declared universes. Full black-core connectivity and a
smaller unrestricted structural counting recurrence remain **OPEN**.

## The results and their scope

- [Closed endpoint-four theorem](PORT_FRONTIER.md): every closed blue
  marked path with endpoint outdegrees (n-5,4), n>=15, has a literal
  three-black detour. There is no extra strong-target, mixedness, module,
  disjoint-port or isomorphism filter. The inherited head uses24 complete
  affine systems; the infinite tail is structural. The order-seven
  pure-blue target remains a genuine obstruction at smaller orders.
  A separate [finite certificate](LOW_ORDERS.md), checked by an
  [independent raw-arc solver](LOW_ORDERS_REVIEW.md), also closes order13.
  Its480 affine cells represent4,026,531,840 normalized labelled sources.
  The unchanged ledgers fail in two order12 and four order14 cells;
  these are template failures, not connectivity obstructions.
- [h=3, 1->3](h3/README.md) and [h=4, 3->1](h4/README.md): each new
  branch works for every n>=12 with incoming ports {h,n-3,n-2,n-1} at
  vertex0. The first uses induced degrees along a reversal-product
  permutation; the second has three explicit arc contradictions.
  Both retain actual sources, word legality, blackness and the exact
  path-edge XOR. [h3 review](h3/INDEPENDENT_REVIEW.md),
  [h4 review](h4/INDEPENDENT_REVIEW.md).
- [Integral adaptive recurrence](counting/README.md): the already stored
  defects give b=max ceil(D(d)/(d-1)). The normalized moments
  Y_r=O(u) 2^(-rb) W_r have integer transfer matrices on an explicit
  product lattice. All shifts are nonnegative; destination gauges make
  conversions and mergers commute. The same sparse/dense state graph
  and scalar rank are retained. [Independent review](counting/INDEPENDENT_REVIEW.md).
- [Fixed-precision moment quotient](counting/MODULAR.md): modulo2^K,
  every channel Y_r with r>=K is the same tail value. The allowed
  transitions preserve this module, reducing adaptive storage to
  min(h,S,K) weight coordinates per nonterminal fibre. Ordinary states
  were already scalar; the reduction affects clone-excess banks.
  [Independent review](counting/MODULAR_REVIEW.md). Absolute K-bit
  precision and normalized K-bit precision remain different requests.
- [Parity and positive interpolation](parity/README.md): the user's
  odd-action identity extends through exact canonical clone inflation.
  If c_N counts pair-free classes, c_N is odd precisely for positive odd
  N, with c_0=1. The positive interpolation R_N is 1+u modulo2 at odd
  N and u+u^2 at positive even N. [Independent review](parity/INDEPENDENT_REVIEW.md).
- [Fixed binary precision](precision/README.md): after removing the
  proved baseline power of two, K binary digits of A000568 need only
  O(sqrt(K)) permutation cycles, independently of n. This filters the
  arithmetic Burnside sum; it does not enumerate literal tournament
  cores. It also proves an even-order totient leading-bit formula and
  normalized modulo-four formulas. [Independent review](precision/INDEPENDENT_REVIEW.md).

## Inheritance and the three requirements

The anchor is exact port transport. The niche is integral adaptive
counting. The wildcard recovers the April
[THM-305 valuation theorem](../../../01-canon/theorems/THM-305-two-adic-valuation-tournament-count.md)
and June
[THM-514 unit-tail kernel](../../../01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md)
as a route to a smaller modular observable. No historical novelty claim
is made for the recovered results.

The closest prior mechanisms are the unprotected-port and finite-horizon
moment proofs in the
[preceding packet](../unprotected_ports_cancellation_20261025/README.md).
Its missing branches become the first two rows above. The corrected near
miss is the h=3 fixed insertion that repeats a vertex or requests a
backward adjacent arc; the later insertion used here changes the cut.
The canonical counting hostile is a genuine past5-cycle whose ungauged
moment transition has a negative binary exponent. The least-used sidecar
gcd(c,d)>=1 supplies the needed bound on the hidden cycle count.

| Live concept | Integral coordinates | Preserved operations | Positivity / lost information |
|---|---|---|---|
| Marked port repair | Literal arc bits and three path words | Every K move and exact labelled target | No clone quotient is taken; no counting reduction follows merely from connectivity |
| Adaptive moments | Odd-factorial and defect-normalized moments | Transfers, stagewise mergers and decreasing-horizon continuation | Signed excess remains; individual histories and tournament objects are forgotten |
| Odd subset action | Integral orbit polynomial P and correction B | Constant XOR descends only for invariant gap sets | B>=0 supplies weighted counts, but C5 forbids the proposed two-point equivariant quotient |
| Canonical clone inflation | Actual weight-orbit counts before reduction | Isomorphism and canonical module reconstruction | The binomial transform is valid modulo2; its integer version first fails at order4 |
| Fixed binary precision | Valuation plus odd denominator modulo2^K | Addition and multiplication at the declared precision | High-valuation terms vanish only for that observable; restoring more digits needs more terms |

This makes the user's three requirements separate, testable statements.
The coordinate change fixes integrality; it does not manufacture a
positive recurrence. The parity collapse is a smaller observable; it
does not supply independent toggles. Positive coefficients alone do not
construct an action-preserving geometric or combinatorial quotient.

## Concurrent connection and next boundary

Incoming commit f093b00695 contains the independently audited
[all-odd-action parity and denominator-precision interface](../session63_integral_observers_root_arrangements_20261025/cancellation/RESULT.md).
Its Theorem1 includes the user's general modulo-two law, which is
recovered here as an input, not claimed as a new discovery. Its Theorem2
explains why reducing an averaged numerator before paying its denominator
loses information. Our integral gauge addresses the analogous issue
at every adaptive transition; the fixed-precision computation strips
the even denominator before reducing and inverts only its odd part.
The incoming synchronized-action spectrum also reinforces the retained
operation boundary: modular depth is not independent-toggle rank.

The [live-merge connection ledger](INCOMING.md) also records the incoming
Fourier positivity criterion, sparse-carry conductors and regular-action
depth formulas, retaining each source's status and operation boundary.

The next port target is the still-unsettled low-order head and other
endpoint/open-path layers, not the exhausted closed n>=15 family.
For counting, the constant-tail theorem now supplies one exact
precision-dependent continuation quotient. Further reduction would need
the remaining congruences in the actual moment module. A
[reverse valuation potential](counting/REVERSE_POTENTIAL.md) also gives a
proved construction for normalized ordinary digits: a reverse graph
pass pays the baseline before residue arithmetic. Its implementation in
the production engine remains open; the extra pass is not free. No
connection to Collatz, RH, HN or LRC is asserted from these tournament
results.

A precise next compression probe is the Newton basis
F_j(X)=X product_(0<=i<j)(X-2^i) on the moment nodes X=2^(m-b).
Test its divisibility filtration and the transported batch operations
together. Additional congruences between the lower coordinates cannot
be discarded merely because high channels become redundant. This is
an OPEN follow-up, not an additional compression claimed here.

## Reproduction

From the repository root:

```sh
python -B 05-knowledge/results/port_residual_integral_20261025/replay.py
python -B 05-knowledge/results/port_residual_integral_20261025/verify.py
python -B -O 05-knowledge/results/port_residual_integral_20261025/verify.py
```

The replay runs producers and independent controls under normal and
optimized Python, compares deterministic receipts, and performs no
mutation without --refresh. [AUDIT.md](AUDIT.md) records exact universes,
independence, inherited proof inputs and remaining limits.
