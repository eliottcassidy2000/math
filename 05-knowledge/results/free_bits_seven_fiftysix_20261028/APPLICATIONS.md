# Applications, surviving constraints and new helper questions

## What is actually unified

Use a **labelled constraint carrier**: variables naming source choices, algebraic
relations among them, legal-domain restrictions, exact coefficients/weights,
and named boundary ports. Elimination produces a relation or weighted function
on the remaining variables. A coordinate is free to omit only after an exact
certificate shows it is determined, irrelevant to the declared target, or a
transported gauge. This unifies the operations used here; it does not assert
that every problem has a small representation or a short certificate.

Incoming Session85 proves Boolean-ideal certificates for nonlinear legality,
rank certificates for affine constraints, and integer Smith obstructions for
signed incidence. Our extensions add weighted box transforms and the explicit
56 matching/Frobenius map. The organizing conjecture worth investigating is:
**native proof families may have low-complexity boundary relations even where
no low-complexity scalar invariant exists.** That is a research programme,
not an established universal mathematical theorem.

## Collatz: an actual integer-clock test, not an analogy alone

Let T(n)=n/2 for even n and(3n+1)/2 for odd n. For each n in1,...,128 retain
its seven-step parity word b and endpoint T^7(n). The parity-word map is a
bijection to{0,1}^7. A standard inductive proof is that each extra parity bit
selects one of the two lifts of the previous residue class; the affine slope
has odd numerator. The checker also verifies the whole finite bijection.

Use columns1,...,7 of the Fano check matrix on the **chronologically ordered**
parity bits. All eight cosets give identical distributions of any three parity
coordinates. Yet their source-specific descent counts differ:

| syndrome | sources | first descent by step7 | endpoint below original source | sum of endpoints |
|---:|---:|---:|---:|---:|
|0|16|14|14|2484|
|1|16|15|12|1450|
|2|16|14|12|957|
|3|16|13|12|858|
|4|16|14|12|1187|
|5|16|14|12|1301|
|6|16|15|12|1031|
|7|16|15|12|1160|

**FINITE-EXACT**, with all source integers in receipts.json. This is not the
ROOT measure or an asymptotic theorem. It is a native hostile to replacing a
source-labelled packet by its complete three-bit marginal deck. Permuting
parity coordinates generally changes the integer clock and is not an allowed
Collatz symmetry. Equal marginals cannot pay a marked source's missing charge.

The new box construction explains analytically how such failures can persist:
a low-order deck may miss an entire extreme orthant. The useful positive route
is to carry syndrome-indexed **actual weighted** messages and then prove a
bound for each source/terminal fibre. Small state rank without bounds on message
values is insufficient.

## Hadwiger–Nelson: the orientation bit is observable by an exact unit-distance relation

For the oriented Fano basis(1,2,4), the matching from FIFTYSIX.md is

    M={{0,7},{1,3},{2,6},{4,5}},
    M_reverse={{0,7},{1,5},{4,6},{2,3}}.

Place labelled port1 at(0,0), port3 at(1,0), and each other port i at(10i+10,0).
The induced unit-distance graph has exactly the edge{1,3}. Let R6 be its proper
six-colouring relation on all eight ports. Contract R6 with the equality
kernel of M or M_reverse. The exact counts are

    <R6,B_M>=0,     <R6,B_M_reverse>=6^3*5=1080.

The first matching equates the adjacent ports; the second leaves four colour
variables with one inequality. This is an authenticated geometric **positive
control**: the orientation bit affects a real six-colour boundary relation.
It is deliberately not a six-chromatic gadget and proves no new lower bound
for the plane. A quotient retaining only the direction quadrangle and marked
point would identify these two matchings and give the wrong extension result.
All unary colour marginals alone are also blind to this distinction.

For actual forcing, compose the orientation kernels with the current core's
full boundary relation, preserving shared ports and exact unit distances.
A Petersen minor, a matching count or an abstract parity gadget is not a
substitute for that extension relation.

## Square sums and graceful labels: recover the stronger prior result first

[Session42](../session42_shared_lattices_square_colors_20261023/SQUARE_WORD_RECONSTRUCTION.md)
already proves: an ordered sum word of a permutation1,...,N has one possible
realization for odd N and at most two for even N. The second, when present,
is an alternating divisor-block shift. A marked vertex at a declared position
removes the ambiguity; an unordered mark need not. First square-sum twins
occur at36, while the stronger palindromic-clock threshold is48. These are
inherited results, not new claims of this packet.

For general connected bipartite graphs, the same first-moment argument works
when the two shores have different sizes; for equal shores the second moment
gives at most two candidates. The all-size reconstruction is elementary,
but realizability still needs positivity, distinctness and every square edge.
The new checker replays the36 example and finds further twins at48,56,64.
It treats capped search failures as unknown, never as nonexistence.

For graceful labels, retain the signs of edge differences and the absolute-
value collision test of THM-2761. Session85 proves that suppressing a two-edge
sum path produces a difference equation; a graph minor without that signed
offset loses the original integral label relation. Arbitrary-arity equations
can have Z/k torsion, so “one parity bit” is not a general replacement for the
integer obstruction group.

## Planar Jacobian: valuation codes need a global completion

The exact closest mechanism is
[THM-2685](../../../01-canon/theorems/THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate.md).
For squareclasses on a normal affine complex variety, divisor-valuation parity
is a binary code, but its kernel fits into

    0 -> units/(units)^2 -> unramified squareclasses -> Cl(X)[2] -> 0.

Thus erased boundary rows and hidden unit/class coordinates are distinct.
On X=G_m, t has zero valuation at every divisor **inside X** but is a nonsquare
unit; its double cover persists. On A1_C, every rational function with even
valuation at every finite prime is a square, since units are nonzero complex
constants and the class group is zero. The omitted boundary changes the
conclusion. This is the geometric analogue of missing parity constraints,
with a global sidecar that ordinary finite graph cycles do not supply.

Our inverse weighted-box transform shows that a generic added probe can be
informationally complete. For Jacobian applications the probe must additionally
be a permitted polynomial observable on the original source and preserve the
Jacobian condition. The known rational-carrier exclusions cannot be bypassed
by freely choosing a measurement unavailable to the problem.

## Four concrete open helper problems

1. **Orientation-selector extension of the existing six-colour core (OPEN).**
   For its actual eight labelled ports, compute the56 contractions indexed by
   oriented Fano triangles. Find an exact unit-distance attachment that excludes
   one surviving orientation over a fixed marked quadrangle while retaining the
   other, without identifying ports or assuming their palettes. The one-edge
   control above tests the compiler; a useful solution must act on the core's
   current surviving relation. Cheapest test: compare the two indexed extension
   supports before doing any geometric search.

2. **Source-weighted separator complexity (OPEN).** For a specified native
   first-reset ROOT packet, retain source, endpoint and integer-clock labels and
   compute the least interface dimension/rank needed for its weighted continuation
   tensor across each refinement cut. Does a natural changing-modulus refinement
   have a uniform or logarithmic bound on that complexity **and** a uniform bound
   on the corresponding message norms? A rank bound alone is explicitly not the
   target. First test: reproduce the original charge exactly on a small packet,
   then exhibit a same-marginal hostile before trying deeper scales.

3. **Stable inverse orthant encoding (OPEN quantitative problem).** Fix n and
   allowed widths in[1,1+eta]. Optimize the minimum separation of their subset
   sums and the conditioning of the finite exponential-sum decoder. How many
   moments or Fourier samples recover a sign law to prescribed error? Theorem6
   settles exact uniqueness; it leaves precision and conditioning. An elementary
   pigeonhole bound already prevents uniformly separated2^n heights inside an
   interval of length O(n). This is a place where “a bit for free” can be priced
   rather than asserted.

4. **Arithmetic/gauge extension without invisible global bits (OPEN by family).**
   For a chosen current planar-Jacobian normalization, augment its divisor parity
   matrix by explicit generators of units/(units)^2 and Cl[2], then determine
   which proposed source observables kill the remaining classes. Require the
   observables to extend polynomially on that same source. Test G_m and A1 first;
   they must give opposite answers. This is an executable family question rather
   than a claim that valuation-zero implies triviality everywhere.

A separate arithmetic target inherited from Session85 is realizing arbitrary
binary matroids by distinct prime atoms. The finite15/16-prime constructions
are proved; their unbounded extension requires prime values of d k²-1, and
multiple zero columns are obstructed because only prime3 has square p+1.
