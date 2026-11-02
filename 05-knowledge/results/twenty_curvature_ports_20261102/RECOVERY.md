# Arithmetic and unit distance connections with their retained data

The useful common structure is an operation together with its realization
data. This session recovered several genuine maps involving the requested
numbers. Repeated cardinalities alone do not supply such a map.

## Twenty from arithmetic, geometry and path incidence

**PROVED, inherited:** write s(n) for the sum of proper divisors. The
factor-pair compression

    ab=n,   s(n)=1+a+b,   a,b positive integers

holds exactly when tau(n)=4 or s(n)=n+2; the alternatives are disjoint.
Among numbers having five proper divisors, its unique solution is 20.
See P1-P3 of [the divisor overlap classification](../aliquot_divisor_overlap_20261028/THEORY.md).
Here the compression uses (a,b)=(1,20), while the literal proper divisors
are 1,2,4,5,10. Thus s(20)=22 and their sum after removing 1 is 21.
The unit-pair compression is not a new literal divisor list.

Sixteen fails this condition: s(16)=15, and no factor pair of 16 has sum 14.
Prime cubes such as 8 satisfy the compression, even though a convention
giving prime powers priority would put them in the first bin. Preserve the
overlapping predicates rather than reading the user's lists as disjoint
intrinsic classes.

The prime-guarded row n=4p gives s(n)=3p+7 and excess 7-p:

| p | n | s(n) | Excess |
|---:|---:|---:|---:|
|3|12|16|4|
|5|20|22|2|
|7|28|28|0|
|11|44|40|-4|
|13|52|46|-6|

All five divisor posets have shape a three-element chain times a
two-element chain. Only 20 passes the compression. Exponent shape, or
the Mobius/Liouville labels on that unweighted shape, therefore cannot
replace actual prime weights.

Twenty is also near-perfect with redundant divisor 2:
20=1+4+5+10. The terminology and the broader family
2^a(2^(a+1)-3), when the odd factor is prime, are classical;
650 is a redundant-divisor-2 example outside that family.
See [Pollack and Shevelev, definition and Proposition 1](https://www.pollack-math.net/pnp-4.pdf).

Incoming [Session141's aliquot classification](../session141_metric_aliquot_square_paths_20261102/aliquot_twenty/PROOFS.md)
adds the complete fibre sigma^-1(42)={20,26,41}. We replay it exactly:
sigma(n)>=n+1 for n>1 bounds the search by41. Thus20 and26 share their
total divisor sum42, while their aliquot sums22 and16 differ. The same
incoming packet extends the excess-two classification through tau<=17;
that extension is credited there and is not needed for the current proof.

**CITED:** the finite Euclidean unit-distance extremum has u(10)=20,
u(20)=54 and u(21)=57. These count maximum edges, not rooted-ball
vertices; see [Alexeev, Mixon and Parshall, Theorem 1 and Table 1](https://arxiv.org/html/2412.11914v2).
Their theorem also reports 72<=u(25)<=78 and 76<=u(26)<=84.
These are the paper's bounds, not a claim to have surveyed every later
improvement.

**PROVED / FINITE-EXACT, inherited HN route:** [THM-4661 / odd-degree carriers](../../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md)
gives a twenty-vertex analytic lower bound for a minimal exterior that
blocks an old six-colouring, under its field and odd-extension hypotheses.
Its separately audited finite SAT result raises the necessary lower
bound to 22; it is not proof-checked and does not construct a failed
exterior. The low-vertex blocks are edges and odd cycles of length at
least five. The actual colour lists, old neighbours and field are
essential attachment data. MISTAKE-622 records the corrected near miss.

**FINITE-EXACT, this packet:** 20 directed square-sum paths occur at
25 vertices. At 27 one three-path compatibility component instead has
20 common fixed arcs; its variable incidence gives a 16-state controller.
Neither occurrence is the aliquot equality in disguise. The shared
research operation is to retain the incidence or weights that distinguish
otherwise identical coarse structures.

**Arithmetic correction:** 704-694=10, not 20. The cycle-free inverse
Collatz basin at 5 still exceeds a honeycomb radius-21 ball. At radius
26 the corresponding comparison is 2263>2107 for the triangular lattice.
Both are replayed from actual inverse edges in verify.py.

The [current Session140 synthesis](../session140_collatz_hexagonal_geometry_20261102/README.md)
also proves proper induced Euclidean unit-distance realizations of
countable locally finite pseudoforests, including ordinary Collatz.
They permit arbitrarily close nonadjacent vertices and crossings. Thus
unit-distance realization is compatible with failure of bounded local
lattice routing; the graph's chromatic number remains three independently
of convergence. This is not a high-chromatic witness for the plane.

## Factor pairs and Lorentz coordinates

There is an elementary exact map

    (a,b) -> (X,Y)=(a+b,b-a),   X^2-Y^2=4ab.

For n=20 its unordered factor pairs map as follows:

| (a,b) | (X,Y) |
|---|---|
|(1,20)|(21,19)|
|(2,10)|(12,8)|
|(4,5)|(9,1)|

After division by 2 sqrt(n), these lie on the unit Lorentz hyperbola.
Writing b/a=exp(2t) gives (X,Y)/(2sqrt(n))=(cosh t,sinh t).
The guarded factor transfer (a,b)->(pa,b/p), when p divides b,
changes t by -log p and is a rational Lorentz boost.
This is a genuine multiplication/addition coordinate bridge.
It does not preserve the aliquot sum without the whole divisor sidecar.

The pairs (4,11) and (4,13) map respectively to (15,7) on norm 176
and (17,9) on norm 208. They also appear in a different signed Collatz
construction, so the carrier and operation must be specified.

## Two distinct signed greedy maps

Recover the definitions from [the 3-adic G-map study](../collatz_mod6_20260917_three_adic_g_map.md):

    G_plus(m)=(2^k m-1)/3,
      k least nonnegative with 2^k m in {4,7} modulo 9;
    G_minus(m)=(2^k m+1)/3,
      k least nonnegative with 2^k m in {2,5} modulo 9.

The checked edges with clocks are

    G_minus: 4 --k=3--> 11 --k=0--> 4,
    G_plus: 13 --k=0--> 4 --k=0--> 1.

These explain 3*4-1 and 3*4+1 in that particular setting. Ordinary
Collatz has C(4)=2. No signed-map cycle is being transferred to C.
The geometric 11/13 circle identity and 3/5/15/17 boost composition are
derived separately in [GAPPED_GEOMETRY](GAPPED_GEOMETRY.md).

## Seven, fifteen, twenty-four and twenty-six

The [twin Pythagorean connection board](../twin_pythagorean_124_20261026/CONNECTION_BOARD.md)
already distinguishes the literal sum s(8)=1+2+4=7, the exponent orbit
{1,2,4} of multiplication by 2 modulo 7, and the Paley7 representation
split 1+3+3. There is no representation-theoretic 1+2+4 split asserted.

The [7/56 atlas](../free_bits_seven_fiftysix_20261028/ATLAS.md) supplies
several older, more informative maps: seven nonzero F2^3 parity-check
columns, fifteen nonzero words of the four-dimensional Hamming code,
four-port matchings, and a norm-seven Eisenstein torus with its seam.
The norm-seven quotient's K7 is not a planar Euclidean unit clique.
Likewise a faithful integer-addressed honeycomb drawing retains more
than an unlabelled hexagon picture; nearest-neighbour capacity has
already failed for the actual Collatz inverse tree.

The orders 15=16-1, 24=25-1 and 26=27-1 are also the multiplicative
group sizes of F16, F25 and F27. Frobenius acts on exponent labels by
multiplication by the characteristic, after a generator is chosen.
That is a typed algebraic use of these integers. It does not identify
those groups with the path banks. In this packet the actual directed
square-sum path counts at orders 25 and 26 are 20 and 24, while the
minimal controller sizes are 33 and 71. Their meaning comes from
incidence and operations, rather than an equality of displayed counts.

## The earlier 12/13 and Zagier question

The source-preserving trace repair in [TRACE_CARRIER](../cursor_cm_phase_20261028/TRACE_CARRIER.md)
and [the independent replay R1](../scaled_count_fft_carriers_20261030/THEORY.md)
is reproduced here. For u^2-a*u-1=0, the two primitive traces
(Tr x, Tr(ux)) have determinant a^2+4.

- a=5 modulo13: u^12=1; the sampled trace clock has 13 values on
  169 sources; the primitive pair distinguishes all 169.
- a=2 modulo11: u^12=-1; the sampled history has 11 values on
  121 sources; the primitive pair distinguishes all 121.
- a=5 modulo29: the determinant vanishes; even the primitive pair
  has only 29 values on 841 sources.

The [Ramanujan and Zagier recovery](../cursor_cm_phase_20261028/RAMANUJAN.md)
separates this finite-ring loss from the reciprocal Fricke two-sheet
quotient. G2 now gives another explicit loss and repair: adjacent dyadic
boundary rays collapse unbounded tree distances, while a positive gap
recovers a coarse inverse bound. These are parallel proof mechanisms
with different maps and fibres; neither establishes a new Zagier identity.

## What transfers to the open targets

For Hadwiger-Nelson, the next decisive objects are boundary colour-list
relations with realizable coordinates. A hyperbolic placement does not
keep Euclidean unit lengths. For LRC(14), the analogous obligation is
to preserve the actual marked times and their weights; a finite symbolic
quotient alone gives no all-source loneliness bound. For Collatz, G6
handles a selected return family, leaving coverage of every integer open.
For tournament metagraph enumeration, C1 is an exact compiler for a
specified operation bank, while arbitrary new ports require refinement.

The actionable synthesis is therefore a two-stage calculation: first
retain the geometry or signed incidence that makes operations legal,
then minimize the resulting observer with source weights attached.
The 16-versus-17 triangle and the dyadic-versus-gapped drawing show
exactly why the first stage cannot be replaced by an unweighted graph.
