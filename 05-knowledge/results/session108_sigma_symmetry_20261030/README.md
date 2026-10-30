# Five bits, three bits, and19: symmetry with its observer retained

Session108, started2026-10-29 and completed2026-10-30 UTC.
**PROVED elementary deductions / explicitly CITED classical classification /
independently audited finite computations.** No literature-priority claim,
Collatz proof, or general classification of square divisor sums.

Here sigma(n) includes every positive divisor; s(n)=sigma(n)−n is the aliquot
sum. A squareclass bit records the parity of a prime exponent in a divisor
sum, not a binary digit of the input integer.

## 1. The requested numbers have different exact roles

| Input p | sigma(p) | Its squareclass | s(p) |
|---:|---:|---:|---:|
|3|4|1, the zero vector|1|
|5|6|the 2-bit plus the 3-bit|1|
|19|20|the 5-bit|1|

Direct squareclass observation of(3,5,19) has rank3, sigma observation has
rank2, and their joint direct/sigma observation has rank4. The joint system
recovers the extra2-bit from direct3 and sigma5. Aliquot values alone identify
all three prime inputs at1.

There are two distinct five/three statements. In the earlier six-label bank,
direct characters needed five bits and sigma characters three. In the
zero-augmented seven-prime bank below, five is a **selection-kernel dimension**
and three is a **response-image dimension**. The same numbers do not make
these the same map.
[Code proof and exact observer ledger](codes/SYMMETRY.md).

## 2. A genuinely free3 toggle, and the correct extended code

The seven-prime bank(2,5,7,19,29,59,89) has Hamming selection kernel and
simplex response image. Adding the actual prime3 adds a zero sigma column:

    selection kernel [8,5,1],       response image [8,3,4].

Thus there are32 square-sigma selections and8 character response patterns.
This is not yet extended Hamming. Coprime multiplication by3 multiplies sigma
by4, preserving its squareclass; repeated multiplication across the same prime
port need not preserve it.

Adding prime-count parity as another observer gives the genuine[8,4,4]
extended Hamming code. The free3 choice uniquely corrects that parity.
The coordinate group changes from GL3(2), order168, to AGL3(2), order1344,
with exact joint-syndrome transport. This larger zero-fiber group does not
preserve the original sigma syndrome on every odd-parity selection; it is a
changed observer, not an unexplained gain in symmetry.

The19 port is a Fano point with stabilizer S4 of order24. Its three incident
lines are{19,2,59}, {19,5,29}, {19,7,89}; each gives a square product of sigma
outputs. Seven lines give21 point-line flags. These are explicit incidence
sets, not an identification with unrelated7/21 counts.

## 3. Arithmetic symmetries beyond the formal code

Using the cited classical square case of the Nagell–Ljunggren equation,
the only prime powers n>1 with square sigma(n) are

    3, 81, 343, with sigma values4, 121, 400.

Their prime ports are not independent:3 and81 use the same prime. On
{3^a7^b u : a∈{0,1,4}, b∈{0,3}, gcd(u,21)=1}, replacing exponent states
gives a faithful S3×C2 symmetry preserving sigma squareclass. It does not
preserve sigma's numerical value or an aliquot trajectory.

The complete neutral prime-power list is a cited theorem application, not
inferred from a finite search. Nonzero charges can also cancel, for example
sigma(5·19·29)=3600, so it does not classify all square-sigma composites.
Raw character queries need their unit guards: sigma(3) and sigma(81) are
both squares, but their Legendre symbols modulo11 are+1 and0.
[Classification, shared-prime correction and arithmetic actions](arithmetic/ARITHMETIC.md).

### Why19 is selected by an exact arithmetic return

Let Phi6(x)=x²−x+1. Among integers x>1 with **prime-power Phi6(x)**,

    sigma(Phi6(x)) is square iff x=2 or19;
    s(Phi6(x))=3x iff x=19.

The multiplier3 comes from the identity Phi6(x)+3x=(x+1)². [The June30 recovery](history/COVERING_BOUNDARY.md) already recorded Phi6(19)=7³; the present extension
is the audited divisor-sum selection and return proof. The surviving example is

    19 —Phi6→ 343=7³ —aliquot sum→ 57 —divide by3→ 19.

These are three different operations, not a Collatz or ordinary aliquot cycle.
The full prime-power claim uses the cited classification. In the prime-cube
subdomain an elementary proof suffices: eliminating x gives
(p−7)(p³−1)=0, so p=7 and x=19. No assertion is made without the prime-power
guard. The polynomial identity at composite base18 is not a divisor-sum
identity: sigma(18²)=847, not343.

## 4. A symmetry must declare which operations it preserves

On the squarefree cube of primes(3,5,19), sigma-squareclass fibers admit
16 arbitrary observational symmetries and8 affine ones. Requiring commutation
with every named prime toggle leaves only2: identity and the uniform3-toggle.
For a general rank-r binary observer on m bits, these named-update symmetries
are exactly translations by its(m−r)-dimensional kernel.

The same attachment can be free for one observer and source-revealing for
another. For3∤n,

    sigma(3n)=4sigma(n),       n=s(3n)−4s(n).

Therefore any representation answering both exact aliquot probes must retain
enough information to reconstruct n. The examples5 and19, both with aliquot
sum1, separate after attachment to s(15)=9 and s(57)=23.
[Operational symmetry groups and exact reconstruction](operations/OPERATIONS_AND_OBSERVERS.md).

## 5. Relations become objects in the tournament recursion

Quotienting a color matrix Z by GL_r retains its entire binary rowspace W,
not just its rank. Quotienting by affine palette transformations retains
span(1,rows Z). Marking a port retains(W,w), where w is the actual port-answer
word, with an essential rank guard when w=0. These are exact classifications
of the declared palette quotients. They are not permutations of arithmetic
primes preserving all their named values.

This turns colored tournaments into tournaments decorated with binary codes.
The joint Burnside formula uses Fix(g^l) for each vertex cycle; division by
the palette-group order is wrong. Source deletion retains the source-color
stabilizer, not the whole palette group. These give exact all-order formulas
and independently checked counts:

| Palette quotient | a_0,…,a_4 | h_1,h_2,h_3 |
|---|---|---|
|GL3(2)|1,2,5,24,188|1,3,19|
|marked19 port|1,3,11,70,756|2,8,59|
|affine palette|1,1,2,8,48|0,1,6|

The GL three-vertex count is16 transitive decorations plus8 cyclic ones;
hence h_3=24−5=19. These are typed class counts and harmonic dimensions,
not a proof that this19 is the arithmetic prime19.

The linear palette quotient retains a1/n color-count relaxation mode. The
affine quotient removes those first modes and has gap2/n for n>=2; its
first observable counts equal-color pairs. This accelerates the **quotient**,
not the named process or a number-theoretic trajectory. Equal-rank states
already separate under canonical source deletion, with explicit two- and
three-vertex controls.
[Joint quotient, spectrum, source stabilizers and relation objects](tournaments/PALETTE_SYMMETRY.md).

## Verification and scope

- [Arithmetic/code audit](audit_arithmetic_codes/AUDIT.md): primary citation,
  all256 selections, exact coordinate groups, parity transport, Phi6 uniqueness.
- [Operation/palette audit](audit_operations/AUDIT.md): independently generated
  relation-space orbits through degree4, source counts, refresh matrices,
  group hierarchy and coprime reconstruction probes.

Run `python3 verify_packet.py` for isolated normal/optimized replays and
manifest checks. SymPy1.14.0 is required by the palette operator checker;
other checkers use the standard library. See `verification.json` for exact
completed scope. Proofs establish infinite statements; finite tests guard
implementations. Shared primes, named ports, character zeros and the declared
update operations remain part of the carrier.
