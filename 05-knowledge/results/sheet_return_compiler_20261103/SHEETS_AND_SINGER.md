# From two converse sheets to an actual three-valued scale phase

**PROVED scoped maps and finite field certificates.** Graph-cover and Singer
machinery are classical and inherited. The new bridge is the explicit map
from the complemented tournament's ordered-two-bin response to the old chart.

## 1. Recovering a discarded 343/57 operation, not only its numbers

Session35's [SINGER_LINE_FLIPS.md](../session35_dimensions_dualities_20261023/arithmetic/SINGER_LINE_FLIPS.md)
uses `F343=F7[t]/(t^3+t^2+1)` with primitive element `g=t+5t^2`.
Its trace-zero projective line in Z57 is

    D={8,15,19,38,48,50,51,56}.

The current response recurrence has root

    eta=4+4t+t^2=g^146,
    eta^3-2eta^2-eta-2=0, eta^57=2.

The other roots are g^314 and g^338, its Frobenius conjugates. This is an
explicit field isomorphism/choice of phase between the two charts.

For a joint state s=(S_n,S_(n+1),S_(n+2)), define

    Psi(s)=(eta^2-2eta-1)S_n+(eta-2)S_(n+1)+S_(n+2). (1)

Its coefficient transformation is triangular with unit diagonal, so it is
an F7-linear isomorphism. The actual recurrence gives `Psi(next s)=eta Psi(s)`.
Put `gamma=1/(3eta^2-4eta-1)=4+2t+4t^2`. Then

    S_n=Tr(gamma Psi(s)).                             (2)

One proof is to test the basis1,eta,eta^2: division by the derivative of the
minimal polynomial followed by trace extracts the eta^2 coefficient. The
checker verifies all three basis identities directly in the retained old
field. Consequently this map preserves both the growth operation and the
specific predicate `S_n=0 mod7`.

For A at n=4, Psi=(4,6,1) and the projective exponent of gamma*Psi is12;
for B they are (6,4,1) and15. Therefore, for r=n-4,

    S_n=0 mod7 iff a+146r mod57 belongs to D,
    with a=12 for A, a=15 for B.                      (3)

This classifies *every* n, not a finite sample:

| Family | Zero n residues modulo57 |
|---|---|
|A|6,7,11,13,23,41,44,55|
|B|2,4,14,32,35,46,54,55|

There are eight because a projective line of PG(2,7) has eight points.
Multiplication by eta induces a regular57-cycle on its projective points:
eta has order171 and its intersection with F7* has order3. The quotient
therefore has order57. This is the old Singer action, now driven by a literal
tournament attachment response. It is not a claim that the tournament has
57 vertices, 343 automorphisms, or that its vertex set is this projective plane.

The older aliquot edge343->57 has the same prime-cube/geometric-sum formula
`7^3 -> 1+7+7^2`. Here that sum is the norm exponent and number of projective
points. The valid transferred operation is (1)--(3); no aliquot iteration is
identified with tournament growth. The old source1001 and its uniqueness
theorem stay in [Session82](../session82_decimal_aliquot_square_sum_20261028/README.md).

## 2. Projectivization loses an observable trit

After57 bulk steps the S state returns projectively but is multiplied by2.
After171 it returns literally modulo7. Each projective state thus has three
actual scale states on this orbit. Dropping the scale preserves zero tests
but changes nonzero counts. For example the first A response has residue6;
57 steps later it is5, and another57 later it is3.

The H triple has literal period48; its sixteenth power is4I modulo7.
The joint projective clock has period `lcm(16,57)=912`. After912 steps, H
has returned literally and S has been multiplied by2. The full joint phase
returns after2736. Hence there is a literal map of controlled cycles

    Z/2736 -> Z/912,     r -> r mod912,

intertwining one bulk insertion with addition by1. Its fiber has size3, and
one trip around the base cycle acts nontrivially on that fiber. The missing
coordinate is a trit, not a free choice of a tournament isomorphism class.

This directly parallels the *form* of
[Session166's reflection-voltage classification](../session166_planar_metagraph_20261103/reflection_voltages/REFLECTION_VOLTAGES.md):
after a spanning tree fixes local labels, a closed walk can retain transport
that no local relabeling removes. There the coefficients are C2 and the object
is the actual black-line quotient with its individual incidences. Here the
coefficients are C3 and the object is this finite operation clock. In both,
retaining quotient vertices without the transport around cycles loses data.
The two covers are not identified and their realization languages differ.

## 3. The switch and converse retain different data

If two labelled tournaments differ in one arc e and p is Hamiltonian in both,
p cannot use e. Their native complements along p still differ in exactly e.
Applying the same newborn mask preserves this commuting square on literal
relations. This justifies the common-path search in README; equality of the
old golden observers does not force equality after this operation.

For converse there is a separate exact covariance:

    native(T^op, reverse(p)) = native(T,p)^op,
    birth(T^op, complement(mask)) = birth(T,mask)^op.

Thus a sheet change must reverse the selected path and complement the newborn
mask together. Rooted starts and ends also exchange. Neither field Frobenius
nor the neutral source switch can silently stand in for this transport.

The open construction problem is now sharper: a return must match the literal
source relation and ports while converting both response cubics to the target
source's dyadic law, with the ternary scale and F49 phase aligned. Sheet signs
alone do not supply that conversion. State-dependent path changes are the
next promising operation; the constant two escape templates have been settled.
