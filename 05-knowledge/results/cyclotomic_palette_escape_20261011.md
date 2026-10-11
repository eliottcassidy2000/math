# Cyclotomic residue palettes, and a unit that escapes a local colouring

2026-10-11. **PROVED, scoped and self-audited:** C1--C4 below.
**FINITE-EXACT:** finite-field spectra and colour maps. **INHERITED:**
residue colourings, THM-4558 Lemma R. **OPEN:** the chromatic number of the
full field Q(zeta_55), a forcing six-colour context, and the plane target.
No lower bound is transferred backwards through reduction.

## Inheritance and the decisive hostile

[THM-4558, residue colourings of field planes](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md)
already supplies the known local colouring mechanism and stronger bounds
than the initial dyadic tests in our [quartic comparison](pentagon_context_number_fields_20261011.md).
The Moser field has chromatic4; the Polymath field
Q(sqrt(-3),sqrt(-11),sqrt(-15)), and its Cartesian enlargement
Q(sqrt3,sqrt5,sqrt11)^2, have chromatic5. Adjoining the golden quadratic
there is an inherited successful increase, not a new conjecture here.

Incoming commit08590cec09 was read during this session. Its
[all-length field-escape family](unit_field_escape_recursion_20261011f.md)
is a stronger hostile: for every finite prime bank it supplies actual
four-chromatic geometric graphs whose fields escape all corresponding
fixed-prime tests. Its [defect profiles](unit_defect_frontier_20261011f.md)
retain complete boundary support. Evading a test alone is insufficient.

This note's niche is a small explicit cyclotomic compositum, its exact
first available residue size, and an infinite unit-direction generator
certifying how reduction fails. Board: **unit norm, conjugation, prime
decomposition, retained boundary partition, full incidence, colour witness**.

## C1. A residue graph can need six colours while its source needs three

Let p be one of3,5,11,13. Two generates (Z/pZ)^*, so reduction of Q(zeta_p)
at2 gives F_(2^(p-1)). Put q=2^((p-1)/2). Conjugation reduces to z->z^q.
Unit-vector residues lie in the norm-one group of order q+1. Its additive
Cayley graph has these exact certificates:

| p | q | Vertices | Degree | Minimum eigenvalue | Chromatic bounds |
|---:|---:|---:|---:|---:|---|
|3|2|4|3|-1|4|
|5|4|16|5|-3|4|
|11|32|1024|33|-11|4|
|13|64|4096|65|-15|6 through8|

These bounds concern the **finite residue graphs**. Every prime cyclotomic
field Q(zeta_p), p odd, has chromatic number3: at its unique prime above p,
conjugation fixes the prime and acts trivially on F_p. Unit vectors are
local units reducing to plus or minus1. The residue graph is C_p, whose
three-colouring pulls back by Lemma R. Conversely, partial sums of
1,zeta_p,...,zeta_p^(p-1) give a regular odd unit p-cycle.

Thus **Q(zeta_13) has chromatic3 although its dyadic residue graph needs
at least6**. A hard finite quotient need not lift to a hard metric graph.

For the certificates, label field elements by binary polynomials modulo
Phi_p. The eigenvalue at additive character a is
lambda_a=sum_(u^(q+1)=1) (-1)^bitdot(a,u). The replay independently computes
all character sums and the Walsh transform, and checks two spectral moments.
For p=11 the complete spectrum, as(value,multiplicity), is

    (-11,33),(-7,165),(-3,330),(1,165),(5,165),(9,165),(33,1).

For p=13 it is

    (-15,195),(-11,520),(-7,585),(-3,390),(1,780),
    (5,780),(9,455),(13,390),(65,1).

For a d-regular graph on N vertices with minimum eigenvalue -s, an
independent-set indicator of density rho yields
0>=d rho^2 N-s rho(1-rho)N by splitting off its constant component.
Hence its size is at most Ns/(d+s). This gives256 at p=11 and768 at
p=13, proving the respective four- and six-colour lower bounds. At p=5,
all65536 subsets instead certify the sharper independence number5.

Colour a field element by its binary dot products with these masks:

    p3: (1,2); p5: (1,14); p11: (511,654);
    p13: (1,726,3448).

For every norm-one difference at least one dot product is1, checked exactly.
These prove the upper bounds4,4,4,8. The last exact chromatic number is
not determined here.

## C2. Fifth and eleventh roots have incompatible dyadic involutions

The completions at2 of Q(zeta_5) and Q(zeta_11) have degrees4 and10,
with conjugation the Frobenius powers2 and5 respectively. In their
compositum L=Q(zeta_55), a fixed prime at2 would require simultaneously

    a=2 mod4,       a=5 mod10.                         (1)

Parity makes this impossible. Equivalently ord_55(2)=20 but no power
of2 is -1 modulo55. The two primes over2, each of residue degree20,
are exchanged by conjugation. This only defeats the specified local test.

The failure has an explicit unit witness. Set

    f(X)=X^20+X^18+X^17+X^15+X^14+X^12+X^9
          +X^7+X^6+X^5+X^2+X+1,
    g(X)=X^20 f(X^-1).

Modulo2, Phi_55=f g and gcd(f,g)=1. All irreducible factors of Phi_55
have degree ord_55(2)=20, so these are exactly its two irreducible factors.
For zeta=zeta_55 put

    eta=f(zeta)/f(zeta^-1)=zeta^20 f(zeta)/g(zeta).       (2)

Neither factor vanishes, since its degree is below40.
Then eta conjugate(eta)=1. At the f prime eta has positive valuation
and reduces to **zero**; at the other it has a pole. The real unit edge
from0 to eta collapses under reduction at the first prime. The exact
missing hypothesis is conjugation fixing that prime.

The resultant gives

    Norm(f(zeta))=295823077277696=2^20 *991*284681.

Since the f prime has residue degree20 and g does not divide f modulo2,
eta's two dyadic valuations are exactly1 and -1. Thus eta is not an
algebraic integer, contrasting with the integral Salem generators.

## C3. The first available local norm graph has q = 109

At a prime ell not dividing55, the decomposition group is generated by
ell modulo55. Conjugation fixes a prime iff -1 belongs to that group.
Modulo5 an element of order4 requires an even exponent to reach -1;
modulo11 any element reaching -1 requires an odd exponent, since its
order is2 or10. The exact condition is therefore

    ell=-1 mod5, and ell is a quadratic nonresidue mod11. (3)

The order is2 or10, with conjugation its half power. The local norm graph
uses

    q=ell if ell=-1 mod11, and q=ell^5 otherwise.        (4)

At5 the unramified factor has Frobenius5 modulo11 of odd order5, so no
prime is conjugation-fixed. At11 its unramified factor has Frobenius1
modulo5, again excluding a fixed prime.

The only primes through109 satisfying(3) are19,29,79,109. Their q values
are2476099,20511149,3077056399,109. This list is exactly checked; every
larger prime has q>=ell>109. Consequently **109 is the minimum q over
all conjugation-fixed finite primes of L**. None of the small-residue
colourings in the inherited table applies directly.

This does not rule out other colouring methods or establish high chromatic
number. Incoming four-chromatic field-escape examples are explicit hostiles
to that inference. Compared with their arbitrary finite-prime-bank theorem,
this note gives a compact degree40 example and an exact all-place minimum,
not a stronger universal escape theorem.

## C4. A recurring star, and the context it still needs

The nonzero valuation makes eta a non-root of unity, so its powers are
distinct unit directions. If |eta^i-eta^j|=1 for i!=j, eta^(i-j) would
have real part1/2 and unit norm, hence be a primitive sixth root of unity.
That would make eta itself a root of unity. Therefore

    {0,eta^a1,...,eta^ar},       all a_i distinct,

is an induced unit star for every r. Its exact extension weight at k
colours is k minus the number of colours used on its leaves, as in the
inherited star gadget. This is a different arithmetic realization of an
existing relation, not a chromatic increase.

The private-leaf compiler can use these powers around a golden pentagon
in L. At any stage coincidences and unwanted unit distances to other
existing vertices exclude finitely many powers: two distinct unit circles
meet in at most two points. Choose another power, record its exponent,
and authenticate every pair. This produces controlled geometry, with
the twenty leaf colours still unconstrained.

The [pentagon context theorem](pentagon_context_number_fields_20261011.md)
requires a six-chromatic external context under its stated interface.
Incoming marked defect antichains give the exact composition state:
retain every equality partition of the twenty ports, compose compatible
profiles, and test whether all surviving proper-colouring profiles use
the same four colours in all five leaf groups. A proper six-colouring of
the context must also be retained, preventing a vacuous implication.
A spectrum, field degree, or total extension count loses this information.

**OPEN helpers.** Can a mixed-field overlap change this full boundary
support while preserving its metric realization? Can integral directions
in a larger field do so, or do opposite local valuations offer a useful
construction coordinate? Can the antichain compiler expose a small exact
obstruction before expanding the point set? The induced stars are positive
controls, not evidence of forcing.

[Replay](../../04-computation/experiments/marked_context_transfers_20261011.py)
and [output](marked_context_transfers_20261011.json) retain norm-one
differences, spectra, masks, reciprocal factors, resultant and prime list.
The all-prime conclusion follows from(3)--(4) and q>=ell, not an unbounded
computer search.
