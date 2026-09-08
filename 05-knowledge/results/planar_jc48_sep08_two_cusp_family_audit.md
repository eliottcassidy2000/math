# Independent audit of the two-ordinary-cusp family and its two cyclic strata

**Status: INDEPENDENT ANALYTIC + EXACT SOURCE + REPLAY AUDIT PASS.**
The [primary](planar_jc48_sep08_two_cusp_family.md) exhausts the declared
birational (4,6) class into infinity types7,9,11, and proves cyclic
affine complements for the entire good type7 and type9 strata. It
retains the nonempty type11 stratum as OPEN. Its optional numerical
words are not dependencies of this proof.

## 1. Complete parameters and actual normalization

I read the complete proof and source. Sending the two distinct cusp
preimages to ±1 and normalizing the two leading coefficients uses only
affine changes. Integrating their common derivative factor gives all
four initial parameters. The literal target shear V-3bU/2 removes b
and changes c to c+bs, preserving both the degree pair and the actual
affine complement. No normalization by a possibly zero s±1 is used.

The gcd condition P(s)!=0 makes the common derivative exactly t²-1.
The two determinant conditions distinguish ordinary cusps from higher
ones, including s=±1: there V supplies the quadratic coordinate even
when U does not. A quadratic intermediate function field would give an
affine involution fixing the unique pole at infinity. Both coordinates
would be polynomials in its squared centered parameter, forcing an
odd-degree derivative gcd. The degree-two gcd excludes this. The
remaining field degree is one, since it divides both4 and6. The monic
quartic makes the normalization finite; the homogeneous degree-six
map is basepoint-free. These checks justify the actual rational sextic
and its one infinity preimage before any good-locus inference.

## 2. Exhaustion, including the vanishing even coefficient

I checked the role of each source identity. The initial infinity
coordinates have orders2 and6, with nonzero leading coefficients. The
odd7 coefficient cuts out the type9 plane; the odd9 coefficient cuts
out the type11 graph. On that graph the odd11 coefficient is a nonzero
factor times P(s), so a higher infinity type violates the exact
critical-locus condition. This is a coefficient obstruction, not merely
a genus estimate. The genus equation then fixes the node counts5,4,3.

The even coefficients are subtracted as local coordinate changes,
without being inverted. Thus their possible vanishing does not remove
a member of a stratum. The three initial infinity blowups separate
the strict infinity line from the branch point Z3=1. Subsequent even
translations fix the exceptional axis and its marked intersections.
The final odd coefficient is a unit within each stratum, which pays
the constant marked resolution type used below.

## 3. Proper incidence and marked transport

I recovered the corrected proof in
[higher-odd classification §4](planar_jc48_sep07_higher_odd.md) and
[the (2,5) family transport](planar_jc48_sep07_twofive_sextics.md).
The current family supplies both cusp sections explicitly. Dividing
the two fibre equations by (t-e)² loses no extra preimage because the
ordinary-cusp determinant makes at least one divided value nonzero at
t=e. The first divided polynomial is monic quadratic, so the incidence
is finite over the parameters and its image closed. This also removes
coincident cusp targets. The exclusion precedes the normalization lift.

The fixed ordinary-cusp and infinity resolutions then really lift the
proper normalization: each center pulls back to a fixed parameter power
times a unit. The U'' and V'' charts cover the repeated-U' cases and
glue. After resolution the map is unramified and separated. Its diagonal
is open and closed, so the off-diagonal double incidence and distinct
triple incidence are proper. Their tangency and coincidence images are
closed. This pays openness of the complete good locus; no escaping
unresolved pair incidence is used.

Each nonempty good locus is open in one irreducible affine space, plane,
or line, hence path connected. A complex line through two good points
is not contained in the bad algebraic set and meets it at finitely many
points. Blowing up the entire finite étale node multisection avoids an
unpaid global labeling choice. The resulting proper relative normal
crossing pair admits a stratum-tangent flow along compact paths. All
centers were removed, so this identifies the actual affine complements,
with positive curve meridians preserved up to access conjugacy.

## 4. Openness across strata is used with the correct direction

The type9 supplier is the previously certified literal curve. Its four
configuration certificates impose finitely many strict inequalities.
The resultant coefficients vary polynomially with the three parameters;
therefore those very same rational loops and center polygons certify
the same words throughout a Euclidean neighborhood. The old common
base-fibre identification remains valid through its endpoint disks.

The exact type7 representative proves that its good locus is nonempty,
and therefore dense in the full parameter three-space. It meets that
certificate neighborhood, supplying one actual cyclic type7 member.
Connected marked transport then reaches every good type7 member. This
does not assert equisingularity across the type9/type7 boundary, and
does not reverse a specialization quotient.

No such inference is made for type11. It is a separate closed graph
away from the old supplier, has an explicit good three-node point, and
needs its own actual braid certificate. The arbitrary-group implication
of its five numerical words is useful only after those paths are paid.

The inherited winding argument makes each cyclic group infinite, with
a positive meridian generator. Its actual retained-sheet consumer then
excludes whole Keller support for the type7 and type9 classes in every
mapping degree. The abstract A6 passport and one-cusp Euler bounds are
not substituted for these global relations.

## 5. Exact reproduction

Independent normal and optimized runs match all 915 frozen bytes and
pass77 always-active gates. The two representative checks retain every
divided-difference branch, cusp diagonal, exceptional denominator,
tangency, and possible triple image. The constant-leading cubic test
excludes triple images because a quadratic remainder cannot vanish at
three distinct roots. The source recomputes each full resultant and
vertical discriminant. These controls support the representatives;
they do not extrapolate openness or topology from a finite sample.

    python3 -B 04-computation/planar_jc48_sep08_two_cusp_family.py
    python3 -B -O 04-computation/planar_jc48_sep08_two_cusp_family.py

Source SHA256: `50daef4d2689fe62ca23e66ba4fb9f8306c46afd50983ad83477321ab58eaba7`.

Output SHA256: `1ec6e9e83b3aac30adaf90e91ced31448c54c7a3dffa80d4c142acf313ba0332`.

The primary may be promoted with precisely its two-stratum exclusion
and three-stratum classification. JC(2) and the type11 group remain OPEN.
