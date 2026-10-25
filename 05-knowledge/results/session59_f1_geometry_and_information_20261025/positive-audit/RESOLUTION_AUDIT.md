# Independent second-wave audit: resolution and Frobenius recovery

**Verdict: PASS** for TheoremsR1–R3 and PropositionR4 in the geometric
`RESOLUTION.md`, read independently at2026-10-25 15:33UTC. This is separate
from the first-checkpoint `AUDIT.md`. No integer-square existence, additional
root-cover toricity, characteristic2 resolution, or broad F1 inference is
approved.

## 1. Direct local resolution

At each geometric node the local equation over its residue field L is
xy=z². The residue field contains i when needed; in the nonsplit case the
four geometric nodes comprise two degree-two closed points. Blow up the
maximal ideal(x,y,z). In the x-chart, write y=xY,z=xZ; the strict transform
is Y=Z² and is smooth. The y-chart is identical; the z-chart is XY=1.
These charts cover the blowup. Its exceptional conic in P²_L is XY=Z²,
which has an L-point and is explicitly P¹_L under[s:t]→[s²:t²:st]. Thus no
hidden nonsplit exceptional conic is present. Restricting the tautological
normal line bundle to this degree-two conic gives O_P¹(−2).

Blowing up the reduced singular subscheme is defined over the base field;
base extension to its residue fields gives exactly these local models.
Everything else was already smooth, so the result is a smooth resolution.
There are precisely four geometric exceptional curves, all−2; none is−1,
so this is the minimal resolution. This argument requires characteristic≠2,
as the draft states.

## 2. Independent global construction and map check

The parametrization

    (c,a,b)=(s²+t²,s²−t²+2st,s²−t²−2st)

satisfies a²+b²=2c². Its three quadratics have coefficient determinant−8,
a unit over every field of odd characteristic, so it is the standard
Veronese conic after an invertible linear change. In particular C≅P¹ over
the base field, with two infinity points conjugate precisely when−1 is
nonsquare.

For C×C→B the proposed five coordinates have common zero locus exactly
c1=c2=0. Around one such point, use x=c1/a and y=c2/d; a,b,d,e are units.
The base ideal is(x,y), so one ordinary point blowup resolves the rational
map. The map is an isomorphism on c1c2≠0, hence birational. This supplies a
proper smooth resolution independently of the node-blowup description.

To identify all contracted curves without relying on a fan assertion:

*if c1=0,c2≠0, the map is the fixed node[0:a:b:0:0]; the analogous claim
holds with the factors exchanged;
*each of the four such boundary fibers passes through two blown-up centers,
so its strict transform has self-intersection0−2=−2;
*at a blown-up center the exceptional P¹ of slopes[x:y] maps by the linear
coordinates[0:a0 y:b0 y:d0 x:e0 x] isomorphically to a boundary line of B;
these−1 curves are not contracted;
*every curve mapped into c=0 lies among these boundary components, while
outside the boundary the map was an isomorphism.

Therefore the four−2 strict fibers are exactly the exceptional locus and
this resolution is the same minimal resolution as in R1. This proves the
claimed identification rather than merely matching point counts.

The linear system has class2H1+2H2−ΣE_j=−K and square8−4=4. It is not a
(1,1) system after identifying the conics with P¹. The draft retains this
important degree distinction. The four blown-up centers and the four nodes
are different labelled sets, even though their Frobenius permutations have
the same two-transposition cycle type.

## 3. Counts by two independent geometric mechanisms

Let Q=q^n and ε=χ_q(−1). The singular surface has(Q+1)² points. Exactly
2(1+ε^n) nodes are rational over F_Q; each is replaced by P¹_F_Q, addingQ
points. Therefore

    #Y(F_Q)=Q²+(4+2ε^n)Q+1.

Alternatively C×C has(Q+1)² points. Its blowup-center scheme has
(1+ε^n)²=2(1+ε^n) rational points, each addingQ. This yields the same
formula from the second construction. Nonrational closed centers contribute
no rational points to either source or exceptional fiber. In particular,
adding4Q unconditionally would be wrong.

The zeta denominator follows by directly exponentiating this sequence:

    (1−T)(1−q²T)(1−qT)^4(1−εqT)^2.

For q=3 the resolution has16 points, but over F9 it has136; the distinction
appears at the extension level even when blowing up changes no F3 points.
No cohomological trace theorem is needed for this count or zeta derivation.

## 4. Independent Picard calculation

Over an algebraic closure the blowup construction gives the Picard basis
H1,H2,E_00,E_01,E_10,E_11. Both ruling classes are fixed by Galois. In the
nonsplit form, conjugation sends each sign pair(σ,τ) to(−σ,−τ), exchanging
the exceptional classes in two pairs. Thus the permutation representation
on the four E-classes is triv²⊕sign², and

    Pic(Y_kbar)⊗Q = triv⁴⊕sign².

This proves R4 by a route independent of the draft's toric-divisor sequence.
Its trace is2 in the nonsplit odd-extension case and6 after splitting,
matching4+2ε^n. There is no torsion or omitted ruling class in this count.

The fan computation was also checked: the eight cyclic rays alternate axis
and diagonal directions, adjacent determinants are1, self-intersections are
−2 on axes and−1 on diagonals, and central inversion pairs opposite rays.
The character sublattice is saturated of rank2, with Smith factors1,1;
the quotient has rank6. It removes precisely two sign directions from the
four-trivial/four-sign boundary permutation representation, agreeing with
the independent blowup basis.

This genuinely restores an arithmetic sidecar that the singular total count
had cancelled. It does not recover that sidecar from the old count alone:
the marked node scheme or explicit arithmetic form is essential.

## 5. Independent reproduction

`check_resolution.py` imports no producer implementation. It checks the
conic identity/determinant, fan determinants and self-intersections, the
two Picard models, anticanonical degree and contracted-fiber intersections,
direct conic/infinity enumeration for ten odd prime fields3 through31 as
listed in its JSON, six extension identities per field, and the first five
zeta power sums for each splitting sign. Both ordinary and optimized Python
runs must pass and produce identical output; the precise gate count is in
`resolution_certificate.json` and hashes in `resolution_manifest.json`.

No mathematical repair is requested. Preserve the odd-characteristic scope,
the distinction between the two labelled four-point sets, and the absence of
any assertion that rationality solves the additional square-root conditions.
