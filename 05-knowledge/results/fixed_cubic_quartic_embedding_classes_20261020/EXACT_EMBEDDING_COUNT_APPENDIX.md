# Exact one-versus-two embedding count

20 October 2026. **PROVED; independently reviewed PASS.**

See `review/AUDIT.md` for the independent proof and exact-control audit.

This bounded appendix refines, without changing, the at-most-two theorem in
`FINITE_EMBEDDING_CLASSES_AND_FIXED_CYLINDERS.md`. It uses its already
proved generation criterion and the previously audited triple-point
inequivalence witness.

## The exact count

Let

    A=C[h,k] ⊂ C[y],
    h=y³+a y,   k=y⁴+b y²+c y,   a,b,c∈C.

The number of plane embeddings of Spec A up to ambient polynomial
equivalence is

    two, if c≠0 and 2a=3b;
    one, otherwise.

On the two-class locus the representatives are

    (h,k) and (k,q),  q=h²+(3c/2)h.

The main theorem already proves that these are the only possible classes
and that the second pair generates precisely on that locus. It remains
only to prove inequivalence there.

## 1. Every intrinsic automorphism extends in the standard embedding

An automorphism of A extends to its normalization C[y], so has the form
y↦αy+β, with α≠0. This preserves polynomial degree. The filtration identities

    A_≤3=C⊕C h,
    A_≤4=C⊕C h⊕C k

therefore imply

    h↦λh+e,
    k↦μk+νh+f,

where λμ≠0. The same formulas define an affine automorphism of the ambient
plane. Consequently, any ambient equivalence of the two embedded images
could be adjusted by an ambient lift of the induced intrinsic automorphism
to give the exact parametrized lift

    Φ(h,k)=(k,q).                                (1)

It suffices to rule out (1). The normalization maps themselves are not
embeddings of A¹ when they identify distinct parameters; the embedded
object throughout is Spec A.

## 2. The singular-point determinant obstruction when a≠0

Assume b=2a/3 and c≠0. For each root s of

    L(s)=s³+(4a/3)s−c,                           (2)

let r,t be the roots of

    z²−sz+(s²+a)=0.

They satisfy h(r)=h(t) and k(r)=k(t). Their squared separation is

    (r−t)²=−3s²−4a.                             (3)

It is nonzero: if it vanished, then s²=−4a/3 and (2) would reduce to
−c=0. Equivalently the resultant of (2) and the right side of (3) is
−27c². Thus every s in (2) yields a genuine two-parameter collision.
Its common h-value is

    h(r)=h(t)=−s(s²+a)=−c+(a/3)s.                (4)

The cubic L has at least two distinct roots. Otherwise it would be
(s−s₀)³; its zero quadratic coefficient would force s₀=0 and then c=0.
When a≠0, equation (4) therefore gives at least two distinct singular
image points with different h-values. They are singular because a smooth
point of the plane curve is normal and has only one normalization preimage.

Let F(U,V) be an irreducible equation for the standard embedded curve.
If an ambient automorphism Φ satisfied (1), then its coordinate differences
from

    T(U,V)=(V,U²+(3c/2)U)

would be multiples of F. At a singular point both F and its two first
partial derivatives vanish, so DΦ=DT there. Consequently the constant
Jacobian determinant of Φ would have to equal

    det DT=−2U−3c/2=c/2−(2a/3)s                (5)

at each of the singular points furnished above. These values differ for
distinct s when a≠0. This contradiction rules out (1), and Section 1
upgrades it to inequivalence of the embedded images.

## 3. The scaled triple-point case when a=0

Now b=0 and c≠0. Choose d∈C* with d³=c and set y=dt. Then

    h=c t³,
    k=cd(t⁴+t),
    q=c²(t⁶+(3/2)t³).

An affine scaling of the plane coordinates transforms the two pairs into

    (t³,t⁴+t),
    (t⁴+t,t⁶+(3/2)t³).

The second coordinate of the latter is affinely related by
z=8(t⁶+(3/2)t³)+3 to the independently audited triple-point witness in the
[boundary packet](../cubic_quartic_globalization_embedding_boundaries_20261020/CUBIC_QUARTIC_GLOBALIZATION_AND_EMBEDDING_BOUNDARIES.md).
That witness proves unparametrized ambient inequivalence by incompatible
cotangent determinants and collision periods. Affine source reparametrization
and invertible affine target scalings preserve equivalence or inequivalence.
Hence the two classes remain distinct for every c≠0.

This completes the exact count. It strengthens the description of fixed
curves only. It supplies no additional uniform lifting or moving-family
obstruction beyond the main fixed-cylinder corollary.

## Exact controls

Run from this packet directory:

    python3 verify_exact_embedding_count.py
    python3 -O verify_exact_embedding_count.py

These symbolic controls verify the collision equations, nonzero resultant,
common h-value, singular-point determinant formula, triple-root coefficient
criterion, affine filtration map, and triple-point scaling. They complement,
but do not replace, the universal argument above.
