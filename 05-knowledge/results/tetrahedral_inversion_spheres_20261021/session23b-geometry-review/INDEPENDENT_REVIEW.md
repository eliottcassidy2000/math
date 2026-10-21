# Independent review: seven tetrahedron-derived points

## Verdict and scope

PASS. I independently derived the obstruction, exceptional-locus classification,
center, radius, and coincidence structure before reading the producer's report.
I then reviewed `session23b-geometry/geometry_report.md` and its exact checker.
The mathematical claims in that report are correct with the stated assumptions:
a genuinely noncoplanar tetrahedron, nonzero centroid relative to the original
circumcenter, a face circumcenter different from the original circumcenter, and
A defined as the second intersection of the line through V0 and the centroid
with the original circumsphere.

The universal seven-point sphere assertion is false even when every point is
finite. Its exceptional solutions have exactly five distinct points, and those
five points are not coplanar, hence are not concyclic. The review makes no
historical novelty claim and does not validate a Fano-plane or tournament
identification merely from the presence of seven labels.

## Independent derivation

Normalize the original sphere to center 0 and radius 1. Let v be V0, let s be
the sum of the other three vertices, let w be their circumcenter, and put
u=v+s, q=|u|^2, k=|w|^2. The face plane is x·w=k, so s·w=3k. The construction is

    G=u/4, C=v, D=u, E=4u/q, B=w/k, F=4w, H=s−2w.

The centroid is strictly inside the original sphere. In particular A and C
are distinct, G lies strictly between them, and 0≤q<16. The case q=0 is
excluded because E is not finite; k=0 is excluded because the face tangent
planes have no common finite point.

If D and E are distinct, a hypothetical common sphere would have contradictory
powers at G:

    Power along AC = |G|^2−1 < 0,
    Power along DE = (D−G)·(E−G) = 3(1−|G|^2) > 0.

The second equation cannot be applied when D=E. Thus every exceptional
solution requires q=4, precisely D=E.

Write a candidate sphere as |x|^2−2z·x+c=0. Since A,C are on the original
sphere, z is perpendicular to A−C and hence to s−3v. Subtracting the C and D
equations when q=4 gives z·s=3/2. Consequently z·v=1/2 and c=0. The new sphere
passes through the original circumcenter.

The B equation gives z·w=1/2. The F equation then requires k=1/4, precisely
B=F. Finally the H equation requires |s|^2=3. Together with q=4 this yields
v·s=0. The center constraints are exactly

    z·v=1/2, z·s=3/2, z·w=1/2.

Let p=v·w. The Gram matrix of v,s,w has determinant

    det Gram(v,s,w) = 3(1/16−p^2).

The Gram matrix is positive semidefinite, so −1/4≤p≤1/4. At p=1/4, v=4w−s
lies in the opposite face plane, contrary to tetrahedron nondegeneracy. At
p=−1/4, v=s−4w, and the center constraints force both z·v=1/2 and z·v=−1/2.
The negative endpoint is a real, nondegenerate tetrahedron, so it cannot be
discarded merely by assuming a tetrahedron. This is why linear independence,
or equivalently the strict inequality p>−1/4 under the other assumptions,
is an essential part of the theorem.

For −1/4<p<1/4, the unique solution is

    z=(v+4p s+4w)/(2(1+4p)),
    radius^2=|z|^2=(5+12p)/(4(1+4p)).

Substitution verifies all seven labels, so the conditions are sufficient as
well as necessary. This matches the producer's formula exactly.

## Degeneracy and coincidence audit

1. A=C cannot occur in a genuine tetrahedron. The centroid is strictly inside
   the sphere and therefore the line from a vertex through the centroid is
   a secant, not a tangent.
2. G=0 makes E undefined as a finite Euclidean inversion point. This occurs,
   for example, in a regular tetrahedron and must not be silently assigned O.
3. w=0 makes the three face tangent-plane equations inconsistent, not merely
   nonunique. Their normals have rank two, but affine independence of the
   three face vertices makes the augmented system rank three.
4. The negative rank endpoint p=−1/4 satisfies all three scalar equalities
   and can have nonzero tetrahedral volume, but has no finite common sphere.
5. The positive rank endpoint p=1/4 is a coplanar original configuration. If
   such inputs were allowed, it would admit nonunique spheres and additional
   coincidences. It is correctly outside the theorem's assumptions.
6. On the valid locus, A=(s−v)/2. The only coincidences among the seven labels
   are D=E and B=F. A,C,H have norm 1; B,D have norm 2. A=H or B=D would force
   v=4w−s, the excluded positive endpoint; C=H would put V0 in the face plane.
   The remaining equalities are ruled out by these norms or |s|^2=3.
7. H has norm 1 and lies on the face circle. A nondegenerate triangle's
   orthocenter lies on its circumcircle exactly when the triangle is right,
   in which case H is the right-angle vertex. Thus the geometric description
   of the exceptional face is correct.

## Stronger exact nonconcyclicity certificate

The producer's radical-plane proof is valid. There is also a direct
noncoplanarity identity:

    det(A−C, H−C, D−C) = −3 det(v,s,w) ≠ 0.

It follows by substituting A=(s−v)/2, C=v, H=s−2w, and D=v+s. Thus four of
the five distinct points are affinely independent. This establishes both
nonconcyclicity and uniqueness of the sphere without relying on a drawing.

## Exact controls and results

The independent `exact_controls.py` imports no producer code and uses
explicit checks that remain active under Python's optimization flag.
It verifies:

- Rational counterexample V0=(0,0,−1), V1=e1, V2=e2, V3=e3. The sphere
  coefficient and augmented ranks for A,C,D,E are 3 and 4, respectively;
  the planar circle determinant is −56/9. H,B,F are three distinct collinear
  points, disproving the universal assertion independently of A's wording.
- The stated algebraic exceptional example. Every residual vanishes, the
  sphere ranks are 4 and 4, the center is (0,1/2,1), radius squared is 5/4,
  and precisely the pairs D=E and B=F coincide.
- The negative endpoint using the same face and V0=(sqrt(3)/2,0,−1/2).
  Its tetrahedral determinant is 3/2, all three scalar equalities hold,
  and its sphere system nevertheless has ranks 3 and 4.
- The positive endpoint using V0=(−sqrt(3)/2,0,1/2). Its tetrahedral
  determinant is zero and sphere system has ranks 3 and 3. Two different
  sphere centers are exhibited to check the nonuniqueness outside scope.
- G=0 and w=0 controls, the symbolic Gram determinant, center and radius,
  residual factorization, and noncoplanarity identity.
- The valid nearby inversion theorem A,C,O,E concyclic on the rational
  counterexample, with center (1,1,−1/2) and radius squared 9/4 in x=y.

Normal and optimized runs are recorded separately. The outputs contain
SHA-256 hashes pinning the producer artifacts reviewed. No Git operations
were performed as part of this independent review.

## Additional audit: the generic six-point inversion sphere

The producer subsequently added a positive construction, which I separately
checked. If v,s,w are linearly independent, the sphere with constant term 4
and center z determined by

    2z·v=5, 2z·s=q−1, 2z·w=1+4k

contains B,C,D,E,F and 4C. Its squared radius is |z|^2−4≥9/4, so it is a
genuine real sphere. The cross-product reciprocal-basis expression in the
producer's report is the exact solution to these three scalar equations.

For any sphere through D,E, substitution yields

    Q(E) = (q−4)(c−4)/q.

Thus D≠E forces c=4; likewise B≠F forces c=4. Independence then gives the
claimed uniqueness. When both pairs coincide, only three distinct labels
remain among B,C,D,E,F, so the stated uniqueness caveat is necessary. The
constant-4 sphere is invariant under inversion x↦4x/|x|^2, since

    Q(4x/|x|^2) = 4Q(x)/|x|^2.

It is consequently orthogonal to the radius-2 sphere. The inverse pairs are
(D,E), (B,F), and (C,4C). Neither A nor H lies on this sphere: their exact
residuals are μ(16−q)>0 and 6−2v·s>0, respectively. These are strict for every
nondegenerate tetrahedron in the independence domain.

The six-point construction has constant term 4; the exceptional seven-label
construction has constant term 0. They must not be conflated, including on
the exceptional locus where both constructions can coexist.

The independent exact controls now also solve for the generic center directly
from the linear equations and compare it with the cross-product expression,
verify all six point residuals, verify both strict-outside residuals and the
inversion identity symbolically, and reproduce the producer's rational core
example with center (−43/24,37/6,−5/2), radius squared 25049/576, and residuals
Q(A)=1178/225, Q(H)=38/5. This additional theorem passes the independent audit.
