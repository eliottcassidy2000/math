# Seven-point tetrahedron sphere claim: independent audit

## Verdict

The proposed universal theorem is false. This is not just an undefined-point issue: there are rational, nondegenerate tetrahedra for which every named point is finite and three of the named points are distinct and collinear. A complete classification is possible: with the natural second-intersection interpretation of A, the seven labels can be cospherical only on a special locus where two pairs of labels coincide. There are then exactly five distinct points, not seven. They are not all concyclic.

Everything below is an original calculation for this audit. No historical novelty or literature priority is claimed.

## 1. Normalization and formulas

Translate and scale so the tetrahedron's circumsphere is the unit sphere centered at 0. Write

- v = V0;
- s = V1 + V2 + V3;
- t = v + s = 4G;
- w = O'; k = w·w;
- T = t·t; p = v·w.

The face plane is w·x = k, and s·w = 3k. If k > 0, the tangent-plane intersection is B = w/k. The remaining points are

C = v, D = t, E = 4t/T, F = 4w, H' = s − 2w.

Let q = t − 4v. The second intersection is

A = v − 2(v·q)q/(q·q).

Every vertex has norm 1, and tetrahedral nondegeneracy implies ||G|| < 1, hence T < 16. In fact v·G < 1, so q is nonzero and the second intersection A is always distinct from v. Tangency at v is not an exception for a genuine nondegenerate tetrahedron: G is interior to the circumsphere.

If G = 0, E is undefined as a finite Euclidean point. If w = 0, the opposite face passes through O and the three tangent planes have no common finite point; B is undefined. These are genuine flaws in an assertion made for every tetrahedron, but the counterexample below avoids both.

## 2. An exact rational counterexample

Take

V0 = (0,0,−1), V1 = (1,0,0), V2 = (0,1,0), V3 = (0,0,1).

The circumsphere is the unit sphere, and the tetrahedron is nondegenerate. Direct calculation gives

G = (1/4,1/4,0), O' = (1/3,1/3,1/3),
A = (4/9,4/9,7/9),
B = (1,1,1),
C = (0,0,−1),
H' = (1/3,1/3,1/3),
E = (2,2,0),
D = (1,1,0),
F = (4/3,4/3,4/3).

The three distinct points H', B, F are collinear. A Euclidean sphere intersects a line in at most two points. This alone disproves the claim, even if A's wording is interpreted differently.

There is also an independent four-point certificate. A,C,D,E lie in the plane x=y, and the determinant with rows (x,z,x²+y²+z²,1), in that order, is −56/9, not zero. They are not concyclic in that plane and therefore not cospherical in three dimensions.

The same three-collinear-point obstruction works for any equilateral opposite face with 0 < ||O'||² != 1/4: then H' = O', while B = O'/||O'||² and F = 4O' are distinct points on one line.

## 3. A dimension-free four-point obstruction

Suppose all four points A,C,D,E lie on a sphere S, and D != E. Because G lies between A and C on their chord, G lies strictly inside S. But D=4G and E=G/||G||² lie on the same ray from G, and a ray from an interior point meets a sphere only once. This is impossible.

Equivalently, the directed power of G computed along AC is

Pow_S(G) = ||G||² − 1 < 0.

Computed along DE, it is

Pow_S(G) = (D−G)·(E−G) = 3(1−||G||²) > 0.

One must retain D != E: repeated labels do not give two intersection points and cannot be used in a secant power formula.

For an algebraic certificate, write the sphere equation as

Q(x) = ||x||² − 2u·x + c = 0.

Since A,C are on the unit sphere and their line contains G,

2u·t = 4(2u·v), and 2u·v = 1+c.

The equation Q(D)=0 yields c=(T−4)/3. Then

Q(E) = (T−4)(T−16)/(3T).

For a nondegenerate tetrahedron with E defined, 0<T<16. Therefore Q(E)=0 forces T=4. This is exactly the coincidence D=E.

## 4. Complete exceptional-locus theorem

Assume the tetrahedron is nondegenerate, G != O, O' != O, and A is the second intersection. Then all seven labels lie on a finite sphere if and only if

(i) ||t||² = 4;
(ii) ||w||² = 1/4;
(iii) ||s||² = 3;
(iv) v,s,w are linearly independent.

Equivalently, conditions (i)–(iii) hold and p=v·w > −1/4. Nondegeneracy already excludes p=1/4.

When these conditions hold, the unique sphere has center

u = (v + 4p s + 4w) / (2(1+4p))

and radius squared

rho² = ||u||² = (5+12p) / (4(1+4p)).

It also passes through O. In unnormalized coordinates, apply the same formulas to (Vi−O)/R and (O'−O)/R, then multiply u and rho by R and translate the center by O.

### Necessity

Section 3 forces T=4 and c=0. The equations at C and D become

2u·v=1, 2u·s=3.

The equation at B=w/k gives 2u·w=1. The equation at F=4w then gives 16k−4=0, so k=1/4. Finally,

||H'||² = ||s−2w||² = ||s||²−8k,
2u·H' = 2u·s−4u·w = 1.

Consequently Q(H')=0 forces ||s||²=3. Since T=||v+s||²=4, these conditions imply v·s=0.

To settle the rank issue, let r=w−s/4. Then r·s=0 and ||r||²=1/16. Since v·s=0 and ||v||=1,

|p| = |v·r| <= 1/4.

At p=1/4, v=4w−s belongs to the opposite face plane, so the tetrahedron is degenerate. At p=−1/4, v=s−4w. The required equations would then give

2u·v = 2u·s − 4(2u·w) = 3−4 = −1,

contradicting 2u·v=1. Hence p is strictly between −1/4 and 1/4, and v,s,w are independent.

### Sufficiency and explicit center

For the displayed u, use

||v||²=1, ||s||²=3, ||w||²=1/4,
v·s=0, v·w=p, s·w=3/4.

Direct substitution gives

2u·v=1, 2u·s=3, 2u·w=1.

These equations establish Q=0 at all seven labels, with c=0. Independence gives uniqueness. Expanding the numerator of ||u||² gives

5+32p+48p² = (1+4p)(5+12p),

which proves the radius formula.

### Geometric interpretation and repetitions

Conditions (ii)–(iii) give ||H'||²=1. Since H' lies in the face plane, it lies on that triangle's circumcircle. For a nondegenerate triangle this happens exactly when the triangle is right-angled, with H' equal to the right-angle vertex. Thus an equivalent description is:

- The opposite face is a right triangle.
- Its plane is at distance R/2 from O.
- The tetrahedron's centroid is at distance R/2 from O.
- The rank-collapse endpoint p=−1/4 is excluded.

Always D=E and B=F. The three distinct points A,C,H' lie on the original circumsphere; B and D have norm 2 and are distinct from those three and from each other. Therefore precisely five distinct points remain.

For completeness, A != H': otherwise G,v,H' would be collinear, and because s=H'+2w, the vectors v,s,w would be dependent. Also B=D would imply v=4w−s, the excluded coplanar case.

These five points are not all concyclic. The three points A,C,H' define the radical plane of the original unit sphere and the new sphere, namely 2u·x=1. But B and D satisfy 2u·x=4, not 1. Thus no nondegenerate tetrahedron with all definitions finite satisfies the stronger literal reading that all seven labels lie on a single circle.

## 5. Exact nontrivial exceptional example

Take

V0 = (0,1,0),
V1 = (sqrt(3)/2,0,1/2),
V2 = (0,sqrt(3)/2,1/2),
V3 = (0,−sqrt(3)/2,1/2).

Then

A = (sqrt(3)/4,−1/2,3/4),
C = (0,1,0),
H' = V1 = (sqrt(3)/2,0,1/2),
B = F = (0,0,2),
D = E = (sqrt(3)/2,1,3/2).

They lie on the sphere with center (0,1/2,1) and radius sqrt(5)/2. Every displayed sphere residual is exactly zero in the attached verification script.

A convenient one-parameter family keeps V1,V2,V3 fixed and, for nonzero real a, sets

V0 = (−sqrt(3)(1−a²)/(2(1+a²)), 2a/(1+a²), (1−a²)/(2(1+a²))).

The common sphere has center (0,a/2,1) and radius sqrt(1+a²/4). The limit a=0 is a coplanar tetrahedron; the limit |a|→infinity reaches the rank-collapse configuration with no finite common sphere.

## 6. A corrected generic five-point theorem, naturally six-point

There is a true generic inversion construction inside the false claim. Suppose delta = v·(s cross w) != 0. Define

u = [5(s cross w) + (T−1)(w cross v) + (1+4k)(v cross s)] / (2 delta).

Then the sphere

Q(x) = ||x||² − 2u·x + 4 = 0

contains B,C,D,E,F. Its radius squared is ||u||²−4, which is at least 9/4 because u·v=5/2 and ||v||=1. If D!=E or B!=F, this is the unique sphere through those five labels. If both pairs coincide, it is the unique sphere with constant term 4 through the three remaining distinct labels, but other spheres through those three points exist.

Indeed the reciprocal-basis formula gives

2u·v=5, 2u·s=T−1, 2u·w=1+4k.

Substitute the five point formulas to verify Q=0. The two claimed extra points are not merely generically absent; they are strictly outside this sphere:

Q(H') = 6−2v·s > 0,
Q(A) = mu(16−T) > 0,

where A=v+mu(t−4v) and mu=−2v·(t−4v)/||t−4v||²>0. The first inequality uses v·Vi<1 for each other vertex, and the second uses T<16.

This sphere is invariant under inversion of radius 2 centered at O, because

Q(4x/||x||²) = 4Q(x)/||x||².

That inversion exchanges D with E, B with F, and C with 4C. Therefore 4C is an additional natural sixth point on the same sphere. The construction is a generic six-point sphere formed by three inverse pairs, not the claimed seven-point theorem. The sphere is orthogonal to the radius-2 sphere since ||u||²=rho²+4.

An exact generic example is

V0=(0,0,−1), V1=(1,0,0), V2=(0,1,0), V3=(3/5,0,4/5).

Here delta=−4/15, T=18/5, k=4/9. The common five-point/six-point sphere has center (−43/24,37/6,−5/2) and radius squared 25049/576. Exact residuals are 0 at B,C,D,E,F,4C, while Q(A)=1178/225 and Q(H')=38/5.

The independence condition matters. For example, the simple equilateral-face counterexample in Section 2 has s parallel to w, and even B,C,D,E,F are not cospherical there.

## 7. Generalization and a valid related inversion theorem

For an n-dimensional simplex, n>=2, let m=n+1 and D=mG, while A,C,E retain their meanings. The same ray/power argument proves that no hypersphere contains A,C,D,E unless D=E, which is equivalent to ||G||²=1/m after normalization. The obstruction is dimension-free: all four points already lie in the plane spanned by O,G,C.

More generally, with any interior nonzero point g in place of the centroid and D=lambda*g for lambda>1, the obstruction is identical unless D=E. Algebraically it also excludes lambda!=0 whenever the two radial points are distinct. The lambda>1 version has the simplest ray explanation.

There is a clean, valid nearby four-point theorem: if O,G,C are not collinear, then A,C,O,E are concyclic. Inversion in the original circumsphere fixes A,C, sends E to G, and takes a circle through O to a line. The required line is ACG. This is a related corrected construction, not evidence of which correction the original author intended. If O,G,C are collinear, A,C,O are three distinct collinear points and the finite-circle statement fails.

No adjustment of terminology from 'concyclic' to 'cospherical' repairs the universal seven-point assertion.

## 8. Implication for combinatorial analogies

Seven names are not seven geometric points: whenever the claimed cosphericity actually holds, D=E and B=F. A Fano-plane analogy therefore cannot literally identify these labels with its seven distinct points. Any useful analogy must instead specify an incidence structure, a field or characteristic, or an explicitly formal role for the seven labels, independently of this false universal geometric assertion.

## 9. Verification artifact

The companion `verify_geometry.py` uses exact SymPy expressions, checks the rational counterexample, checks the classification formulas symbolically, checks a one-parameter exceptional family, checks the generic six-point inversion theorem and strict outside residuals, and audits random rational unit-sphere tetrahedra. Its output is saved as `verification.txt`.
