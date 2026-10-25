# Session 60: independently audited upper bound for the outer six-line K3

## Exact result

Let X/Q be the minimal resolution of the double plane

    w² = E(E²−x²)(E²−z²)(E+x+z).

Then **rho(X over Qbar) <= 19**. This is unconditional: it uses good reduction at
7 and 13 and the proven Tate/Artin–Tate theorems, not a generic-rank assertion.
Together with an independently established rank-19 lower bound, it gives equality.
This packet does not itself reprove that lower bound or the identification of X
with the symplectic quotient of the outer multiquadratic K3.

## 1. Arrangement and a smooth integral model

The six line normals, in coordinates (E,x,z), are

    (1,0,0), (1,−1,0), (1,1,0), (1,0,−1), (1,0,1), (1,1,1).

There are precisely nine ordinary double points and two ordinary triple points.
The triples are [0:0:1] and [0:1:0]. The supplied checker enumerates every
pairwise intersection and every triple determinant. The nonzero determinants
have absolute value 1, 2, or 3. Thus the incidence pattern and transversality
persist over R=Z[1/6]. All centers below are disjoint rational sections.

Blow up the two triple points of P²_R, including the exceptional divisors in the
new branch divisor (odd branch multiplicity). There are now fifteen ordinary
nodes of the branch divisor: the original nine doubles and three on each new
exceptional divisor. Blow up these fifteen nodes, excluding their exceptional
divisors from the branch (even multiplicity). The resulting branch divisor is
a disjoint union of eight smooth rational curves. Its double cover is smooth
and proper over R, because 2 is invertible and the branch is relatively smooth.
This proves good reduction at every p>=5, including 7 and 13.

For completeness, on the base Y write the branch as 2L. Initially K+L=0.
At a blowup of branch multiplicity m=2 or 3 the new line bundle is
L'=pi*L−floor(m/2)E and K'=pi*K+E, so K'+L'=0 persists. The cover has trivial
canonical bundle. Since Y is rational, the splitting
f_*O_X=O_Y plus O_Y(−L)=O_Y plus O_Y(K_Y), and Serre duality, give h¹(O_X)=0.
The cover is geometrically connected (the branch is nonempty), hence its
fibers are K3 surfaces. The cover is the minimal resolution: all exceptional
curves described next have self-intersection −2.

## 2. Resolution correction and eighteen rational classes

Over an ordinary node the resolution has one P¹. Indeed, the preimage of the
base exceptional P¹ is a degree-two cover branched at its two rational
intersection points with the strict branch. It is a genus-zero curve with a
rational ramification point, so it is P¹ over the base field.

Over either triple point, the first exceptional divisor is itself a branch
component; its preimage is a central P¹. Resolving its three intersections with
the strict branch creates three more P¹'s by the preceding node construction.
Their graph is D4. All four components and all three intersection points are
rational over Q, and over every F_p for p>=5. In particular, no unidentified
quadratic splitting field is being assumed.

Consequently an A1 fiber contributes q extra points over its contracted point,
and a D4 tree contributes 4q extra points: its union has 4(q+1)−3=4q+1 points.
The correction from the singular double plane to X is exactly **17q**.

The pullback h of a general line has h²=2 and is orthogonal to all contracted
curves. Their Gram lattice is

    <2> direct-sum A1(−1)^9 direct-sum D4(−1)^2,

of rank 18 and determinant 2(−2)^9*4²=−16384. Thus these eighteen divisor
classes are independent and defined over Q and F_p. No assertion that nineteen
classes are rational, and no elliptic-fiber splitting assumption, is needed.

## 3. Exact point counts

For q odd, put

    S_q = sum_{x,z in F_q} chi((1−x²)(1−z²)(1+x+z)).

The chart E=1 gives q²+S_q points on the singular cover. On E=0 every point is
branched, giving q+1 more points. The resolution calculation proves

    #X(F_q) = q² + 18q + 1 + S_q.

The independent checker uses F_(p²)=F_p[u]/(u²−d), where d is the smallest
positive nonsquare greater than 1. For every element it checks the quadratic
character both by Euler exponentiation and by explicit square-set membership;
it counts every affine (x,z), and cross-checks the singular cover point total.
Normal and optimized Python runs have identical outputs.

| p | S_p | S_(p²) | #X(F_p) | #X(F_(p²)) |
|---|---:|---:|---:|---:|
| 5 | 2 | 4 | 118 | 1080 |
| 7 | 4 | 100 | 180 | 3384 |
| 11 | 14 | 196 | 334 | 17016 |
| 13 | 4 | 484 | 408 | 32088 |

These agree with the separately written counter of the other Session60 worker;
that implementation was not consulted in constructing this one.

## 4. Reconstructing Frobenius without guessing its sign

Use geometric Frobenius F on untwisted H². Its eigenvalues have absolute value
p. The eighteen classes split off a nondegenerate F-stable subspace with
characteristic polynomial (T−p)^18. On the orthogonal complement V of dimension
four, F/p is orthogonal. Let a=tr(F|V), A=tr(F²|V), and b=(a²−A)/2.

The characteristic polynomial on V is reciprocal or anti-reciprocal with
multiplier p². If det(F/p)=−1, the middle coefficient b must be zero. Explicitly
its polynomial would be T⁴−aT³+0T²+ap²T−p⁴. Thus nonzero b forces determinant
+1, and the polynomial is

    R_p(T)=T⁴−aT³+bT²−ap²T+p⁴.

The point counts give (a,A,b)=(4,100,−42) at 7 and (4,484,−234) at 13. Therefore

    charpoly(F at 7)  = (T−7)^20 (T²+10T+49),
    charpoly(F at 13) = (T−13)^20(T²+22T+169).

The normalized two remaining roots have respective rational sums −10/7 and
−22/13. If either root were a root of unity, this sum would be a rational
algebraic integer, hence an integer; it is not. These roots are therefore not
roots of unity. The geometric Picard ranks of both reductions are at most 20.
The counts at 5 and 11 have b=0, so we do NOT claim that their first two traces
alone determine the reciprocity sign or their characteristic polynomials.

## 5. Tate, Artin–Tate, and the discriminant obstruction

Charles's Corollary 2 proves the Tate conjecture for K3 surfaces over finite
fields of characteristic at least 5. Hence both reductions have geometric
Picard rank exactly 20. To avoid a hidden cohomological semisimplicity
assumption: the rational Neron–Severi lattice is nondegenerate, so it splits
off orthogonally in H². If eigenvalue p remained on its complement, that
complement would contain a p-eigenvector, contradicting the Tate invariant
statement. Thus the multiplicity of p is exactly the rank over F_p; applying
the same argument over extensions gives the geometric assertion. All algebraic
Frobenius eigenvalues are p, so Frobenius acts trivially on the rational
Neron–Severi space: its action on this lattice has finite order and therefore
is semisimple. The Picard group over F_p thus has the same rank and
discriminant square class as the geometric one.

[Charles, *The Tate conjecture for K3 surfaces over finite fields*, Cor. 2](https://arxiv.org/pdf/1206.4002).

Here is the exact normalization of Artin–Tate used. For a K3 surface over F_q,
let Phi(T) be the characteristic polynomial of F/q, and let rho be the rank
over F_q. Then

    |disc NS| * #Br(X) = q * lim_(T→1) Phi(T)/(T−1)^rho.

The Tate conjecture implies this formula; this implication is Milne's theorem.
For K3, chi(O_X)=2 and Pic^0=0, so the exponent of q is exactly one.
[Milne, *On a conjecture of Artin and Tate*, 1975, Theorem 4.1](https://math.stanford.edu/~conrad/BSDseminar/refs/MilneAT.pdf);
[author's corrections and clarification](https://www.jmilne.org/math/articles/1975a.html).
The displayed K3 normalization is also stated explicitly in
[Elsenhans–Jahnel, §4.2–4.3](https://wwwuser.gwdguser.de/~jjahnel/Arbeiten/pic_dst_cam2a.pdf).

The finite Brauer group has square order, by
[Liu–Lorenzini–Raynaud, Theorem 1](https://dinolorenzini.franklinresearch.uga.edu/sites/default/files/inline-files/papers/LLR2InventionesPublished.pdf).
Their [2018 corrigendum](https://arxiv.org/abs/1804.11158) corrects an intermediate
formula but expressly preserves this square-order conclusion.

For rank 20 the intersection form has signature (1,19), so its determinant is
negative. At 7 the remaining normalized polynomial is T²+(10/7)T+1, whose
value at 1 is 24/7. Hence the discriminant square class is −24, equivalently −6.
At 13 the value is 48/13, giving square class −48, equivalently −3.

If rho(X_Qbar)=20, good-reduction specialization embeds NS(X_Qbar) as a
finite-index sublattice in each rank-20 reduction. Since specialization preserves
intersection numbers, both reduction discriminants must have the same square
class as NS(X_Qbar). Their ratio 2 is not a rational square, a contradiction.
Thus **rho(X_Qbar)<=19**.

This is the usual two-reduction discriminant argument; see
[van Luijk, *K3 surfaces with Picard number one and infinitely many rational points*, §2–3](https://arxiv.org/pdf/math/0506416).

## Boundaries

* This proof is specific to the displayed fixed arrangement, not a genericity claim.
* The upper bound is complete; the exact-rank claim needs the separately proved
  nineteenth class. No conditional Tate or unverified resolution hypothesis remains.
* The smooth-model proof is geometric; the Python checker checks its finite
  incidence and arithmetic inputs, not a formal proof in a proof assistant.
* Only the square classes of the reduction discriminants are determined here.
