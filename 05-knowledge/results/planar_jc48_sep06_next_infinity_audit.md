# Independent audit: the next infinity branch and its marked survivor

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS; FINITE-EXACT replay
PASS.** Auditor: `three_ray_geometry`, 2026-09-06. The audited objects are
[the proof note](planar_jc48_sep06_next_infinity.md), its standalone source,
and its saved output. The positive result is an actual curve and a
surjective representation of its marked infinity-link group. Its extension
to the full affine complement, and any Keller realization, remain OPEN.

## Geometry and inventory

The curve is the actual polynomial image
\[
 U=t^4+t^3+t^2,\qquad V=16t^6+24t^5-19t^3-19t^2.
\]
I independently checked the parameter-selection calculation: the
coefficient of \(t^3\) in
\(V-cU^{3/2}-hU\), where \(h=-3a^2c/8-3bc/2\), is
\(a(7a^2c+12bc+16d)/16\) in the declared family. The formal fractional
power selects a parameter; it is not used as a polynomial coordinate
change or as the proof of the geometry.

The exact finite shear gives orders (2,5) and leading fifth-order
coefficient -14. The derivative gcd is t, so there is no other ramified
normalization point, and gcd(U,V)=t² excludes another preimage of the
cusp image. This is a single finite (2,5) cusp.

I checked the divided-difference elimination and its denominator:
\(p=s+t,q=st\) satisfy
\((2p+1)q=p(p^2+p+1)\), and \(p=-1/2\) is impossible.
The second divided difference becomes
\(-p^2(4p^3+10p^2+15p+14)\). The cubic is separable and coprime to the
pair-discriminant numerator and the denominator; it yields three
non-diagonal pairs. The tangent-determinant ideal is supported only at
the cusp diagonal, so all these pairs are transverse. A three-point
fibre would force the displayed fixed-degree cubic remainder to divide
\(U-a\); the three resulting coefficients generate the unit ideal.
Thus no two pair collisions share a third branch. The finite singularity
inventory is exactly the cusp and three ordinary nodes.

Finiteness of the normalization map follows already from monicity of U;
the finite collision list implies generic injectivity. Its projective
parametrization has degree six with no base point at infinity, so the
image is a rational sextic. This supplies the types needed for the genus
check, without assuming birationality from the parameter degrees alone.

At infinity the exact coordinate cancellations give
\[
 Z-256X^3+15360X^4-753664X^5
       =-\frac{17}{1024}z^{11}+O(z^{12}),\qquad \operatorname{ord}X=2.
\]
The nonzero first odd order is eleven, hence the branch has type (2,11).
I checked all eight chart order pairs and the final nonzero translation.
The original infinity line leaves the branch at the third blowup,
\(\rho=-256\), whereas the branch passes through \(\rho=0\).
The final parameter \(\eta=-17408z+O(z^2)\) is transverse; the arrow
point \(X/\eta^2=1/4848615424\) lies away from the other divisor
intersections. The simultaneous line/curve resolution therefore uses
finite multiplicities (2,2,1,1) and infinity multiplicities
(2,2,2,2,2,1,1). Their delta sum plus the three nodes is 10, the sextic
arithmetic genus, while \(D^2=4\) and \(D^2-2N=-2\). No strict Nori
conclusion is available.

## Independent ordinary-resultant controls

As a check using coordinates different from the producer's symmetric
pair elimination, I independently reconstructed
\(F(u,v)=\operatorname{Res}_t(U-u,V-v)\). It is irreducible, has total
degree six and is monic of vertical degree four. Its exact discriminant is
\[
\operatorname{disc}_vF=-u^5(256u^2+11u+12)
 (1183744u^3+1870400u^2+974463u+168070)^2.
\tag{A}
\]
For the ordinary two-parameter divided differences N(s,t), M(s,t),
elimination of s gives
\[
\operatorname{Res}_s(N,M)=t^4(
1088t^6+2720t^5+4136t^4+2968t^3+1827t^2+686t+686).
\tag{B}
\]
These are independent exact controls, not a numerical inference of
transversality or of the singularity inventory. The source's critical,
tangent and three-point tests still pay those separate obligations.
Equations (A)--(B) also retain concrete inputs for a future certified
braid computation, which has not been performed for this curve.

A direct reconstruction of these controls is:

```python
import sympy as S
s,t,u,v=S.symbols('s t u v')
U=t**4+t**3+t**2
V=16*t**6+24*t**5-19*t**3-19*t**2
F=S.resultant(U-u,V-v,t)
print(S.factor(S.discriminant(F,v)))
N=S.cancel((U.subs(t,s)-U)/(s-t))
M=S.cancel((V.subs(t,s)-V)/(s-t))
print(S.factor(S.resultant(N,M,s)))
```

## Marked boundary and group witness

The seven centres are exactly
\(L;L\cap E_1;L\cap E_2;E_3;E_4;E_5;E_5\cap E_6\).
I independently reconstructed the tree and squares. Compared with the
preceding equality curve, the extra free point of multiplicity two
lengthens the central chain before the final free/satellite pair. The
weights in the order \((L,E_1,\ldots,E_7)\) are
\((-2,-2,-2,-2,-2,-3,-2,-1)\), with the arrow on \(E_7\).

The positive-complex-fibre convention of the audited predecessor gives
the eight displayed vertex relations and all crossing commutators.
Because the graph is a tree of rational components, the circle-bundle
presentations have no genus generators or graph stable letters. Gluing
them gives the complete marked boundary presentation, so a solution of
all these relations does define a representation of the actual infinity
exterior. This sufficiency is needed for the present positive witness;
checking only an arbitrary list of necessary relations would not suffice.

The marked abelianization is
\((-6,-4,-8,-12,-10,-8,-7,-14)\mu\). I checked it independently by
starting with \(L=-6\mu\) and adding the incident valuations plus strict
curve multiplicities \(2,2,2,2,2,1,1\) under the seven blowups. It matches
the negative intersection matrix and fixes the arrow sign. In particular
this remains an actual infinity calculation, not the projective quotient
after killing the line meridian.

I independently reconstructed the permutation witness from l,a,e using
composition \((pq)(i)=p(q(i))\): set \(b=a^2\),
\(d=(lb)^{-1}\), \(g=d^2\), \(\mu=(ge)^{-1}\), and \(c=f=1\).
The resulting d,g,mu agree with the displayed raw tuples; d has order
five and mu is the single cycle (012). All vertex equations, adjacent
commutators and parity conditions hold.

An independent breadth-first closure, multiplying on the left and
including inverse generators, gives:

| Generators | Order | Orbit of 0 | Common fixed labels |
|---|---:|---|---|
| l,a | 60 | all six | none |
| g,e | 60 | 0,1,2,3,4 | 5 |
| l,a,e | 360 | all six | none |

The last group consists of even permutations and has order 360, hence is
\(A_6\). The former tetrahedral support argument has a precise failure:
when \(c=f=1\), the new relations allow \(d^5=1\), and the first
order-60 piece has a faithful transitive six-letter action. Sharing the
five-cycle with the second piece does not force a common fixed label.
The witness therefore refutes extending that support obstruction to this
new boundary. It does not supply a representation of the affine group:
the latter is a quotient with potentially additional relations.

## Frozen source and scope

I read every source check. The universe is one exact polynomial curve
and one raw labelled representation. It is not a census, and no
unrecorded search count is needed by the conclusion. The two symbolic
ideals, chart tests, graph construction, all original relations and
commutators, parity and generated orders are always-active checks.

Independent normal and optimized replays both equal the saved output
byte-for-byte: **63 gates, 1,565 bytes**.

```text
python3 -B 04-computation/planar_jc48_sep06_next_infinity.py
python3 -B -O 04-computation/planar_jc48_sep06_next_infinity.py
source SHA256 4e39a2fecef78b7b312d4cbec3606bdef3d9d5fbd25f3630d3ce6c4c461cf4e2
output SHA256 4b6d18e523c5163f538a84c645395cf37cfe9539aec1c069a743f4803bf3fa35
semantic SHA256 ed764a10ed4c3bc549a97b7b7f15d76d32aa262c584406e69f29bf7dd33853dc
```

The complete actual geometry, boundary presentation and finite witness
pass. The stated open boundary is correct: neither a full affine
monodromy representation nor a finite algebraic cover nor a Keller map
has been constructed. No correction is required.
