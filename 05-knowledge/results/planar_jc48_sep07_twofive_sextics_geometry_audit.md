# Independent geometry and marked-group audit of the (4,6) cusp class

**Status: INDEPENDENT ANALYTIC + SOURCE AUDIT PASS for the scopes below;
FINITE-EXACT replay PASS.** Auditor: `three_ray_geometry`, 2026-09-07.
This audits the normal-form and local-resolution parts and the complete
infinity-(2,13) marked-group argument in
[the primary proof](planar_jc48_sep07_twofive_sextics.md), together with
its exact source. The openness, connected-family and simultaneous
resolution transfer are assigned to a separate referee; this note does
not substitute a local-geometry audit for that global argument.

## Normal forms and local types

I checked the finite cusp normalization, including the exceptional chart.
Translations of parameter and target, a nonzero scaling of U, and
subtracting the fourth-degree coefficient times U from V are lawful
polynomial affine operations. They give
\(U=t^4+at^3+bt^2\), \(V=ct^6+et^5+ft^3+gt^2\), with \(c\ne0\).

For \(b\ne0\), cancelling the quadratic term forces
\(f=ga/b\) if the singularity is to have type (2,5), rather than (2,3).
Putting \(d=g/b\), the local shear
\(V-dU+(d/b^2)U^2\) has fifth coefficient \(e+2ad/b\); its
nonvanishing is the exact (2,5) condition. The degree-raising local shear
is correctly kept out of the degree-preserving normal-form operations.
For \(b=0,a\ne0\), a nonzero quadratic coefficient of V yields type
(2,3), while a zero one makes the multiplicity at least three. For
\(a=b=0\), one instead uses V as the order-two local coordinate;
\(V^2-g^2U\) has fifth coefficient \(2fg\), so the exact criterion is
\(fg\ne0\). These branches exhaust the declared finite cusp.

I independently checked the infinity coefficient logic and the
normalizing scalings. If the seventh coefficient vanishes with
\(b\ne0\), then \(e=3ca/2\); \(a=0\) would make both polynomials
even and contradict birationality. Substituting \(t=a\tau\) and
dividing U by \(a^4\), V by \(ca^6/2\), gives the stated normalized
\(\beta=b/a^2\), \(\delta=2d/(ca^2)\) family. The successive even
translations and odd coefficients in the source agree exactly.
The exceptional \(U=t^4\) case gives infinity seven when \(e\ne0\)
and infinity nine when \(e=0\), because \(f\ne0\).

The arithmetic-genus calculation is correctly typed: birationality and
the basepoint-free degree-six homogenization give a rational sextic.
The finite cusp contributes delta two, an infinity (2,m) branch
contributes \((m-1)/2\), and each ordinary node contributes one.
Thus \(N=(17-m)/2\); in the stated \(N\ge2\) class, \(m\le13\)
already follows from genus. No assertion about all higher-infinity
coefficient specializations is needed for this exhaustion.

## Why zero even coefficients do not change the marked tree

This is the load-bearing local point for transporting the infinity-nine
plumbing result across its normal forms. At the unique infinity point,
\(X\) has order two and \(Z\) has order six. After the first three
blowups write
\(\rho=Z/X^3-c^2\). The branch meets \(\rho=0\) on the new
exceptional component. The original line meets at \(\rho=-c^2\ne0\),
and the previous exceptional component is at infinity in this chart.
They do not pass through the subsequent branch centres.

Suppose the present chart is \((X,Y)\), its exceptional component is
\(X=0\), and the next even term is \(kX\). Blowing up the branch point
and using \((X,Y/X-k)\) removes it. The old exceptional component is
at infinity in the second coordinate. This remains true if \(k=0\):
zero is an ordinary finite point of the new exceptional component, not
its intersection with the old one. If k vanishes, the order of Y may
increase from two to three or more, but the branch still has multiplicity
two and the blowup centre remains free. Thus a zero even coefficient
changes an expansion order, not this proximity relation.

Eventually the chart is \((X,Y)\) with orders (2,1). The curve is then
smooth and tangent to the exceptional component \(X=0\). Blow up that
free point to get \((X/Y,Y)\), with orders (1,1). The two exceptional
axes and the strict curve pass through the origin, with distinct
curve tangent from both axes. Blowing up their intersection gives the
final exceptional component. The strict curve meets it at the nonzero
finite coordinate \(X/Y^2\) and is transverse because Y has order one.
This final free/satellite pair is uniform.

Consequently the infinity-seven, nine, eleven and thirteen branches
have respectively zero, one, two and three free double-point centres
after the initial three blowups, followed by this final pair. This is
exactly the marked-tree sequence used in the proof.

Two exact hostile controls make the zero-coefficient issue concrete:

* In the normalized family, \(\beta=-1/4,\delta=0\) gives \(k_4=0\)
  and \([z]\tau=4\ne0\). The post-third-blowup chart has orders (2,3),
  followed by (2,1), then the same final free/satellite pair. The
  infinity-nine tree remains the audited one.
* For \(U=t^4,V=2t^6+t^3+t^2\), also \(k_4=0\). Here
  \(X=z^2/(2+z^3+z^4)\),
  \(\rho=(2+z^3+z^4)^2-4\), so \(\rho\) has order three and
  \(\rho/X=8z+O(z^2)\). The exceptional normal form gives precisely
  the same infinity-nine tree.

The identical argument covers \(k_5=0\) in the infinity-eleven family.
For example \(\beta=-1/3,\delta=-3/8\) gives \(k_4=2,k_5=0\),
with \([z]\eta=29/12\ne0\). These are local tree controls; no assertion
about a particular control's entire finite singularity inventory is
needed here.

## The actual infinity-thirteen graph

At \(\beta=3/20,\delta=-11/10\), the normalized coefficients are
\(k_4=-48/5,k_5=87/5,k_6=-7069/250\). With
\(\phi=\eta/X-k_6\), one has
\(\phi=(3/100)z+O(z^2)\) and \(X=(1/2)z^2+O(z^3)\).
The final transverse coordinate is therefore
\(X/\phi^2=5000/9+O(z)\), away from the exceptional intersections.
The seven intermediate chart transitions give exactly the eight centres
stated in the primary proof, and the curve arrow is on E8.

Reconstructing the squares from L of square +1 gives, in the order
\((L,E_1,\ldots,E_8)\),
\((-2,-2,-2,-2,-2,-2,-3,-2,-1)\). The edges are
\(L-E_3-E_2-E_1\) and \(E_3-E_4-E_5-E_6-E_8-E_7\).
The marked abelianization also reconstructs independently: start with
\(L=-6\mu\) and add incident valuations plus the strict curve
multiplicities \(2,2,2,2,2,2,1,1\) at the eight centres. This gives
\((-6,-4,-8,-12,-10,-8,-6,-5,-10)\mu\), exactly as stated.
It agrees with the negative intersection matrix of determinant -1.
The fibre relations therefore have the same positive complex normal
orientation and arrow sign as the already audited predecessor.

As an extra nonvacuity control, the normalized infinity-thirteen pair
actually has just the required finite cusp and two ordinary nodes.
Its derivative gcd is t. Symmetric pair elimination gives
\[
 M=-\frac{p^2(2p+1)(20p^2+40p+21)}{80},
 \qquad q=\frac{p(p^2+p+3/20)}{2p+1}.
\]
The apparent factor \(2p+1\) is not a collision: the first pair
equation is nonzero there. The remaining quadratic gives two pairs.
An ordinary ordered-pair resultant is
\[
-\frac{t^4}{160000}
 (1200t^4+2400t^3+1580t^2+420t+63).
\]
The tangent-determinant ideal has Groebner basis
\(s+400t^3/9+20t^2/3+t,\ t^4\), supported only at the cusp diagonal.
The three-fibre remainder test gives the unit ideal. Thus neither a
hidden tangency nor a shared third branch is being used as a false
positive example. This control is supplementary to the conditional
class proof, which assumes the finite inventory.

## Analytic A6 elimination

I checked the unsimplified vertex equations and the eliminations in
(5). The commutator \([c,d]=1\) justifies the powers of c in the
expressions for g,h,f; it is not an unmentioned commutativity assumption.
The common square/cube c is identity, a double transposition, or a
five-cycle. The double-transposition branch gives an order-four f,
which cannot be a square in \(A_6\).

For \(c=1\), \(f=d^7\). The order-one, two, four and five cases
respectively force a meridian of order at most two, order four, an
impossible square f, or a meridian of order five. If d has order three,
its square root in \(A_6\) is uniquely \(e=d^2\): an alternative
would have order six, which \(A_6\) does not contain. This holds for
both three-cycle types before the meridian condition is applied.
Then \(\mu=d^2\), so the required single three-cycle meridian forces
d to be a single three-cycle. The group \(\langle l,a\rangle\) is a
quotient of the (2,3,3) triangle group; all other fibres lie in it.
A faithful transitive six-letter \(A_4\) action has order-two
stabilizers and no fixed point for an order-three element. The smaller
quotient and orbit cases give support at most four when a single
three-cycle is present. The image therefore fixes labels.

For five-cycle c, \(d=1,h=c^3,f=1\), and \(e\) is an even
involution. Let j be the fixed label of h. If the single three-cycle
\(\mu\) moves j, write it as \((j\ q)(j\ p)\), with p,q distinct
on the five-cycle. Multiplying h first by \((j\ p)\) joins its
five-cycle and singleton into one six-cycle. Multiplication by
\((j\ q)\) then splits this into exactly two cycles. Thus \(\mu h\)
has two cycles on all six labels. But \(\mu h=e^{-1}\) is an even
involution, which has four or six cycles, a contradiction. Therefore
\(\mu\) fixes j, and the equation forces e and all fibres to fix j.
This proves the stated marked obstruction without a census.

## Independent finite controls and frozen replay

I reconstructed the census using **unsimplified successive vertex
elimination**: enumerate all raw pairs (l,a) in \(A_6^2\), retain
\(l^2=a^3\), then compute
\[
 b=a^2,\quad d=(lb)^{-1}c^2,\quad g=c^{-1}d^2,\quad
 h=d^{-1}g^2,\quad f=g^{-1}h^3.
\]
For every square root e of f, recover \(\mu=(he)^{-1}f\) and test
all original vertex equations and crossing commutators. This does not
assume the simplified h or f formulas during the reconstruction.
It independently gives **1,120 assignments**, with 40,360,720 having
respectively three,two,one fixed labels; 400 have c=1 and 720 have c a
five-cycle. A separately generated positive boundary image has order
60 and fixes exactly one label, so abelianity would be an invalid upgrade.

I also checked every one of the **2,880** pairs consisting of a five-cycle
h and a single three-cycle moving its fixed label. The product \(\mu h\)
has types (4,2), (3,3), (5,1), with counts 1440,720,720. Every case has
exactly two cycles, independently confirming the cycle-cut argument.

The entire producer source was read. Its universe and retained filters
are complete, negative exponents are lawful because the exponent of
\(A_6\) divides 60, and every represented candidate is checked against
the original relations before counting. Independent normal and optimized
runs are both byte-for-byte equal to the saved output:
**3,385 always-active gates, 1,890 bytes**.

```text
python3 -B 04-computation/planar_jc48_sep07_twofive_sextics.py
python3 -B -O 04-computation/planar_jc48_sep07_twofive_sextics.py
source SHA256 02ef96f5a39fee7382fa81d0970fd06285c7bb3e180362bff234d4e2abad1e75
output SHA256 2713e10ba1715d4ae3b9d5d0095be91fa901ea3f63a16c9d53778e15b072e472
assignment SHA256 b35120fec36b3ded751aec092972d2e870466fb4d7defe0a193d41cc8412373f
```

No mathematical correction is requested in these audited parts. The
actual infinity-thirteen representation obstruction, the same-tree
infinity-nine transfer across zero even coefficients, and the exact
source all pass. The separate referee must still supply the assigned
connected-family/topological-transfer audit before this note can count
as an audit of the whole broader-class theorem.
