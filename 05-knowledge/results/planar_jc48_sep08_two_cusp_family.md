# Two ordinary cusps in the birational (4,6) class

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED; coefficient and
representative controls FINITE-EXACT.** The exclusion below is
only for infinity types \((2,7)\) and \((2,9)\). The genuine
\((2,11)\) stratum remains OPEN; its numerical scout is not a supplier.

## 1. Exact class and the three possible infinity strata

Consider an irreducible complex affine curve with a birational polynomial
normalization \((U(t),V(t))\), of degrees 4 and 6, whose finite
singularities are exactly two ordinary \((2,3)\) cusps and otherwise
ordinary nodes. The cusp images and their normalization preimages are
distinct. The classification has exactly these possible
infinity types and node counts:

| Infinity branch | Ordinary finite nodes | Affine complement consequence |
|---|---:|---|
| \((2,7)\) | 5 | PROVED: \(\pi_1=\mathbb Z\); excluded as whole Keller support |
| \((2,9)\) | 4 | PROVED: \(\pi_1=\mathbb Z\); excluded as whole Keller support |
| \((2,11)\) | 3 | Actual nonempty residual stratum; group and Keller exclusion OPEN |

Infinity type at least thirteen is impossible in this class. This is an
algebraic coefficient obstruction, not an extrapolation from the genus
bound. In particular all the curves in the declared class have at least
three nodes.

The closest proved supplier is the
[literal two-cusp cyclic-complement certificate](planar_jc48_sep08_two_cusp_braid.md)
and its [independent audit](planar_jc48_sep08_two_cusp_braid_audit.md).
The family transport follows the corrected proper-incidence and marked
resolution mechanism of
[the higher-odd classification](planar_jc48_sep07_higher_odd.md), §4,
with its hypotheses paid afresh for two cusp sections. The canonical
hostile is the local/Euler/three-generator \(A_6\) passport in
[the two-cusp passport](planar_jc48_sep08_two_cusp_passport.md): its
missing coordinate is the common global access gauge. The corrected
near miss in this family would be dividing by \(s-1\) or \(s+1\)
and discarding an ordinary cusp where U has order three. The least-used
sidecar is **openness of a finite strict rational braid certificate**,
which gives a nearby generic-stratum supplier without asserting infinity
isotopy across different types.

The five-concept board is parameter completeness, derivative common
zeros, infinity coefficients, connected good strata, and finite strict
configuration certificates. The passage from one stratum to another in
§5 uses only coefficient openness of four actual loops. The passage
within a stratum uses a proper marked resolution. These are distinct
maps with distinct hypotheses; neither is replaced by a bare limit.

## 2. Complete parameters and the repeated projection-critical cases

An affine parameter change sends the two cusp preimages to \(t=\pm1\).
Independent nonzero scalings of U and V make their leading coefficients
one, and target translations remove constants. Since both coordinate
derivatives vanish at both cusps, the complete forms are
\[
\begin{aligned}
U'&=4(t^2-1)(t-s),\\
V'&=6(t^2-1)(t^3+a t^2+b t+c),\\
U&=t^4-\frac{4s}{3}t^3-2t^2+4st,\\
V&=t^6+\frac{6a}{5}t^5+\frac32(b-1)t^4
       +2(c-a)t^3-3bt^2-6ct.
\end{aligned}\tag{1}
\]
The actual affine target shear \(V\mapsto V-3bU/2\) removes b and
replaces c by \(c+bs\). Thus it is enough to use the single irreducible
coefficient space \(\mathbb A^3_{s,a,c}\), with
\[
 U=t^4-\frac{4s}{3}t^3-2t^2+4st,\qquad
 V=t^6+\frac{6a}{5}t^5-\frac32t^4+2(c-a)t^3-6ct.
\tag{2}
\]
All operations in this normal form are genuine source affine changes
or target affine automorphisms. Subsequent local infinity subtractions
are singularity calculations, not claimed global target automorphisms.

Put \(P(t)=t^3+a t^2+c\). The derivative gcd is exactly
\(t^2-1\), up to a nonzero scalar, if and only if
\[
 P(s)\ne0.
\tag{3}
\]
For \(e=\pm1\), the ordinary-cusp determinant is
\[
 (U''V'''-V''U''')(e)
 =192D_e,\qquad D_e=(e-s)P'(e)-P(e).
\tag{4}
\]
Thus require \(D_+D_-\ne0\). These conditions include both
\(s=1\) and \(s=-1\): at \(s=e\), U'' vanishes but
\(D_e=-P(e)\ne0\), so V'' is nonzero and supplies the order-two
coordinate. For example \((s,a,c)=(1,0,1/6)\) has
\(U''(1)=0,V''(1)=14\), and determinant -224. This is a local
exceptional-coordinate control, not a claim about that point's remaining
global singularities.

The only repeated-root possibilities for U' are these double roots:
\(\operatorname{disc}_t(U')=1024(s^2-1)^2\). A triple root is
impossible because U' already contains the two distinct roots 1 and -1.
No division by \(s\pm1\) is used in the family argument.

Conditions (3)--(4) also provide birationality before any good-locus
argument. The degree of \(\mathbb C(t)/\mathbb C(U,V)\) divides both
4 and 6, so is one or two. In the degree-two case the nontrivial field
involution is a Möbius transformation fixing the unique pole of U, hence
an affine involution \(t\mapsto\beta-t\). Both U and V would then
be polynomials in \((t-\beta/2)^2\). Their derivative gcd would have
odd degree, of the form
\((t-\beta/2)R((t-\beta/2)^2)\), contradicting its degree two in
(3). Thus the degree is one. Finiteness follows from the monic quartic
equation for t over \(\mathbb C[U,V]\). The normalization is
\(\mathbb A^1\), and the basepoint-free degree-six projective
parametrization has rational sextic image with one infinity preimage.

The excluded point \((s,a,c)=(0,0,0)\) is a sharp control for this
step: both U and V are even and their derivative gcd has degree three.

## 3. Exact infinity exhaustion

At \(z=1/t\), set \(X=U/V,Z=1/V\). In all of (2),
\(X\sim z^2,Z\sim z^6\), with contact six along the original
infinity line. The first odd coefficient is
\[
 [z^7](Z-X^3)=\frac45(3a+5s).
\tag{5}
\]
Its nonvanishing gives infinity type \((2,7)\). On its zero plane,
\(a=-5s/3\), put
\[
 k_4=3-\frac43s^2.
\]
Then all earlier coefficients vanish after the indicated even removal,
and
\[
 [z^9](Z-X^3-k_4X^4)
 =4\left(c-\frac43s-\frac{14}{27}s^3\right).
\tag{6}
\]
Nonvanishing gives infinity \((2,9)\). On the remaining graph
\[
 a=-5s/3,\qquad c=4s/3+14s^3/27,
\tag{7}
\]
the next even removal and odd coefficient are
\[
\begin{aligned}
 k_5&=(128s^4-216s^2+1053)/108,\\
 [z^{11}](Z-X^3-k_4X^4-k_5X^5)
 &=-16s(s-3)^2(s+3)^2/81.
\end{aligned}\tag{8}
\]
But on this graph
\[
 P(s)=-4s(s-3)(s+3)/27.
\tag{9}
\]
Consequently (3) makes (8) nonzero, so infinity at least thirteen is
impossible. The ordinary cusp conditions on (7) are
\[
 D_e=-2(s-3e)^2(7s-3e)/27\quad(e=\pm1).
\tag{10}
\]
The starting infinity-eleven line is therefore
\(s\notin\{0,3,-3,3/7,-3/7\}\), before the further good-node
conditions. In particular it retains \(s=\pm1\).

The zero value of the intermediate even coefficient \(k_4\) at
\(s=\pm3/2\) is retained as well. It does not force an additional
infinity degeneration. After the third marked infinity blowup, the
subsequent even-coordinate changes fix the exceptional axis, as
explained below. The source checks the lower coefficients as well as
the displayed leading coefficients; no division by an even coefficient
is involved.

For the declared ordinary-node inventory, the rational-sextic genus
formula is
\[
 10=1+1+(m_\infty-1)/2+N.
\tag{11}
\]
This gives exactly the node counts in §1.

## 4. Intrinsic good loci are open and connected

The three starting spaces, subject to (3)--(4), are the following
irreducible open parameter spaces:

* \(T_7\): the open part of \(\mathbb A^3\) with \(3a+5s\ne0\).
* \(T_9\): the open part of the plane \(a=-5s/3\) with the
  coefficient in (6) nonzero; coordinates are \((s,c)\).
* \(T_{11}\): the open line (7), with the exclusions following (10).

Let \(G_m\subset T_m\) be the locus where the two cusp images are
distinct, have no extra preimages, and all other finite singularities
are ordinary nodes. We justify that each \(G_m\) is Zariski open,
using the corrected incidence mechanism of the audited supplier.

First exclude a second preimage of either cusp image. For each
\(e=\pm1\), the simultaneous equations
\[
 \frac{U(t)-U(e)}{(t-e)^2}=0,\qquad
 \frac{V(t)-V(e)}{(t-e)^2}=0
\tag{12}
\]
give exactly those extra preimages: at \(t=e\), at least one of
U'' and V'' is nonzero by (4), so no diagonal solution remains. The
first polynomial in (12) is monic of degree two in t. Thus this
incidence is finite over the parameter space, and its image is closed.
It includes coincidence of the two cusp targets by taking \(t=-e\).
The unique infinity preimage causes no extra finite cusp preimage.

After this exclusion, resolve both prescribed cusp sections and the
marked infinity branch in the proper \(\mathbb P^2\) family. The
normalization genuinely lifts: the pullback of the initial cusp centre
ideal is \((t-e)^2\) times an ideal containing a unit. Subsequent
ordinary-cusp centre ideals have the same fixed-order property because
(4) is a unit on the base. The charts with U'' or V'' nonzero cover and
glue, so the cases \(s=e\) require no new component or division.

At infinity the first three charts are
\(Z=XZ_1\), \(Z_1=XZ_2\), \(Z_2=XZ_3\). The curve meets the
last exceptional line at \(Z_3=1\), while the strict infinity line
is at \(Z_3=0\) and the other boundary intersection is at infinity
in this chart. With \(Y=Z_3-1\), the remaining branch has orders
\((2,m_\infty-6)\), after removing its lower even terms. Those
removals replace Y by Y minus a polynomial in X, fixing the exceptional
axis X=0 and its marked boundary points. Hence vanishing an intermediate
even coefficient does not change the marked resolution tree. The final
odd coefficient (5), (6), or (8) is a unit on its own stratum, so the
remaining fixed sequence of centres gives a simultaneous marked
resolution. This uses constant infinity type within one stratum, not
across their boundaries.

The lifted proper normalization map is an immersion: its only initial
finite critical parameters were the two resolved cusps, and the
prescribed infinity branch has also been resolved. Equivalently one
may remove the closed parameter image of its nonimmersion locus. Its
diagonal in the fibre product is therefore open (unramified) and
closed (separated). The off-diagonal double incidence is consequently
proper over the base. Pairwise-distinct triple incidence is proper for
the same reason. Nontransverse double pairs and triple incidences have
closed algebraic parameter images. Removing those images leaves exactly
\(G_m\). There is no nonproper pair incidence silently assumed closed.
Birationality has already been proved by the derivative argument.

Every nonempty \(G_m\) is path connected. It is a nonempty Zariski
open subset of an affine three-space, plane, or line. Two points in
such an open subset can be joined inside their complex affine line
while avoiding its finitely many bad points. All good repeated-U'
parameters and vanishing intermediate-even coefficients belong to the
same open locus as the other good points.

The ordinary nodes form a finite étale multisection. Blow up that
entire multisection without choosing global labels. Together with the
prescribed cusp/infinity resolutions this gives a smooth proper family
with reduced total curve-plus-infinity divisor relative normal crossing.
Along a compact path in \(G_m\), a tangent field in divisor charts can
be patched tangent to each stratum; properness gives its flow. It
identifies the actual affine complements, because every blowup centre
lies in the removed divisor. The complex normal orientation preserves
positive curve meridians up to their access conjugacies. Thus the actual
marked affine group is constant on each nonempty \(G_m\).

## 5. Two proved suppliers and the exclusion

The literal audited cyclic-complement point is
\[
 (s,a,c)=(0,0,1/6)\in G_9.
\tag{13}
\]
Hence the connected-family result gives the same cyclic group for
every member of \(G_9\).

For \(G_7\), use the following general feature of the finite rational
certificate. Its four fixed rational base loops and fixed rational
centre strands satisfy finitely many strict Rouché and endpoint
inequalities at (13). The monic quartic resultant's coefficients are
polynomial functions of \((s,a,c)\). Each inequality is continuous
in those coefficients. Therefore all of them remain strict in some
Euclidean neighborhood \(\mathcal U\) of (13), using exactly the
same rational centres, radii, paths, and signed crossing words. The
actual roots still lie in the certified disks, so the same common-base
marked relations hold. Every good parameter in \(\mathcal U\)
therefore has cyclic complement by the same arbitrary-group argument.
This does not assert an infinity isotopy at (13).

The open set \(G_7\) is nonempty. A small exact representative is
\[
 (s,a,c)=(0,1,1/6),\qquad
 U=t^4-2t^2,\quad
 V=t^6+\frac65t^5-\frac32t^4-\frac53t^3-t.
\tag{14}
\]
For \(p=x+t,q=xt\), its first divided difference is
\(p(p^2-2q-2)\). On \(p=0\), the second is
\((18q^2+25q-15)/15\), giving two distinct off-diagonal
parameter pairs. On \(q=(p^2-2)/2\), it is
\[
 -(p^2-4)(15p^3+18p^2-22)/60,
\]
giving three further off-diagonal pairs in addition to the two cusp
diagonals. The relevant discriminants and resultants are nonzero.
The exact tangent ideal has Gröbner basis
\(\{x+t^3-2t,(t^2-1)^2\}\), and the constant-leading cubic
remainder test gives triple ideal one. Both cusp fibres have single
parameter preimages. Thus (14) has exactly two ordinary cusps and five
ordinary nodes, with infinity seven.

Since \(G_7\) is a nonempty Zariski open subset of \(\mathbb A^3\),
it meets \(\mathcal U\). That supplies one good infinity-seven
cyclic-complement member, without needing a numerical epsilon threshold
or another braid production. Connected marked transport then gives
\(\pi_1=\mathbb Z\) for all of \(G_7\).

Finally the all-degree cyclic-complement consumer in the audited
[literal braid theorem](planar_jc48_sep08_two_cusp_braid.md), §5,
excludes every curve in \(G_7\cup G_9\) as a whole irreducible
Keller nonproperness curve. A hypothetical generic-degree-d cover has
at least one retained sheet, but a transitive cyclic monodromy generated
by one positive meridian has none when \(d>1\). No one-cusp ledger is
imported into this two-cusp class.

## 6. The genuine infinity-eleven residual

The remaining good locus \(G_{11}\) is real, rather than a parameter
artifact. At \(s=2\), (7) gives
\[
\begin{aligned}
 U&=t^4-\frac83t^3-2t^2+8t,\\
 V&=t^6-4t^5-\frac32t^4+\frac{548}{27}t^3-\frac{368}{9}t.
\end{aligned}\tag{15}
\]
Here the complete first pair equation gives
\[
 q=\frac{3p^3-8p^2-6p+24}{2(3p-4)};
\]
its denominator cannot vanish on the equation. The second divided
difference is
\[
 -(p^2-4)(9p^3-60p^2+112p-144)/36.
\]
The cubic has three distinct roots, avoiding the denominator and all
diagonal conditions. The exact tangent ideal is again
\(\{x+t^3-2t,(t^2-1)^2\}\); the triple ideal is one; both cusp
fibres have a single parameter. This proves exactly three nodes,
two ordinary cusps, and infinity eleven. Its full vertical discriminant
is retained in the source/output for a later actual braid supplier.

The connectedness of \(G_{11}\) does not identify its group. Nor can
one infer its group from perturbing into \(G_9\): the strict
certificate neighborhood of (13) can be chosen disjoint from the
infinity-eleven graph, since (6) is nonzero at (13). This is the precise stopping point of the proved family
exclusion.

A separate cheap sidecar at (15) used two numerical meshes, 256 and 512
subdivisions per edge, retaining all six critical projection values.
Five scout words suggest another cyclic collapse. **They are HEURISTIC
ONLY and are not dependencies of any theorem above.** The optional
`--scout` mode reproduces all six numerical words; a new rational
certificate and independent audit are required to close \(G_{11}\).

## 7. Exact controls and reproduction

The [source](../../04-computation/planar_jc48_sep08_two_cusp_family.py)
checks the complete coefficient normal form and shear; ordinary jets
without dropping repeated U' roots; all infinity coefficients and
lower-order vanishings; the excluded even double-cover control; the
vanishing even-coefficient control; and the two complete finite
representatives (14), (15). It does not infer good-locus openness or
topological transport from finite parameter tests.

```bash
python3 04-computation/planar_jc48_sep08_two_cusp_family.py
python3 -O 04-computation/planar_jc48_sep08_two_cusp_family.py
# Optional, explicitly HEURISTIC; not part of the exact output:
python3 04-computation/planar_jc48_sep08_two_cusp_family.py --scout
```

The source and [output](planar_jc48_sep08_two_cusp_family.out) are frozen.
Normal and optimized replays are byte-identical: **77 always-active gates**,
915 output bytes. The [independent analytic/source audit](planar_jc48_sep08_two_cusp_family_audit.md)
passes, including fresh normal and optimized replays.

* Source SHA256: `50daef4d2689fe62ca23e66ba4fb9f8306c46afd50983ad83477321ab58eaba7`.
* Output SHA256: `1ec6e9e83b3aac30adaf90e91ced31448c54c7a3dffa80d4c142acf313ba0332`.

The old cyclic braid and all prior frozen artifacts are unchanged.
