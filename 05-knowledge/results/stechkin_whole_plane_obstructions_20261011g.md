# Stechkin's mixed distances, whole-plane coverage, and exact repair obstructions

2026-10-11. **CITED:** Stechkin's mixed-distance construction, published in
Raiskii's 1970 paper. **PROVED:** the interior odd-cycle obstruction, the
rescaling coverage calculation and periodic obstruction, the full unit-path
component statement, and the restricted rational-direction regularity contrast.
**FINITE-EXACT:** the accompanying rational arithmetic controls. **OPEN:** a
six-colouring of the whole unit-distance plane, or a seven-colour lower bound.
No novelty claim is made for the elementary general arguments.

## 1. Inheritance and the missing synchronized distance

The user's reference is **Stechkin**, not Steinhaus. The original paper ends
with Stechkin's six-set tiling: when every tile side is 1, four hexagon classes
avoid distance 2 and two open-triangle classes avoid distance 1. Rescaling
the complete construction by one half gives type

\[
 (1,1,1,1,1/2,1/2).                                      \tag{1}
\]

The assertion is that each class avoids its assigned distance. It does not
say all six classes avoid one common distance. Raiskii's distinct theorem
in the same paper says **some** class of any three-set plane cover realizes
every positive distance. [Original paper, pp.319–323](https://www.mathnet.ru/php/getFT.phtml?jrnid=mzm&paperid=9512&what=fullt),
[publication record](https://doi.org/10.1007/BF01093113).

For the displayed colouring convention, the two triangle classes are the
two triangle orientations. Soifer's authored exposition states this
explicitly and uses the same boundary ownership: triangle interiors are
open, while hexagons own the remaining boundaries with specified vertex
exceptions. [Soifer, *Mathematics Competitions* 16(1), 2003, printed pp.26–27](https://www.wfnmc.org/Journal%202003%201.pdf).
The controls below use only strict interiors, so no vertex ownership is
inferred from a picture.

Closest proved mechanism: the exact residual-bipartite and boundary-parity
criterion in [the six-colour reframe](plane_six_reframes_parity_20261011b.md).
Canonical hostile: an internal monochromatic odd cycle survives even when
the boundary is empty. Corrected near miss: the
[Vitali recovery](vitali_atom_recovery_20261011.md) distinguishes an ordinary
lattice fundamental domain from a transversal of a dense countable subgroup.
Least-used sidecar: **the covered set after applying different dilations to
different colours**. The board is mixed-distance type / covered set / odd
cycle / translation component / regularity. Anchor is an actual repair of
(1), niche is coset extension, wildcard is a periodic dilation obstruction.

## 2. Exact unit triangles inside the two triangle colours

Use the unit triangular-lattice basis

\[
 u=(1,0),\quad v=(1/2,\sqrt3/2),\quad
 \|au+bv\|^2=a^2+ab+b^2.
\]

Hexagon centres form \(\Lambda=\mathbb Zu+\mathbb Zv\). One gap triangle has
vertices \(u/2,v/2,(u+v)/2\), and its centroid is
\(p=(u+v)/3\). Its translates by \(u,v\) are triangles of the same
orientation. Therefore

\[
                 p,\ p+u,\ p+v                            \tag{2}
\]

are strict interior points in **one original triangle colour**, with all
three pairwise distances equal to 1. The opposite orientation has centroid
\(2(u+v)/3\) and supplies a second such witness. Translating (2) by any
lattice vector preserves these properties. These are actual triangles in
the Euclidean plane, not quotient cycles.

Let \(T\) be the union of the two triangle colours. If the four original
hexagon colours are retained, assigning just two remaining colours to
\(T\) is possible exactly when its induced unit-distance graph is
bipartite. Witness (2) proves it is not. In fact **every perfect pairing of
the six whole original classes into three pairs fails**: the pair containing
either original triangle class still contains its unit triangle. The script
checks all 15 pairings and all eight binary assignments on that triangle.

This rules out only these repairs. A successful new six-colouring could
split old classes pointwise and change which points belong to the three
pairs. No assertion that every repartition fails is made.

## 3. Normalizing each colour separately loses coverage

Suppose \((C_i)\) is a partition, and \(C_i\) avoids distance \(d_i>0\).
Then \(A_i=d_i^{-1}C_i\) avoids unit distance. If the sets \(A_i\) cover the
plane, selecting the least available index gives a unit-distance colouring:
each resulting colour is a subset of a unit-safe \(A_i\). Overlap is
harmless. **Coverage is the additional premise.**

In (1), normalization about the hexagon centre 0 gives

\[
 \bigcup_i A_i=(\mathbb R^2\setminus T)\cup2T,
 \qquad \text{uncovered set}=T\setminus2T.                  \tag{3}
\]

Already \(p\) from (2) is missing: \(p\in T\), while
\(p/2=(u+v)/6\) lies in the interior of the hexagon at 0, hence not in
\(T\). In particular \(\|p/2\|^2=1/12\), less than the squared hexagon
inradius \(3/16\).

There is an exact density statement. In lattice coordinates, a fundamental
\([0,1]^2\) cell has two open triangles

\[
 \begin{split}
 T_-&=\{0<a,b<1/2:a+b>1/2\},\\
 T_+&=\{1/2<a,b<1:a+b<3/2\}.
 \end{split}                                               \tag{4}
\]

Each has coordinate area \(1/8\). In the common period \([0,2]^2\),
\(T\) consists of eight such triangles with total area 1. The intersection
\(T\cap2T\) consists, up to boundaries, of exactly \(T_+\) in the
\((0,0)\) cell and \(T_-+(1,1)\), with total area \(1/4\). Thus

\[
 \operatorname{dens}(T\setminus2T)=3/16,
 \qquad \operatorname{dens}(\bigcup_i A_i)=13/16.           \tag{5}
\]

The coordinate Jacobian \(\sqrt3/2\) cancels from these ratios. Exact
rational convex clipping independently computes all 16 pairwise triangle
intersections: two have area \(1/8\), the other fourteen area 0.
Equation (5) concerns this specified hexagon-centred dilation.

**All-centres periodic obstruction.** More generally let \(T\) be a
measurable lattice-periodic set with density strictly between 0 and 1.
For any vector \(w\), it is impossible that \(T\subseteq2T+w\), even
modulo a null set. Both sets are periodic for \(2\Lambda\) and have the
same density, so inclusion would give equality almost everywhere. Writing
\(f(x)=f((x-w)/2)\) for its indicator, every period \(\lambda\) of
\(f\) yields the period \(\lambda/2\). Iteration gives the dense period
group \(\bigcup_m2^{-m}\Lambda\). Convolving the bounded measurable
indicator with a smooth compactly supported approximate identity gives
continuous functions invariant under this dense group, hence constants.
Approximation in local \(L^1\) forces \(f\) constant almost everywhere,
contradicting its density. Therefore changing the common centre of the
normalization cannot remove every hole. This is a restriction on this
specific dilation repair, not on all possible recolourings.

## 4. Why the whole-plane unit graph has no residual coset decomposition

For a prescribed symmetric set of allowed directions \(D\), connected
components of its translation graph are exactly the cosets of
\(H=\sum_{d\in D}\mathbb Zd\). Translating a colouring of \(H\) to all
cosets is valid only because every allowed edge stays within a coset. A
chosen transversal need not be measurable when \(H\) is dense/countable.

For **all real unit directions**, however, \(H=\mathbb R^2\). If
\(0<\|x\|\le2\), let \(Jx\) be the right-angle rotation of \(x\). Then

\[
 x=\left(\frac{x}{2}+\sqrt{\frac1{\|x\|^2}-\frac14}\,Jx\right)
  +\left(\frac{x}{2}-\sqrt{\frac1{\|x\|^2}-\frac14}\,Jx\right),\tag{6}
\]

and each summand has norm 1. Subdivide a longer vector into finitely many
vectors of norm at most 2 and concatenate these two-edge paths. Thus the
whole unit-distance graph is connected. A countable-direction or field
coset extension does not address the omitted real unit directions.

The incoming [cyclotomic unit-component note](cyclotomic_unit_component_20261011.md)
has a different, explicitly typed source: directions in the countable
field \(\mathbb Q(\zeta_{55})\). Its component localization and translated
coset colourings concern that field graph. Equation (6) uses arbitrary real
unit directions; its square roots need not lie in the chosen field. We do
not re-audit the cyclotomic localization theorem here.

## 5. A precise Vitali-style regularity contrast, not a six-colour premise

Take only rational-coordinate unit directions on the whole plane. Every
such direction \((a,b)\) has odd reduced denominators and
\(a+b=1\pmod2\) in \(\mathbb Z_{(2)}\). For example the parametrization
\(((1-t^2)/(1+t^2),2t/(1+t^2))\), \(t\in\mathbb Q\), proves this after
cancelling a possible factor 2; the omitted point \((-1,0)\) also has the
property. Hence the additive parity map is well defined on their generated
subgroup \(H\) and changes on every allowed edge. Choosing one
representative of each coset gives an abstract proper two-colouring of
this restricted-direction graph.

Such a two-colouring cannot be Lebesgue measurable. Any proper binary
colouring is invariant under every even-length direction word. That
subgroup is dense: it contains every difference of two rational unit
directions, whose closure contains all differences of real unit vectors,
and those differences fill the radius-2 disk. If its binary indicator were
measurable, the convolution argument from §3 would make it constant almost
everywhere. A unit translation must instead exchange the two colours,
which is impossible for an almost-everywhere constant binary function.

This exact contrast requires no claim about arbitrary six-colourings being
measurable. It also makes no use of the nonmeasurability of a transversal
as a proof of a chromatic bound. Enlarging from rational unit directions
to all real unit directions introduces the equilateral-triangle edges
that already preclude two colours. Finite and discrete selectors in the
[corrected Vitali canon](vitali_selector_certificate_20261005.md) retain
their separate types.

## 6. Reproduction and next obligation

Run from the repository root:

```text
python -B -X utf8 04-computation/experiments/stechkin_whole_plane_obstructions_20261011g.py
python -B -O -X utf8 04-computation/experiments/stechkin_whole_plane_obstructions_20261011g.py
```

The finite universe is 242 translated centroid triangles, all 15 whole-class
pairings, all eight binary triangle assignments, 16 polygon intersections,
four finite rational grids, 300 rational-circle parameters plus the omitted
stereographic direction, bounded two-unit norm identities, and 13 malformed
or native-guard hostiles. These checks support the displayed exact proofs;
they do not enumerate the plane or decide six-colourability.

Normal and optimized runs agree with the saved output at **3,617 checks**.
LF-normalized SHA256:
`c1e6401398d1283cc078f27876b5cf2a4bb25f2b1ce88c7aa37d5a0565a315c3`.

The useful next target is now precise: construct a **new pointwise**
three-pair assignment for which every pair's entire unit graph is bipartite,
retaining all crossing-edge parity constraints. Permuting the six original
whole classes cannot do it, and separate distance normalization must pay
the coverage condition (3). The companion global seam work owns that
remaining compatibility problem.
