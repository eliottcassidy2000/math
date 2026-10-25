# The same A2 hexagon, two index-seven quotients, and the metric that makes seven colors work

**Status:** PROVED elementary identities and finite classifications; classical coloring/perfect-code facts explicitly attributed. New-to-route coordinate identification with Session60, not a literature-wide novelty claim. **No new lower bound for the chromatic number of the plane.**

## Inheritance, scope, and the concept board

The source arrangement is Session60 `hourglass-symmetry/RESULT.md`:

\[
L_0=E,\quad L_1=E+x,\quad L_2=E-x,\quad
L_3=E+z,\quad L_4=E-z,\quad
L_5=E+x+z,\quad L_6=E-x-z.
\]

The nearest proved mechanism is its exact A2 dual-root identification and order-twelve arrangement symmetry. The inherited HN control is THM-418, `lattice-chromatic-dichotomy-of-the-plane.md`, restricted here to its valid triangular-lattice assertion. Its broader all-lattice wording is not used. The canonical hostile is the graph on a single triangular lattice, which needs only three colors despite having six neighbors per point. The corrected near miss is confusing a seven-set marking with Fano incidence. The least-used sidecar is the **Euclidean metric and cell separation**, not the seven labels.

The four live concepts are (i) the root-label representation, (ii) the conjugate index-seven ideals, (iii) a perfect one-code and graph squaring, and (iv) the Voronoi metric lift. The [torus/Fano/Campedelli continuation](../torus-fano/TORUS_FANO_CAMPEDELLI.md) is developed in a companion note; this note supplies its integral coordinates and does not claim that the torus incidence is Euclidean collinearity.

## 1. One integral root lattice, with dual and primal roles separated

Put

\[
\omega=e^{2\pi i/3},\quad \Lambda=\mathbb Z[\omega],\quad
Q(a,b)=|a+b\omega|^2=a^2-ab+b^2,
\quad G=\begin{pmatrix}1&-1/2\\-1/2&1\end{pmatrix}.
\]

The six units are

\[
U=\{\pm(1,0),\pm(0,1),\pm(1,1)\},\qquad B=\{0\}\cup U.
\]

The map sending the dual affine coefficient vector of `E+ax+bz` to `a+bω` identifies the seven branch **labels** with `B`. It preserves the six-root metric and the actual dihedral permutation action. It does not turn branch lines into Euclidean points, and it does not preserve their incidence as unit-distance edges.

There is also an exact primal metric identification. In the chart `E=1`, put

\[
x=2a-b,\qquad z=-a+2b,\qquad y=a+b\omega.
\]

The inverse is `a=(2x+z)/3`, `b=(x+2z)/3`, and

\[
|y|^2=(x^2+xz+z^2)/3.
\]

Then

\[
H=\{y: |x|\le1,\ |z|\le1,\ |x+z|\le1\}
\]

is precisely the closed Voronoi cell of `Λ` at zero: these are the inequalities
`|<y,u>|≤1/2` for the three opposite unit pairs. Thus the six outer branch lines are the six supporting lines of the actual Voronoi hexagon for this chosen Euclidean metric. The central branch `E=0` is the **line at infinity in this chart**, not the center of the hexagon. It is only the central branch's *dual label* that maps to zero in `B`.

The six listed inequalities suffice for the full Voronoi cell. The six vertices are

\[
(2/3,1/3),(1/3,2/3),(-1/3,1/3),
(-2/3,-1/3),(-1/3,-2/3),(1/3,-1/3).
\]

They have norm squared `1/3`. Any nonunit nonzero lattice vector has norm at least `√3`, since `Q≡0,1 (mod3)` and `Q` is a positive integer. For such `v` and `y∈H`,
`<y,v>≤|v|/√3≤|v|²/2`; therefore its Voronoi inequality is redundant. The circumradius is `1/√3`, and the diameter is `2/√3`.

## 2. Exactly two index-seven sublattices make the flower a residue system

Define

\[
\alpha=2-\omega,\qquad\bar\alpha=3+\omega,\qquad
\alpha\bar\alpha=7.
\]

Their multiplication matrices are

\[
M_2=\begin{pmatrix}2&1\\-1&3\end{pmatrix},\qquad
M_4=\begin{pmatrix}3&-1\\1&2\end{pmatrix}.
\]

Both determinants are seven, `M₂M₄=7I`, and `MᵀGM=7G`. The corresponding ideals are

\[
I_2=(2-\omega)=\ker c_2,\quad c_2(a,b)=a+2b\pmod7,
\]
\[
I_4=(3+\omega)=\ker c_4,\quad c_4(a,b)=a+4b\pmod7.
\]

Indeed evaluation at `ω=2` or `4` is a ring homomorphism because these are the two roots of `t²+t+1` modulo seven. Each displayed ideal lies in the corresponding kernel and both have index seven, so equality follows.

**Classification.** Among all index-seven additive sublattices of `Λ`, exactly `I₂,I₄` have `B` as a complete set of coset representatives.

Every index-seven subgroup contains `7Λ`, so it is the kernel of a nonzero linear functional on `F₇²`, unique up to a nonzero scalar. The eight possibilities are `a+kb`, `k=0,…,6`, and `b`. The latter and `k=0` kill a unit. For the others one must have six distinct nonzero residues

\[
\{\pm1,\pm k,\pm(1+k)\}=\mathbb F_7^*.
\]

`k=1` repeats a value; `k=3` has `1+k=−k`; `k=5` has `1+k=−1`; and `k=6` has `1+k=0`. Exactly `k=2,4` survive. Thus twelve of the 48 nonzero linear functionals work, six color relabelings of each kernel.

In branch order `L₀,…,L₆`, the residue lists are

```
c2: [0, 1, 6, 2, 5, 3, 4]
c4: [0, 1, 6, 4, 3, 5, 2].
```

The complete split residue ring makes the mechanism transparent:

\[
\Lambda/(7)\simeq\mathbb F_7\times\mathbb F_7,
\quad a+b\omega\longmapsto(a+2b,a+4b),
\]
\[
Q(a,b)\equiv(a+2b)(a+4b)\pmod7.
\]

The six units exhaust the affine norm-one conic; their images are exactly `(u,u⁻¹)` for `u∈F₇*`. The projective conic `E²=Q(a,b)` adds the two isotropic points at `E=0`, the two split directions. This is a split-norm explanation of `1+6`, not a lower-bound proof.

## 3. The perfect-code seven and the three-color control are different predicates

Let `T` be the graph on `Λ` whose adjacency differences are `U`.
Since `B` is a residue system, every `v∈Λ` is uniquely `i+r` with `i∈I₂` and `r∈B`. Thus the closed graph-radius-one balls about `I₂` partition `Λ`. The same holds for `I₄`, and each coset is a perfect one-code. This is a **classical construction**, explicitly present in Evans–Goryainov–Konstantinova–Mednykh §4.4; Dejter gives stronger perfect-code classification results.

The seven-point flower is a clique in the **square graph** `T²`, where distinct points are adjacent when their graph distance in `T` is at most two. A common `c₂` color cannot occur at graph distance at most two: a distance-two difference is a difference of two elements of `B`, which have distinct residues. Consequently

\[
\chi(T^2)=7.
\]

However `T²` includes three Euclidean edge lengths `1,√3,2`. It is **not** a unit-distance graph. For `T` itself, `a+b (mod3)` is a proper three-coloring, and the triangle `{0,1,1+ω}` proves

\[
\chi(T)=3.
\]

The induced graph `T[B]` is the wheel with six rim vertices (12 edges); its chromatic number is three, not seven. Graph squaring is the precise operation that creates the seven-clique.

The quotient graph `T/I₂` is `K₇`: every one of the six nonzero color differences is represented by exactly one unit. The map `T→K₇` is locally bijective, a covering of 6-regular graphs with infinitely many sheets. It is **not a six-sheeted covering**. A quotient can increase chromatic number: `χ(T)=3` and `χ(K₇)=7` are compatible because only a coloring of the quotient automatically pulls back, not conversely.

## 4. Rotation descends; reflection exchanges conjugate ideals

Use the matrices

\[
R=\begin{pmatrix}1&-1\\1&0\end{pmatrix},\qquad
C=\begin{pmatrix}1&-1\\0&-1\end{pmatrix}.
\]

They represent multiplication by `1+ω` (60 degrees) and complex conjugation. They satisfy `R⁶=C²=1`, `CRC=R⁻¹`; their twelve products are the entire root-hexagon dihedral group.

On residues,

\[
c_2R=3c_2,\quad c_4R=5c_4,\quad
c_2C=c_4,\quad c_4C=c_2.
\]

Thus rotations retain either ideal and give the multiplicative `C₆=F₇*` action; every reflection exchanges `I₂,I₄`. No reflection induces a well-defined map from the **fixed** lattice quotient `Λ/I₂` to itself. A hostile witness is `α∈I₂` but `Cα=ᾱ∉I₂`.

This does not deny the original `D₁₂` action on the seven branch labels: choosing the section `B` permits *any* label permutation to be transported to a permutation of seven residues. That section-dependent action need not be induced by a map on the full quotient lattice and need not preserve quotient faces. The destroyed coordinate is the ideal/chirality choice.

For the `C₆` subgroup, the bijection from branch labels to residues is an actual intertwiner of seven-dimensional permutation representations (`1` plus the regular representation on six roots). Session60's canonical representation restricts to this same action because rotations have determinant `+1`. For reflections Session60 retains its determinant twist; a bare seven-color permutation model cannot discard that twist or the exchanged quotient.

## 5. The sharp metric sidecar produces the classical plane coloring

At center spacing one, let `H` be the cell above. The least distance between distinct same-color cells `H+i`, `H+j`, for `i−j∈I₂\{0}`, is

\[
\delta=\sqrt{7/3}.
\]

Here is an elementary exact proof. The nearest ideal vectors are `αu` for the six units `u`, with norm squared seven. By rotation it suffices to use `α=(2,−1)`. Since `H−H=2H`, the distance between `H` and `H+α` is the distance from `α` to `2H`. Its nearest point is

\[
p=(2/3,-2/3),\qquad \alpha-p=(4/3,-1/3).
\]

The residual is a positive combination `(4/3)e₁+(1/3)(−e₂)` of the two outward active normals at `p`, proving the projection condition. Its norm squared is `7/3`.

All other nonzero ideal vectors have norm squared at least 21: `Q(αβ)=7Q(β)` and nonunit `β` has `Q(β)≥3`. Their cell distances are at least `√21−2/√3>√(7/3)`. This proves the exact minimum, not merely a finite-box guess. The conjugate ideal has the same distances by reflection.

Scale the lattice and cells by `h>0`. The cell diameter and distinct same-color cell separation are

\[
D(h)=\frac{2h}{\sqrt3},\qquad \delta(h)=h\sqrt{7/3}.
\]

The two **strict** safety inequalities `D(h)<1<δ(h)` hold exactly when

\[
\sqrt{3/7}<h<\sqrt3/2.
\]

At the rational choice `h=3/4`, their squared values are

\[
D(h)^2=3/4<1,\qquad \delta(h)^2=21/16>1.
\]

Color each point by the `c₂` or `c₄` color of a cell containing it, choosing any incident cell on a boundary. Two points given the same color either belong to one closed cell (distance strictly below one), or to distinct same-color closed cells (distance strictly above one). This proves a proper seven-coloring of the entire plane, including all boundaries.

The interval is exact for the **strict diameter/separation test**. No endpoint impossibility theorem for all boundary conventions is asserted. Classical unit-diameter normalization gives `δ/D=√7/2`, exactly the value recorded by de Grey–Parts §1. The global HN lower-bound question is separate.

## 6. Fano: the missing operation is quotienting the triangular 2-cells

The historical output `05-knowledge/results/hexagonal_forbidden_s116m.out` says the centered hexagonal neighborhood *is* the Fano plane. Taken as an identification of existing incidence or symmetry this is false. The root hexagon has three collinear opposite pairs through its center, not seven projective lines. Central inversion has cycle type `1·2³`, which preserves no Fano structure, as Session60 already proves.

There is a repaired positive route: keep the **triangular 2-cells**, then quotient by `I₂`. Positive face translates map to `{t,t+1,t+3}`, and negative face translates to `{t,t+2,t+3}`, for `t∈F₇`. Each family is a Fano design because every unordered vertex pair occurs in exactly one member. The families are distinct, and 60-degree rotation/inversion exchange them. The torus quotient and its link to the Session60 special four-kernel orbit are treated in the [torus/Campedelli note](../torus-fano/TORUS_FANO_CAMPEDELLI.md). This note verifies the displayed local formulas only.

Thus the repair is not to declare the seven Euclidean points a Fano plane. It is to retain and quotient the triangular face incidence, choose an orientation class, and record the ideal/chirality sidecar. This creates a real incidence bridge while explaining exactly why no fixed Fano structure admits central inversion.

## 7. Source pins and non-imports

Sources inspected 2026-10-25:

* Aubrey de Grey and Jaan Parts, [*Tiling the plane with hexagons: improved separations for k-colourings*, arXiv:2206.12635v1](https://arxiv.org/html/2206.12635v1), §1: the classical regular-hexagon seven-color separation `√7/2` at unit diameter. The metric calculation above is an independent rederivation in the hourglass coordinates.
* Rhys J. Evans, Sergey Goryainov, Elena V. Konstantinova and Alexander D. Mednykh, [*A general construction of strictly Neumaier graphs and a related switching*, arXiv:2109.13884v6](https://arxiv.org/pdf/2109.13884v6), §4.4, printed p.9: the ideal `(2−ω)` and its seven cosets form perfect one-codes. **Do not import** the adjacent assertion of only six norm-seven elements in the entire ring: there are twelve, six in each conjugate ideal. Our exact enumeration verifies the distinction.
* Italo J. Dejter, [*Quasiperfect domination in triangular lattices*, arXiv:0903.3685v1](https://arxiv.org/pdf/0903.3685), introduction and §2: prior perfect-code classification in the triangular lattice. This note does **not** reproduce or rely on its stronger nonlinear classification.
* D. Coulson, [*On the chromatic number of plane tilings*, J. Aust. Math. Soc. 77 (2004), 191–196](https://doi.org/10.1017/S1446788700013574), pp.191–192: distinguishes lattice-sublattice Voronoi coloring from arbitrary plane coloring; it cites a seven-color lower bound for the former class. We do not import a universal seven-color lower bound.

The current lower bound for unrestricted HN is [CITED externally](../SOURCE_STATUS.md); no local Lean build or transitive axiom audit is claimed. The exact alternatives six versus seven remain open there. No graph-cover, code, Fano, or root-surface assertion here resolves that remaining choice, and no integer-square hourglass construction follows.

## 8. Reproduction and connection contract

```
python3 verify_hexagonal.py > verification.json
python3 -O verify_hexagonal.py > verification-optimized.json
cmp verification.json verification-optimized.json
```

Standard-library-only exact rational arithmetic. The complete finite universes are the 48 nonzero `F₇²→F₇` linear maps, eight subgroup kernels, twelve dihedral matrices, six Voronoi vertices/short ideal vectors, and fourteen quotient faces. The radius-two difference set and graph checks are exhaustive. A bounded `[-12,12]²` decomposition check is a redundant diagnostic, **not** the proof of the unbounded perfect-code theorem. Norm-seven enumeration is exhaustive because `Q(a,b)≥3a²/4,3b²/4` bounds both coordinates by three.

The script uses explicit exceptions, not optimization-removable `assert`. Hostiles include a bad linear map killing a graph-distance-two vector, conjugation moving a kernel, a too-small scale failing strict separation, and an unscaled cell exceeding the diameter threshold.

| Source → target | Map / preserved predicate | Information lost / required sidecar |
|---|---|---|
| Branch dual labels → root flower | `(a,b)→a+bω`; exact `D₁₂` root metric | Primal incidence, individual square roots, cover characters |
| Primal outer lines → Voronoi cell | `x=2a−b,z=−a+2b`; exact supporting inequalities | Must choose `E=1` and Euclidean metric; central line is infinity |
| Lattice → seven residues | `c₂` or `c₄`; complete closed-neighborhood colors | Lengths, winding, chirality; quotient may increase chromatic number |
| Perfect code → square-graph coloring | Distinct radius-two differences | Adds lengths `√3,2`; not a unit-distance lower-bound carrier |
| Seven cell colors → plane coloring | Assign entire Voronoi cells | Requires strict diameter/separation and boundary ownership |
| Triangular faces → Fano designs | Quotient faces, then choose orientation | Not Euclidean collinearity; 60-degree rotation swaps designs |
