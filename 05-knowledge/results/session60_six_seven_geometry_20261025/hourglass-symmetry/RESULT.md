# The hourglass root hexagon: exact symmetry and a singular Kummer quotient

**Status:** PROVED over `C`, with the classical Shioda-sextic identification
CITED from a primary source; FINITE-EXACT coordinate checks. No novelty claim.
The statements concern the specified line arrangement and covers, not every
surface with the same numerical invariants.

## 0. Inheritance and connection board

The inherited result is Session59,
`toric-shared-roots/TORIFICATION_OBSTRUCTION.md`, §2: the seven-root
hourglass has 24 ordinary double points, and its minimal resolution has
`K²=16`, `p_g=7`, `q=0`; every six-root coordinate deletion resolves to a K3.
Its triple-point calculation is recovered below independently. The inherited
six-root row was one deletion, not a classification of all seven deletions.

The live objects are (i) the root hexagon, (ii) canonical deck characters,
(iii) six-line double planes, (iv) genus-two Jacobians and elliptic factors,
and (v) the repository's Bolza/CM sidecar. The corrected near miss is that
an isogeny or finite map to a Kummer surface does **not**, alone, imply that
the source is Kummer. The hostile arithmetic example is the nontrivial
scalar `27` in the explicit Shioda equation below. The least-used sidecar is
the exact location of the branch points on the dual conic, not merely their
number.

Write

```
L0=E,
L1=E+x, L2=E-x, L3=E+z, L4=E-z,
L5=E+x+z, L6=E-x-z.
```

Let `H` be the projective seven-root surface with coordinates `[r0:...:r6]`
whose squares are the seven displayed forms on a plane. Let `Y` be its
minimal resolution. The deck group is
`G7={±1}^7/{±(1,...,1)} ≅ (C2)^6`. Write `S_j` for the minimal resolution
of the six-root projection omitting `r_j`. We reserve `D0` for the minimal
resolution of the **double-plane**

```
w²=B(E,x,z),
B=(E²-x²)(E²-z²)(E²-(x+z)²).                       (1)
```

`S_0` and `D0` are different presentations and their map has degree 16;
identifying them requires more than the calculation of their invariants.

## 1. The projective arrangement automorphism group is exactly D12

Here `D12` means the dihedral group of **order twelve**.

**Theorem 1.** The projective transformations of `P²_C` preserving the seven
unlabelled lines form `D12 ≅ W(A2)×C2`. All are defined over `Q`, and after
normalizing `E` they act by integral determinant-`±1` matrices on `(x,z)`.
They fix `L0` and are transitive on `L1,...,L6`.

**Proof.** There are exactly three triple points, with incident label sets
`{0,1,2}`, `{0,3,4}`, `{0,5,6}`. Every other intersection is double. Thus
`L0` is intrinsically the only line containing three triple points; each
other line contains one. Every arrangement automorphism fixes `E=0` and
acts affinely on `E=1`.

The three pairs of parallel affine lines are intrinsic. The midpoint lines
of these pairs are `x=0`, `z=0`, `x+z=0`; affine transformations carry such
midpoint lines to midpoint lines. Their common intersection is the origin,
so the translation part vanishes. The linear action on dual coefficients
therefore permutes

```
R={±(1,0), ±(0,1), ±(1,1)}.                          (2)
```

These six points form a regular hexagon for the positive definite quadratic
form `a²-ab+b²`. A linear symmetry is determined by the images of `(1,0)`
and `(0,1)`. Checking the thirty ordered distinct images (discarding
singular pairs) gives exactly twelve. Equivalently, a symmetry preserves
adjacency in the convex hexagon, giving an upper bound of twelve. The
matrices

```
R = [1 -1; 1 0],             S = [0 1; 1 0]
```

realize twelve distinct symmetries and satisfy `R^6=S²=1`, `SRS=R^-1`.
The corresponding base matrices are their transposes. This proves both
completeness and the rational integral realization. The Weyl subgroup has
order six; adjoining the central `-I` gives the full twelve-element group.
∎

**Corollary 1.1.** Pure permutations of the seven root coordinates give a
split subgroup

```
(C2)^6 ⋊ D12  ⊂  Aut(Y),          order 768.           (3)
```

Indeed, the normalized linear forms are permuted with no scalar factors,
so their square roots may be permuted literally. These permutations
normalize independent sign changes, intersect them trivially in projective
space, and lift to the minimal resolution. **This is not a claim that (3)
is the full automorphism group of `Y`.**

**Corollary 1.2.** The seven labelled six-root deletions have two `D12`
orbits: `{0}` and `{1,...,6}`. The model for `S_0` is smooth. Each of the
other six models has eight ordinary double points and the same K3
isomorphism type as every other member of that orbit.

Deleting `L0` removes all three triple points. Deleting a finite line
removes exactly one triple point and leaves two. The inherited square-root
cover node formula is `2^(6-4)=4` nodes above each triple point, hence
`0` or `8` nodes. The group explicitly identifies the six outer deletions.
No assertion of non-isomorphism between the resolved central and outer K3s
follows just from different singularities in these chosen models.

## 2. The canonical representation has a determinant twist

**Theorem 2.** On `H^0(Y,K_Y)`, the pure permutation lift of a base symmetry
`A` acts by

```
det(A) ⊗ (permutation representation on r0,...,r6).  (4)
```

A projective deck sign vector `ε=(ε0,...,ε6)` acts on the `j`th canonical
basis vector by

```
χ_j(ε) = product_{i≠j} ε_i.                         (5)
```

**Proof.** The complete intersection has four quadrics in `P^6`; adjunction
gives `ω_H=O_H(1)`, and its ordinary double points are crepant. Residue
adjunction realizes its seven canonical forms as `r_j` times one common
residue generator. For a coordinate permutation `P`, the residue multiplier
is the ambient determinant divided by the determinant on the four-dimensional
space of quadratic relations. The exact sequence from relations to seven
square-coordinate forms to the three-dimensional base gives that ratio as
`det(A)`. For a deck sign vector, the quadrics are fixed and the ambient
determinant is `∏ε_i`; multiplying by the numerator `r_j` gives (5).
Changing all seven signs changes neither (5) nor the projective action.
∎

The representation is therefore **not** the bare permutation representation.
For example a reflection has determinant `-1`. Formula (5) also exhibits the
seven even-weight-six characters underlying the seven K3 quotients. The
residue argument is independently recovered in the Session60 deck lane.

## 3. Six tangent lines and the exact genus-two branch set

**Theorem 3.** The six finite-line dual points lie on the smooth dual conic

```
a0²-a1²+a1*a2-a2²=0.                                (6)
```

The corresponding primal tangency conic is

```
3E²-4x²-4xz-4z²=0.                                  (7)
```

Their six parameters on (6) are

```
B6={infinity,-1,0,1/2,1,2}.                          (8)
```

Consequently `D0` is a Kummer surface of a genus-two Jacobian over `C`.

**Proof and explicit map.** Parametrize (6) by

```
[a0:a1:a2] = [t²-t+1 : t²-1 : t(t-2)].              (9)
```

In the order of (8), these are the dual points of
`L5,L3,L2,L6,L4,L1`. Nondegeneracy of (6) follows from its matrix determinant
`3/4`; inverse-matrix duality gives (7). Thus no three of these six lines
are concurrent. The standard Shioda-sextic construction identifies the
resolution of a six-tangent-line double plane with the Kummer surface of
the genus-two double cover of the conic branched at those tangency points.
For this classical step see Clingher–Malmendier, §4.1, Lemma 4.2 and Remark
4.3, printed pp.15–16; see the explicit independent equation below. ∎

A convenient rational genus-two model with branch set (8) is

```
C: y²=f(t)=t(t+1)(t-1)(t-2)(2t-1).                  (10)
```

The polynomial is squarefree of degree five, so infinity is the sixth
branch point and its smooth projective model has genus two.

### 3.1 Exact Shioda map and rational-twist boundary

For a point `[E:x:z]`, put

```
q(t)=(E+x+z)t²-(E+2z)t+(E-x),     q2=E+x+z.         (11)
```

This is precisely the pairing of `[E:x:z]` with the conic parametrization
(9). The map from `[E:x:z]` to the coefficients of `q` is an invertible
rational linear map (determinant `-3`). Its evaluations are

```
q(-1)=3(E+z), q(0)=E-x, q(1/2)=3(E-x-z)/4,
q(1)=E-z,    q(2)=3(E+x), q(infinity)=q2=E+x+z.
```

In particular,

```
q2 * product_{b=-1,0,1/2,1,2} q(b) = (27/4) B.      (12)
```

Generically let `t1,t2` be the two roots of `q`, and take points
`(t1,y1),(t2,y2)` on (10). Then

```
W = y1*y2*q2^3,             W² = 27 B.              (13)
```

The coefficient coordinates and `W` generate the invariant function field
of `(C×C)/<swap,(hyperelliptic,hyperelliptic)>`. Indeed they retain the
unordered two `t`-coordinates and the product of the two square roots;
the generic ambiguity is exactly that four-element group. This quotient
is birational to `Jac(C)/{±1}`, by the symmetric-square Abel map. This gives
a function-field proof of the Shioda identification in this case.

**Arithmetic boundary.** Equation (13), not `W²=B`, is the rational Kummer
model supplied by this particular `C/Q` and these coordinates. The original
(1) is identified with it by `W=3 sqrt(3) w`, over `Q(sqrt(3))` and hence
`C`. This argument does **not** assert that the two rational surfaces are
isomorphic over `Q`, or that they are not: the computation exposes a twist,
not a twist-classification theorem. Twisting (10) by a rational scalar
multiplies (13) by a square and does not erase this issue.

### 3.2 The branch hexagon is the sixth-root curve, not Bolza

Let `δ=sqrt(-3)`, `α=(1+δ)/2`, `β=(1-δ)/2`. The Möbius transformation

```
u=(t-α)/(t-β)
```

sends (8) to the six roots of unity, because

```
(t-α)^6-(t-β)^6 = -3δ f(t).                         (14)
```

Thus setting `v=sqrt(-3δ)y/(t-β)^3` identifies (10) over `C` with

```
C6: v²=u^6-1.                                      (15)
```

The original hexagon rotation becomes `t↦(t+1)/(2-t)`, a six-cycle on (8).
No rationality of (15)'s coordinate change is asserted.

This is **not** the repository's Bolza curve `v²=u^5-u`. All crossratios of
four points in (8) are real. The Bolza branch set
`{0,infinity,1,-1,i,-i}` has crossratio `(0,infinity;1,i)=-i` (using
`(a,b;c,d)=(a-c)(b-d)/((a-d)(b-c))`). The property that all four-point
crossratios are real is Möbius-invariant, so the two branch sextuples are
not equivalent and the genus-two curves are not isomorphic over `C`.

The closest repository comparison is
`THM-4012-weighted-leading-face-good-elliptic-factor-observer.md`, §5,
and `THM-4226-dense-weight-thirteen-primitive-cm-bolza-planar-jacobian-exclusion.md`:
there the Bolza Jacobian has CM elliptic factors with field `Q(sqrt(-2))`;
here the factors below have field `Q(sqrt(-3))`. The actual shared mechanism
is extracting elliptic quotients and retaining their CM field, not an
identity of the curves or a transfer of the planar-Jacobian conclusions.

## 4. Elliptic quotients force Picard rank twenty

**Theorem 4.** `Jac(C)_C` is isogenous to `E0×E0`, where `E0` is any
elliptic curve of `j=0`. Therefore

```
rho(Jac(C)_C)=4,       rho(D0)=20,       rank T(D0)=2. (16)
```

**Proof.** On (15) the degree-two maps

```
(u,v) -> (u²,v)                  to V²=U³-1,
(u,v) -> (u^(-2),v/u³)           to V²=1-U³          (17)
```

extend to morphisms of smooth projective curves. Their pulled-back
elliptic differentials are respectively `2u du/v` and `-2 du/v`, a basis
of `H^0(C6,Ω^1)`. Hence the induced homomorphism from its Jacobian to the
product of these two elliptic curves has invertible differential, and is
an isogeny. Both targets have `j=0` and are isomorphic over `C`.

For completeness, `E0` has its order-three complex multiplication `ζ3`.
On `E0×E0` the two fibers, the diagonal, and the graph of `ζ3` have
intersection matrix

```
[0 1 1 1]
[1 0 1 1]
[1 1 0 3]
[1 1 3 0],                    determinant -3.
```

The diagonal/graph intersection number is `deg(1-ζ3)=3`; the other entries
are immediate. These four divisor classes are independent. Since an
abelian surface has `h^(1,1)=4`, its Picard rank is exactly four; an isogeny
preserves the rational algebraic cohomology and thus Picard rank.

Blow up the sixteen two-torsion points of this abelian surface, then divide
by the lifted inversion involution. On the original `H²`, inversion acts
trivially; all sixteen exceptional classes are algebraic. The rational
cohomology of the smooth quotient is the invariant rational cohomology of
the blow-up. Consequently the Kummer resolution has Picard rank
`4+16=20` and transcendental rank `22-20=2`. The geometric identification
in Theorem 3 proves (16). ∎

A K3 of Picard rank twenty is traditionally called a **singular K3**;
this does not mean its smooth minimal resolution has singular points.

**Transfer, explicitly conditional on the separately proved deck theorem.**
The parity-even deck subgroup `(C2)^4` of `S_0` is symplectic and its quotient
resolves to `D0`. The Session60 deck lane proves that pullback identifies
`T(D0)_Q` and `T(S_0)_Q`, scaling the pairing by sixteen. Combining that
result with (16) gives `rho(S_0)=20`. The assertion that `S_0` itself is
Kummer requires its separate diagonal-net/Gale-dual proof; it is not inferred
from this transfer alone.

## 5. Exact verification, scope and stopping condition

Run in this directory:

```
python3 verify.py > normal.json
python3 -O verify.py > optimized.json
cmp normal.json optimized.json
```

Python 3 and SymPy are required. The two streams are byte-identical; 116
explicit checks survive optimization. The universe includes all 21 line
pairs, all thirty ordered distinct possible basis images in the root
hexagon, and all seven six-root deletions. Symbolic checks include conic
nondegeneracy/duality, all six branch labels, the Möbius six-cycle, the exact
scalar (12), identity (14), elliptic-map equations and differential
coefficients, and the nonzero NS intersection determinant. Hostiles include
the central dual point and a mutated sixth line, both outside the conic,
and the nonreal Bolza crossratio. Computation verifies the displayed algebra;
the geometric implications are proved in the text, not inferred from a
finite list of samples.

A second standard-library-only checker, `verify_independent.py`, exhausts all
120 possible ordered images of three branch points, constructs the unique
rational Möbius map, and independently obtains exactly twelve stabilizers,
generated by `t↦(t+1)/(2-t)` and `t↦1-t`. It also enumerates all 360 ordered
four-point crossratios (15 distinct values), and checks the coefficient-map
inverse on 125 integer controls. Its 372 explicit checks give identical
`independent-normal.json` and `independent-optimized.json` in normal and
optimized Python. This is an independent branch-configuration audit, not an
independent proof of the cited geometric Kummer theorem.

The proved map contracts/desingularizes branch information and forgets
individual root signs. To restore a root cover, one needs the labelled
square-root character extension, not just the genus-two curve. A rational
point on the Kummer quotient is not automatically a rational point with all
six individual squares, and a geometric complex isogeny gives no integral
magic-square construction. The integer-square hourglass problem is
untouched. Exact integral transcendental lattices, rational twist
classification, and full surface automorphism groups remain outside scope.

## 6. Primary-source pin

Adrian Clingher and Andreas Malmendier, *Normal forms for Kummer surfaces*,
arXiv:1809.02003v2 (5 August 2019), §4.1, Lemma 4.2 and Remark 4.3,
printed pp.15–16 (PDF indices 14–15):
https://arxiv.org/pdf/1809.02003 . Primary text inspected 2026-10-25.
The imported result is solely the symmetric-square/Shioda identification of
a six-tangent-line double plane with a Jacobian Kummer; the coordinates,
scalar-twist audit, elliptic maps and divisor-rank computation are explicit
above. Published chapter record:
https://doi.org/10.1017/9781108773355.006 .
