# Planar Jacobian counterexample search through dephasing gain graphs and square-row completion

**Research synthesis, 2026-08-18.**  Status labels are local to each claim.
`JC(2)` remains **OPEN**.  No counterexample is claimed here.  The session
produced one repaired canonical proof, a new root-multiplicity-sensitive
constraint on a hypothetical counterexample, an exact low-degree jet
closure, several refuted construction architectures, an exact Catalan
near-counterexample, and several concrete open counterexample programs.

Companion:
[`jc2_dephasing_square_rows_synthesis_codex_20260818.py`](../04-computation/jc2_dephasing_square_rows_synthesis_codex_20260818.py),
with stored output
[`jc2_dephasing_square_rows_synthesis_codex_20260818.out`](../05-knowledge/results/jc2_dephasing_square_rows_synthesis_codex_20260818.out).

## 1. Outcome first

The strongest exact outcomes are:

1. **PROVED correction.**  THM-3025 previously confused the full Keller
   constant `Jac(P,Q)` with the possibly-zero subleading form
   `j=Jac(H,Q_(m-1))`.  MISTAKE-422 records the error.  The conclusion
   `K>=2 => W=0` survives by a `j!=0` / `j=0` case split and is now proved in
   the repaired
   [THM-3025](../01-canon/theorems/THM-3025-W-is-forced-to-vanish-for-every-jc2-counterexample-candidate.md).
2. **PROVED new derived constraint; literature novelty not asserted.**  In
   the first Euclidean chamber `b<a<2b`, if
   `H=prod_L L^e`, then a hypothetical counterexample must satisfy

   ```text
   D_k(H) := prod_(L^e || H) L^ceil((e(k-1)+1)/2)  |  Q_(m-1),
   k=2b-a.
   ```

   This strictly improves the subleading codimension and remains sensitive
   to repeated roots.
3. **PROVED within a stated stratum + FINITE-EXACT replay.**  There is no
   `(deg P,deg Q)=(6,4)` Keller pair with squarefree quadratic leading form
   `H`.  This is a short homogeneous-jet certificate, not a new global
   degree theorem: the corrected Appelgate--Onishi result and Moh already
   cover this degree pair.
4. **PROVED target-pencil spectrum squeeze** (THM-3550's hostile audit is
   still in progress).  By
   [THM-3544](../01-canon/theorems/THM-3544-planar-keller-target-pencil-total-degree-six-floor.md),
   every nonzero linear combination `sP+tQ` in a nonautomorphic complex
   planar Keller pair has total degree at least six; by
   [THM-3550](../01-canon/theorems/THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor.md),
   none has prime degree and a reduced target basis has distinct degrees
   `n<m` with `m>=8`.  Applying the corrected Appelgate--Onishi theorem in
   every target direction strengthens this to at least three prime factors
   per pencil degree; Moh forces the reduced height `m` above 100.  The
   Heitmann/Guccione--Guccione--Valqui restriction further gives
   `gcd(n,m)>=16` and rules out twice-prime gcd.  Finally, the cited 2022
   sub-`125` classification leaves only the reduced pair `(72,108)` (up to
   transpose) below height `125`; this is the live first degree passport.
5. **PROVED strong-dephasing limit.**  A Hermitian hopping matrix with
   `|H_xy|^2=c_xy`, under uniform site dephasing of coherence-decay rate
   `lambda`, has population generator

   ```text
   (2/lambda)L_c + O(lambda^-2).
   ```

   Thus `p(lambda*tau/2) -> exp(tau L_c)p(0)`.  Kirchhoff theory is exactly
   the leading slow limit.  Phases first return through magnetic loop
   corrections: triangles contribute at the next order.
6. **PROVED exact bridge, with a decisive loss ledger.**  The coefficient
   equation `Jac(P,Q)=1`, for fixed `P`, is a complex bipartite gain graph
   whose rows are grouped by equal exponent-vector sums.  Strong dephasing
   keeps squared edge magnitudes but discards the cycle holonomies that
   control cancellation, rank, and range.  A resistor graph is therefore a
   useful pruning shadow, not a lossless JC carrier.  THM-3548 turns this into
   quantitative polygon, Segre-rank, dark-plaquette, shear-frustration, and
   conditional Puiseux gates.
7. **PROVED arithmetic merge.**  The user's triangular square rows and odd
   equal-sum rows are related by the Gaussian--Hadamard map

   ```text
   M(x,y)=(y-x,y+x),  det M=-2,  M^T M=2I.
   ```

   Hence odd-row norms are exactly twice the corresponding triangular-row
   norms.  The filler `3,11,27,51,...` is one more than the transformed row
   maximum.  It is a collision-free scalar sentinel, while the diagonal
   exponent is the functorial support filler.
8. **PROVED + VERIFIED-EXACT near-counterexample and cubic normal form.**  The separated
   thickening `P=v^2+A(w)`, `Q=v^3-v+vB(w)` retains the transverse collision
   `(+/-1,0)->(1,0)` and has constant Jacobian exactly when
   `B=3A/2` and `A-3A^2/4=kappa*w`.  The unique normalized solution is an
   infinite Catalan series, never a polynomial in `w`.  After adjoining its
   square root, however, the map polynomializes and is affinely equivalent
   to the universal depressed-cubic marked-root cover
   `(t,p)->(p,-t^3-pt)`.  This reveals a third simple point in the collision
   fibre and proves that any polynomial repair fixing the cubic ramification
   line pointwise still ramifies at its cusp.  Mixed repairs that move that
   line at order zero remain open.
9. **PROVED + VERIFIED-EXACT punctured near-counterexample.**  On the curved
   collision surface of the fixed three-variable Keller map, explicit
   coordinates reduce the restriction to the finite etale double cover
   `(s,b)->(b,4s^2)` with `s!=0`.  Filling the missing divisor forces
   ramification, and no etale affine-plane filling can preserve that exact
   quadratic function-field extension.
10. **PROVED/FINITE-EXACT torus-repair squeeze; OPEN cell.**  Polynomial
   corrections of the inherited quotient fail one-sidedly in every degree,
   through total correction degree three, and whenever affine in the
   transverse variable.  Importing the 2022 degree list moves the first
   row-dominant `4:5` boundary not already excluded to `(100,125)`, where
   `u^25w^2|P_99` and `u^49w^3|Q_124`; only `73` homogeneous quotient
   coefficients remain in the lower subleading row.
11. **PROVED inverse-cubic owner + FINITE-EXACT positive seed; OPEN
   projection.**  The cusp-square packet is an exact `1+2` root owner: its
   inverse cubic factors into the marked root `-1/v`, which escapes at
   `v=0`, and a quadratic Kummer pair, with discriminant `-4LS^2`.  An
   observable dual cubic marks `y` and carries the same square class `-L`.
   An explicit four-coordinate packet `A^2 -> A^4` is everywhere immersive,
   although none of its six natural two-coordinate projections, nor any
   constant-linear projection, is Keller.  A nonlinear descending,
   integrable, decomposable projection is the strongest self-contained
   counterexample search seed found here.
12. **PROVED + VERIFIED-EXACT mixed-Catalan squeeze.**  In the polynomial
   thickening of `(v^2,v^3-v)`, transverse widths one and two are impossible
   in every coefficient degree.  Width three is impossible through
   coefficient degree five; degree six is merely the first internally
   unclosed recurrence cell.  Combining this with the live height-108
   gate, a genuinely viable width-three candidate needs coefficient cap at
   least `105`, since both total degrees are at most `D+3`.
13. **PROVED all-degree displayed-graph no-go.**  For the fixed cubic Keller
   map, restriction to every displayed polynomial graph `z=h(x,y)` has an
   explicit nonzero positive-degree tangential Jacobian top row.  THM-3546
   therefore forbids descent between displayed source and target graphs in
   every degree.  Nonlinear ambient coordinates and nongraph coordinate
   hypersurfaces remain open.

## 2. Inheritance pass and concept board

### Anchor / Niche / Wildcard

- **Anchor:** squeeze the degree, leading form, subleading layers, and support
  architecture of a hypothetical planar counterexample.
- **Niche:** turn the linear response map `Q -> Jac(P,Q)` into a magnetic
  coefficient graph and determine exactly what survives dephasing.
- **Wildcard:** identify the real structure behind the two triangular arrays,
  the filler sequence, sums of two squares, and `4^a(8b+7)`.

### Inherited objects

- Closest proved JC mechanism:
  [THM-3016](../01-canon/theorems/THM-3016-jacobian-pair-cross-term-rigidity-at-subleading-order.md),
  the leading-form Pluecker identity and first Euclidean block.
- Corrected near miss: the equality version of the valuation tower in
  [HYP-9070](../05-knowledge/hypotheses/HYP-9070-jc2-leading-form-circuit-and-the-euclidean-depth-search-order.md)
  is refuted; only the stated first-block divisibility survives.  The former
  THM-3025 division is separately repaired by MISTAKE-422.
- Canonical hostile:
  [THM-2045](../01-canon/theorems/THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate.md),
  where a unique highest response row leaf cannot cancel.
- Exact response quotient:
  [THM-2230](../01-canon/theorems/THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient.md).
  Its fibres are target-shear orbits in a module, not a quotient ring.
- New proved degree/collision sidecars:
  [THM-3544](../01-canon/theorems/THM-3544-planar-keller-target-pencil-total-degree-six-floor.md)
  gives an all-target-pencil degree-six floor,
  [THM-3550](../01-canon/theorems/THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor.md)
  excludes prime pencil degrees and the height-six/seven box, and
  [THM-3548](../01-canon/theorems/THM-3548-planar-keller-conductance-shadow-gates.md)
  records the exact conductance-side filters,
  [THM-3553](../01-canon/theorems/THM-3553-fixed-cubic-keller-map-polynomial-graph-section-no-go.md)
  excludes every displayed polynomial graph for the fixed ambient cubic map,
  while
  [THM-3545](../01-canon/theorems/THM-3545-catalan-self-intersection-keller-thickening-boundary.md)
  gives an exact formal collision with a nonterminating Catalan obstruction
  to polynomiality,
  [THM-3555](../01-canon/theorems/THM-3555-catalan-thickening-universal-cubic-root-cover.md)
  identifies its algebraic polynomialization with the universal cubic
  marked-root cover and rules out every correction fixing its ramification
  line pointwise,
  [THM-3557](../01-canon/theorems/THM-3557-low-width-mixed-catalan-thickening-no-go.md)
  closes widths one and two and the width-three coefficient box through
  degree five,
  [THM-3554](../01-canon/theorems/THM-3554-punctured-kummer-collision-surface-normal-form.md)
  gives an exact etale collision whose puncture cannot be filled while
  retaining its quadratic function field, and
  [THM-3556](../01-canon/theorems/THM-3556-cusp-square-packet-marked-root-kummer-owner.md)
  identifies the cusp packet's exact marked-root/Kummer factorization and
  nonlinear decomposable-projection target.
- Least-used sidecar: the normalization/order-index distinction in the fixed
  higher-dimensional Keller packets.  It is not a planar theorem, but it is
  the correct hostile test for proposed inverse-cover constructions.

### Six live concepts

1. leading-form Euclidean divisibility and root multiplicity;
2. equal-sum coefficient fibres and magnetic holonomy;
3. Gaussian--Hadamard row transport and its index-two parity sidecar;
4. cusp-square discriminant packets and nonlinear projection;
5. coordinate-hypersurface descent from an ambient Keller collision;
6. connected cubic root covers and ramification-to-infinity surgery.

The session changed the board as follows.  The row transport became an exact
map, but failed the ordinary Keller predicate by a Laurent Jacobian factor.
The resistor analogy became exact at leading dephasing order, but the gain
holonomy became the mandatory sidecar.  The cusp construction lost every
obvious two-coordinate projection, but gained an everywhere-immersive
four-coordinate packet.  The Euclidean lane gained an explicit
root-multiplicity divisor rather than another untyped recurrence guess.
The Catalan tail became a coordinate symptom rather than an isolated series:
its algebraic coordinate is the ramified universal cubic root cover, sharply
locating the missing operation as moving finite branching to nonproper escape.

## 3. The quantum-to-resistor statement, with constants and boundaries

Let `G=(V,E)` be a finite simple undirected graph, with `c_xy=c_yx>0`.  Put

```text
H_xy=sqrt(c_xy) exp(i theta_xy),  H_yx=conj(H_xy),  H_xx=0.
```

Use site projectors `P_x=|x><x|` and the Lindblad equation

```text
rho_dot=-i[H,rho]
        +lambda sum_x(P_x rho P_x - (1/2){P_x,rho}).       (3.1)
```

Under this convention each off-diagonal matrix unit decays at rate exactly
`lambda`.  Let `Pi` project a density matrix to its diagonal, put
`K=-i[H,-]`, and write `Q=1-Pi`.  On diagonal and off-diagonal blocks,

```text
p_dot=Aq,
q_dot=Bp+Cq-lambda q,
A=Pi K Q,  B=Q K Pi,  C=Q K Q.                           (3.2)
```

Eliminating the fast variable gives the bulk slow manifold

```text
q=lambda^-1 Bp + lambda^-2 C Bp + O(lambda^-3),
p_dot=lambda^-1 Pi K^2 Pi p
     +lambda^-2 Pi K^3 Pi p+O(lambda^-3).                (3.3)
```

An exact entry calculation yields

```text
(Pi K^2 Pi p)_x=2 sum_y |H_xy|^2(p_y-p_x)=2(L_c p)_x.    (3.4)
```

Therefore

```text
p^(lambda)(lambda*tau/2) -> exp(tau L_c)p(0),
(L_c p)_x=sum_y c_xy(p_y-p_x).                           (3.5)
```

The convergence is uniform on compact slow-time intervals in this finite
dimensional setting.  Arbitrary initial coherences contribute an
`O(lambda^-1)` initial slip and do not change the limiting populations.

The convention-independent rate rule is

```text
effective physical jump rate
  =2|H_xy|^2/(coherence-decay rate on xy).                (3.6)
```

This is compatible with the established Haken--Strobl site-dephasing picture;
the literature also emphasizes that very large dephasing freezes motion on
unrescaled physical time ([Bressanini--Benedetti--Paris](https://arxiv.org/abs/2204.01836)).

### Electrical dictionary

Equation `(3.5)` is exactly the variable-speed continuous-time walk of the
conductance network.  Thus

```text
L_c h=0  <=>  sum_y c_xy(h_y-h_x)=0,
```

which is Kirchhoff current balance.  Dirichlet energy, equilibrium potentials,
hitting probabilities, Green functions, effective resistance, and commute
times follow without another approximation.  On a connected component `C`,

```text
E_a T_b + E_b T_a = |V(C)| R_eff(a,b).
```

On a disconnected graph each component conserves its mass, isolated vertices
are stationary, and intercomponent resistance/commute time is infinite.

### The first phase sidecar

For an oriented triangle `x->y->z->x`, set

```text
Phi_xyz=Im(H_xy H_yz H_zx).
```

The `lambda^-2` physical correction in `(3.3)` is

```text
(T_theta p)_x=6 Phi_xyz(p_z-p_y),
(T_theta p)_y=6 Phi_xyz(p_x-p_z),
(T_theta p)_z=6 Phi_xyz(p_y-p_x).                        (3.7)
```

It is a skew circulation, not a correction to symmetric conductance.  After
the slow rescaling it appears with coefficient `1/(2lambda)`.  Consequently:

- trees are phase-independent for every `lambda`, because all phases are
  removed by vertex gauge;
- triangle-free graphs have no phase term at the first slow correction;
- bipartite graphs first see magnetic phase through four-cycle holonomy one
  order later;
- a shortest magnetic cycle of length `ell` can first enter at physical
  order `lambda^(-(ell-1))`, subject to flux cancellation.

### Necessary caveats

- Hermiticity requires opposite directed phases to be conjugate.
- With parallel coherent edges, amplitudes add before squaring.  The original
  conductances are recovered only after aggregating the hopping appropriately
  or giving the channels orthogonal environmental labels.
- For nonuniform vertex dephasing `lambda eta_x`, edge `xy` has coherence
  decay `lambda Gamma_xy`, `Gamma_xy=(eta_x+eta_y)/2`, and effective rate
  `2c_xy/(lambda Gamma_xy)`.  A single clock recovers the original network
  iff `Gamma_xy` is constant on occupied edges.  On a connected non-bipartite
  graph this forces uniform `eta`; on a bipartite graph the two color classes
  may carry complementary constants.
- “All phases disappear” means the singular leading limit, not every finite
  correction.

## 4. The exact JC coefficient gain graph

Write

```text
P=sum_a p_a x^a1 y^a2,    Q=sum_u q_u x^u1 y^u2.
```

Then

```text
Jac(P,Q)
 =sum_(a,u) p_a q_u det(a,u)
    x^(a1+u1-1)y^(a2+u2-1).                              (4.1)
```

Fix `P` and a prospective finite support `U` for `Q`.  Make a bipartite gain
graph with input vertices `u`, output-row vertices `w`, and edge

```text
A_(w,u)=p_a det(a,u),  where a+u=w+(1,1).                (4.2)
```

The finite Keller response equation is exactly

```text
Aq=kappa e_(0,0).                                        (4.3)
```

Thus the user's instruction to group pairs by equal sums is literally the
coefficient-row organization of the planar Jacobian, but in the exponent
lattice and with determinant weights.

Hermitianize the response matrix:

```text
H_A = [0  A*]
      [A   0].                                           (4.4)
```

Strong dephasing turns `(4.4)` into the bipartite resistor graph with
conductances `|A_(w,u)|^2`.  It preserves:

- zero/nonzero support;
- magnitudes;
- matching data;
- unique-channel/leaf obstructions.

It destroys:

- determinant signs and coefficient phases;
- cancellation around equal-sum fibres;
- rank and target-range information on cyclic support.

For a connected support graph, magnitudes plus a basis of cycle holonomies
reconstruct the matrix up to row/column phase gauge.  The number of missing
independent phases is

```text
beta_1=|E|-|V|+1.                                        (4.5)
```

If the support is a forest, all phases are gauge-removable.  Every nonzero
minor has at most one matching term, so its rank equals the maximum matching
number for every assignment of nonzero gains.  Cycles are exactly where the
resistor shadow ceases to determine linear algebra.

### Exact hostile: identical resistor shadow, different Jacobian cancellation

For

```text
P_+=x+(1/2)x^2y^2+(1/3)x^3y+(1/3)xy^3,
P_-=x+(1/2)x^2y^2+(1/3)x^3y-(1/3)xy^3,
Q=q_x x+q_y y,
```

the response block on rows `x^2y,xy^2` is

```text
A_+=[-1  1]       A_-=[-1  1].
    [-1  1]           [ 1  1]
```

Every squared edge magnitude is one.  Nevertheless `rank A_+=1`, with
`A_+(1,1)^T=0`, while `rank A_-=2`.  The square holonomies have opposite
sign.  This is an actual Jacobian coefficient block, not a free matrix
example.  It refutes any attempt to decide the Keller equations using only
the resistor shadow.

### Leaf peeling

Every nonconstant output row in `(4.3)` has target zero.  If it has one live
incident edge, its adjacent coefficient must vanish.  Delete that input and
repeat.  If peeling removes every path to the constant row, the proposed
support admits no mate.  This recovers the mechanism of THM-2045: its highest
weighted output is a unique leaf with a nonzero gain.

### Intrinsic conductance gates

[THM-3548](../01-canon/theorems/THM-3548-planar-keller-conductance-shadow-gates.md)
quantifies four things that remain visible without pretending the shadow is
complete.

First, for `A=DF(z)`, let `C=(|A_ij|^2)` and `T=sum C_ij`.  Then

```text
|det C| <= |kappa|T/2,
dist_F(C/T,{rank<=1}) <= |kappa|/T.
```

Along a sequence with `T->infinity`, either the two determinant matchings
become equal in magnitude and phase -- a dark plaquette -- or one large
channel is paired with a vanishing opposite entry.  The Wilson phase still
distinguishes determinant zero from determinant one for identical `C`.  More
exactly, with `r=|ad|`, `s=|bc|`, and
`Phi=arg(ad*conj(bc))`,

```text
|kappa|^2=(r-s)^2+4rs sin^2(Phi/2).
```

Second, put `z_(a,u)=det(a,u)p_aq_u`, `c_e=|z_e|^2`, and let a nonconstant
equal-sum fibre have `m` live channels and energy `E=sum c_e`.  Polygon
closure gives

```text
max_e sqrt(c_e) <= sum_(f!=e)sqrt(c_f),
c_max <= ((m-1)/m)E.
```

Here `c_e=|A_(w,u)q_u|^2` is a channel-flow intensity, not the bare dephased
edge conductance `|A_(w,u)|^2`; it cannot prune `P` without information about
the prospective mate.

But the fibres cannot be solved independently: the unweighted matrix
`W_(a,u)=p_aq_u` must lie on the Segre rank-one locus, so all of its `2x2`
minors, equivalently all alternating cycle products in the full coefficient
table, vanish.  This Segre graph has `supp(P)` and `supp(Q)` as its two vertex
classes; it is not the response graph of `Q`-monomials versus output rows in
`(4.2)`.  The two cycle sidecars must not be conflated.

Third, the pointwise energy-minimizing target shear is controlled by

```text
beta=<grad Q,grad P>/||grad P||^2.
```

If `beta` were a polynomial in `P`, a legal target shear would make the map
affine.  Hence every counterexample has global shear frustration
`beta notin C[P]`.  Finally, any supplied Puiseux escape branch
`x(t)~b t^(-r)`, `F(x(t))=a+q t^m+...` obeys the conditional bound
`m<=r(d-2)`, where `d=max(deg P,deg Q)`.  The existence of such a displayed
branch is an external curve-selection input, not part of the theorem.

The right low-complexity program is therefore not “solve JC by resistance.”
It is:

```text
gain graph -> resistor shadow for pruning
           + cycle holonomy sidecar for cancellation and range.             (4.6)
```

## 5. The two square arrays are one Gaussian--Hadamard object

Define the triangular fixed-second-coordinate row

```text
B_t={j^2+(t+1)^2 : 1<=j<=t}.                            (5.1)
```

The user's displayed rows are exactly `B_1,B_2,...`:

```text
{5}, {10,13}, {17,20,25}, {26,29,34,41}, ...
```

Now apply

```text
M(j,t+1)=(t+1-j,t+1+j),
M=[-1 1; 1 1],  det M=-2,  M^T M=2I.                   (5.2)
```

This bijects `1<=j<=t` with all strict pairs on equal-sum row
`n=2t+1`, because the output coordinates sum to `2t+2=n+1`.  Therefore

```text
{norms on odd sum-row 2t+1}=2B_t.                       (5.3)
```

Up to a sign gauge, `M` is multiplication by `1+i` in the Gaussian integers.
Its determinant, norm scaling, and parity/content tax are the exact
ramified-two operation audited in
[THM-3336](../01-canon/theorems/THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md).

For even sum-row `n=2t`, the same chart lives on the half-integer coset:

```text
(t+1-j,t+j)=M(j-1/2,t+1/2).                             (5.4)
```

This explains the parity split, but `(5.4)` is not a polynomial exponent
lattice before applying `M`.

### Counts and the correct filler semantics

On row `n`, strict positive pairs `x<y`, `x+y=n+1`, number

```text
floor(n/2).                                              (5.5)
```

Two orientations give `n` entries for even `n` and `n-1` for odd `n`.
For odd `n=2t+1`, the missing swap-fixed support point is

```text
(t+1,t+1),  norm=2(t+1)^2.                              (5.6)
```

It is exactly the image under `M` of the omitted boundary point `(0,t+1)`.
Thus `(5.6)` is the canonical filler when the cells represent actual exponent
vectors or oriented pairs.

The proposed scalar filler has a different, equally clean role.  Put

```text
C_t=t^2+(t+1)^2,    Q_t=2C_t+1=(2t+1)^2+2.              (5.7)
```

`C_t` is the last element of `B_t`; by `(5.3)`, `2C_t` is the largest odd-row
norm.  Hence `Q_t` is exactly the first integer above the row:

```text
Q_t: 3,11,27,51,83,123,...                              (5.8)
```

This is the audited sequence of
[THM-3341](../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md).
Use `Q_t` as a collision-free scalar sentinel; do not pretend it is the
missing support point.

The count generating functions make the completion exact:

```text
R(z)=sum_(n>=1) floor(n/2) z^n=z^2/((1-z)(1-z^2)),
2R(z)+z/(1-z^2)=z/(1-z)^2=sum_(n>=1)n z^n.              (5.9)
```

The filler-value series is

```text
sum_(t>=0) Q_t z^(2t+1)
 =z(3+2z^2+3z^4)/(1-z^2)^3.                             (5.10)
```

### `4^a(8b+7)`, the filler, and the exact nearby ladder

Legendre's three-square obstruction is

```text
N is not a sum of three integer squares
  <=> N=4^a(8b+7).                                      (5.11)
```

For odd `n=2t+1`, `n^2=8T_t+1`, so

```text
n^2+1=8T_t+2=n^2+1^2,                       two squares,
n^2+2=8T_t+3=n^2+1^2+1^2,                   three squares,
n^2+6=8T_t+7=n^2+2^2+1^2+1^2,              four needed. (5.12)
```

The user's filler is the middle line.  It is `3 mod 8`, hence allowed by the
three-square theorem, and `3 mod 4`, hence never a sum of two squares.  The
deliberate forbidden companion is only four away:

```text
Q_t+4=n^2+6=8T_t+7,
4^a(Q_t+4)=4^a(8T_t+7).                                  (5.13)
```

Modulo eight alone detects only the `a=0` sheet: the residues are `7,4,0`
for `a=0,1,>=2`.  The correct test repeatedly divides by four and then checks
for residue seven.

Summing all norms on odd row `2t+1` gives another appearance of `8t+7`:

```text
S_t=2 sum_(j=1)^t ((t+1)^2+j^2)
   =t(t+1)(8t+7)/3.                                     (5.14)
```

The visible factor is not itself a Legendre certificate.  For example,
`S_3=124=4*31` is forbidden, while `S_4=260=4*65` is not.  The prefactor's
two-adic valuation and odd residue remain load-bearing.

### Representations versus distinct values

The triangular array enumerates every positive pair `x<y` once, not every
numeric sum once.  The first duplicate is

```text
65=4^2+7^2=1^2+8^2,                                    (5.15)
```

in consecutive rows 6 and 7.  If `r_2(N)` counts all ordered signed
two-square representations, the number of unordered positive distinct
representations is

```text
D(N)=(r_2(N)-4*1_(N square)-4*1_(N=2 square))/8,
r_2(N)=4(d_1(N)-d_3(N)).                                (5.16)
```

Thus the array has fixed row geometry at the cost of numeric repetitions.
A distinct-value array must choose one representation per value and loses
that geometry.

### Why the beautiful map is not a JC transform

The exponent map `(5.2)` corresponds on the torus to

```text
U=Y/X,  V=XY,  Jac(U,V)=-2Y/X.                          (5.17)
```

Therefore

```text
Jac(P(U,V),Q(U,V))
 =(-2Y/X) Jac(P,Q)(U,V).                                (5.18)
```

It preserves additive fibres, scales exponent determinants by `-2`, and is
natural for a log-canonical toric bracket.  It does **not** preserve an
ordinary constant planar Jacobian.  Its index-two parity sheet, negative
exponent, and missing axes are not cosmetic sidecars.

## 6. Squeezing a hypothetical planar counterexample

Take a representative whose sorted degree pair is minimal under polynomial
target automorphisms.  Let

```text
P_n=cH^a,  Q_m=c'H^b,
n=ga,      m=gb,
g=gcd(n,m), gcd(a,b)=1.                                 (6.1)
```

If `a=1`, a target shear cancels the leading term of `Q` by a power of `P`;
if `b=1`, shear the other coordinate.  Hence a target-shear-minimal
counterexample has

```text
a,b>=2.                                                  (6.2)
```

The cited one-place-at-infinity exclusion gives at least two distinct roots
of `H`.  A stronger independent degree restriction, proved in the
Guccione--Guccione--Valqui shape analysis, gives

```text
g=gcd(n,m)>=16,                 g!=2p for every prime p.
```

([primary paper](https://arxiv.org/abs/1401.1784)).  The corrected
Appelgate--Onishi theorem says a
counterexample cannot have a component degree equal to a product of at most
two primes; its proof history includes the later Nowicki--Nakai repairs
([original](https://doi.org/10.1016/0022-4049(85)90099-4),
[lemma repair](https://doi.org/10.1016/0022-4049(88)90069-2)).  Complete any
nonzero pencil member `R=sP+tQ` to a target basis and apply the theorem there.
Consequently **every** pencil degree has at least three prime factors counted
with multiplicity, and is therefore at least eight.

The internal pencil theorems give a complementary structural statement.
THM-3544 and THM-3550 prove that every pencil degree is composite and at
least six, and that one can choose a reduced target basis `(R_n,S_m)` with

```text
6<=n<m,  n and m composite,                             (6.2a)
```

such that `[R]` is the unique degree-`n` direction and every other pencil
direction has degree `m`.  The internal gates alone give `m>=8`; combined
with Appelgate--Onishi they give `n>=8`.  Moh's computation, applied to this
reduced basis, gives

```text
m>100.                                                  (6.2b)
```

([Moh's summary](https://www.math.purdue.edu/~ttm/jacobian.html)).  Thus the
first arithmetic search is not a displayed low-degree pair: it is a unique
low pencil direction of composite three-prime-factor degree below a height
strictly above 100.

The equality boundary of the internal height theorem is still informative.
If one ignores the stronger cited exclusions and sets `m=8`, the reduced
pair must have degrees `(6,8)` and leading forms either `(L^6,L^8)` or,
after `GL_2`, `((xy)^3,(xy)^4)`.  The one-place and Appelgate--Onishi gates
kill these as counterexamples: the one-place gate removes the linear-base
cell, while Appelgate--Onishi removes both because of the degree-six member.
They remain exact hostile controls for any claimed degree-floor argument.
THM-3550 states this with a power-free base.  In the maximal common-base
normalization of `(6.1)`, the first cell has `H_max=L^2`, exponents `(3,4)`,
and one root direction; the second has `H_max=xy`, the same exponents, and
two root directions.  These normalizations must not be identified.

Combining only the Appelgate--Onishi/Nagata, Moh, and 2014 gcd gates gives a
useful but superseded finite screen.  Its first surviving height is

```text
m=105,
(n,m;g;n/g,m/g) in
{(42,105;21;2,5), (63,105;21;3,5),
 (70,105;35;2,3), (84,105;21;4,5)}.                   (6.2c)
```

This is **FINITE-EXACT arithmetic under that named weaker sieve**, not a live
frontier.  The cited 2022 Guccione--Guccione--Horruitiner--Valqui
classification of degree pairs below height `125` excludes all four cells.
They are retained as hostile controls because the new divisor already taxes
three of them.  After swapping outputs to put the larger leading exponent
first, those taxes are

```text
(63,105): rad(H) | F_62,       H^2 rad(H) | G_104;
(70,105): rad(H) | F_69,       H   rad(H) | G_104;
(84,105): H rad(H) | F_83,     H^2 rad(H) | G_104.
```

The omitted stronger input is a cited preprint
([primary](https://arxiv.org/abs/2204.14178)), not proved here.  Up to
transpose, it leaves exactly one reduced pair below height `125`:

```text
(n,m;g;n/g,m/g)=(72,108;36;2,3).                      (6.2d)
```

This is now the live first degree passport.  Put `F=F_72`, `G=G_108` and
swap the outputs when applying the `b<a<2b` divisor theorem.  Then
`(a,b)=(3,2)`, `k=1`, and

```text
rad(H) | F_71,                    H rad(H) | G_107.    (6.2e)
```

In particular, a squarefree degree-36 base forces `H|F_71` and
`H^2|G_107`.  Then `F_71=H C_35`, so the entire subleading pair is controlled
by only `36` coefficients instead of the original `72`; its codimension in
the `108+72`-dimensional pair of binary-form rows is `108+36=144`.  The
former height-`105` wording is corrected in MISTAKE-427.

For completeness, within the obsolete height-`105` screen the last three
cells enter the chamber `b<a<2b`, with `k=1,1,3` respectively.
The remaining `(42,105)` cell has exponent pair `(2,5)` and lies in the next
chamber
`2b<a<3b`; this explains why the weaker hostile screen still tests both the
proved first chamber and its failure boundary.  THM-3550 is now labelled
`PROVED / INDEPENDENTLY HOSTILE-AUDITED`; the stronger degree lists remain
separately cited inputs.

### Subleading codimension

Assume `a>b`.  Repaired THM-3025 gives

```text
P_(n-1)=kappa H^(a-b) Q_(m-1),
kappa=ca/(c'b).                                         (6.3)
```

The admissible pair of subleading forms is parametrized by the `m`
coefficients of `Q_(m-1)`, inside an ambient space of dimension `n+m`.
Thus `(6.3)` has codimension `n=max(n,m)` at this layer.

The proved first Euclidean block in THM-3016/HYP-9070 gives, only through the
first partial quotient,

```text
H^(a-ib) | P_(n-i),  0<=i<=floor(a/b).                  (6.4)
```

At layer `i`, the quotient degree is `i(m-1)`, so divisibility alone imposes
codimension `n-im`.  No indefinite recurrence is asserted: a new face or
renewal is required when this block ends.

### New root-multiplicity divisor in `b<a<2b`

Normalize the leading constants for readability and put

```text
R=Q_(m-1),  D(F)=Jac(H,F),  k=2b-a>=1.                  (6.5)
```

The degree-`n+m-4` Keller equation, after using `(6.3)`, is

```text
aH^(a-1)D(Q_(m-2))
 +(a/b)(a-b)H^(a-b-1)R D(R)
 -bH^(b-1)D(P_(n-2))=0.                                (6.6)
```

After factoring `H^(a-b-1)`, both outer terms are divisible by `H^k`.
Therefore

```text
H^k | R D(R).                                           (6.7)
```

Fix `L^e || H` and put `u=v_L(R)`.  Since `H` has at least two distinct
roots, `0<e<g`.  In coordinates `L=x`, the exact first coefficient of
`D(R)` at `x`-order `e+u-1` is

```text
unit * [g(eb-u)-e].                                     (6.8)
```

It cannot vanish, because it is congruent to `-e mod g`.  Hence

```text
v_L(D(R))=e+u-1.                                        (6.9)
```

Combining `(6.7)` and `(6.9)` gives

```text
2u+e-1>=ek,
u>=ceil((e(k-1)+1)/2).                                  (6.10)
```

Thus

```text
D_k(H)=prod_(L^e || H)L^ceil((e(k-1)+1)/2) | Q_(m-1),  (6.11)
H^(a-b)D_k(H) | P_(n-1).                                (6.12)
```

Let `delta_k(H)=deg D_k(H)`.  The subleading dimension drops from `m=gb`
to `m-delta_k(H)`, and its ambient codimension rises from `n` to
`n+delta_k(H)`.

For squarefree `H`, `(6.11)` becomes

```text
H^ceil(k/2) | Q_(m-1).                                  (6.13)
```

Examples:

- `(a,b)=(3,2)`: `k=1`, so `H|Q_(m-1)`;
- `(5,4)`: `k=3`, so squarefree `H^2|Q_(m-1)`;
- consecutive `(b+1,b)`: the forced squarefree power is
  `ceil((b-1)/2)`.

Repeated roots do not create a resonance escape; they increase the
multiplicity tax according to `(6.10)`.

The translation boundary is also exact.  In the `(3,2)` squarefree chamber,
`Q_(2g-1)=H C_(g-1)`.  Source translation changes `C` only by the
two-dimensional span of `H_x,H_y`.  When `g=2` this kills every `C`; for
`g>2`, at least `g-2` quotient moduli remain.  This is why the compact
quadratic certificate below closes and its higher-`g` analogues do not.

### Exact `(6,4)`, `H=xy` jet closure

Set `P_6=H^3`, `Q_4=H^2`, `H=xy`.  The proof uses only source `GL_2`, source
translation, target scaling/shear, and homogeneous Keller equations.

1. Repaired THM-3016 gives
   `P_5=(3/2)H Q_3`.
2. The next equation implies `H|Q_3 Jac(H,Q_3)`.  Restriction to `x=0` and
   `y=0` gives nonzero squares of the endpoint coefficients, so
   `Q_3=H L` for a linear form `L`.
3. A two-parameter source translation kills `Q_3` and `P_5` simultaneously.
4. Write `Q_2=Ax^2+Bxy+Cy^2`.  The next equation gives

   ```text
   P_4=(3/2)H(Ax^2+Cy^2)+dH^2.
   ```

5. The degree-four equation forces `H|Jac(P_4,Q_2)`.  Its axis values are
   `3C^2y^4` and `-3A^2x^4`, hence `A=C=0`.
6. A target shear kills `dH^2`; then `P_2=eH`.
7. The next odd layer gives `P_3=(3/2)H Q_1`.
8. Finally

   ```text
   Jac(P_3,Q_1)=-(3/2)(alpha x-beta y)(alpha x+beta y)=0,
   ```

   so `Q_1=0`, contradicting the nonzero constant Jacobian row.

This is a useful compact template for the first Euclidean pair `(3,2)`, but
not a novel degree exclusion.  The external degree theorems already subsume
it.

## 7. Counterexample architecture atlas

### A. Exact reflection symmetry: a beautiful factory that cannot fire

In coordinates `u=x+y`, `v=x-y`, `w=v^2`, a map commuting with reflection
has the form

```text
U=A(u,w),  V=vB(u,w).                                   (7.1)
```

Its Keller equation is the first-order polynomial PDE

```text
A_u(B+2wB_w)-2wA_wB_u=c.                                (7.2)
```

After affine normalization, `A(u,0)=u`, `B(u,0)=1`.  If a nonconstant `B`
solved `(7.2)`, it would vanish somewhere with `w!=0`; the two source points
`(u,+sqrt(w))` and `(u,-sqrt(w))` would then have the same image.  Thus a
nontrivial polynomial solution would be an immediate planar counterexample.

Three low transverse boxes close elementarily:

- `deg_w(A-u),deg_w B <=(1,1)` reduces to
  `2aa''=3(a')^2`, whose leading term is impossible for nonconstant
  polynomial `a`;
- boxes `(1,2)` and `(2,1)` have leading obstructions proportional to
  `N(N+2)(N+4)` and `N^2(N+2)(2N+1)`.

More importantly, this entire exact-symmetry lane is already **CITED
REFUTED**: Miyanishi proves that a planar etale endomorphism commuting with an
effective finite group of even order is an automorphism
([paper](https://arxiv.org/abs/2110.06709)).  A counterexample search must
break every such exact involution.  The PDE remains useful as a near-symmetry
hostile and as a guide to the first asymmetric defect layer.

### B. Norm-balanced gain-cycle tiling: OPEN

For every `r>=3`, let

```text
u=(1,8),       u'=(4,7),
v=(r+2,3r+1),  v'=(r-1,3r+2).                           (7.3)
```

Then

```text
||u||^2=||u'||^2=65,
||v||^2=||v'||^2=10r^2+10r+5,
u+v=u'+v',
det(u,v)=-5(r+3),  det(u',v')=+5(r+3).                  (7.4)
```

One equal-sum row can therefore cancel with matched coefficient products.
For `r=4`, the norm pair is `65/205` and the cancelling determinants are
`+/-35`.

The hostile is immediate: the complete product support also creates cross
leaks `u+v'` and `u'+v`, with singleton determinants `+/-5(r-2)`.  The
four-term gadget alone is impossible.  The open architecture is to tile
these norm-balanced parallelograms until every leak lands in another
multi-edge fibre, then solve the determinant-ratio holonomy equations.  All
cross edges must be retained; a hand-selected cancellation square is not a
polynomial support.

### C. Ambient coordinate-hypersurface descent: OPEN and decisive

[THM-3546](../01-canon/theorems/THM-3546-invariant-graph-keller-descent-criterion.md)
proves a clean general criterion.  If a constant-Jacobian map in `n+1`
variables sends one polynomial coordinate graph scheme-theoretically into
another -- equivalently,
`F_w(x,h(x))=g(F_y(x,h(x)))` as a polynomial identity -- its restriction is
Keller in `n` variables; any collision on the graph descends.  Over an
infinite field, containment on all rational points implies this identity;
over a finite field, pointwise containment alone is not enough.

For an ambient three-variable collision, the search becomes:

1. find a source coordinate polynomial `rho_s` vanishing on two colliding
   preimages;
2. find a target coordinate `rho_t` vanishing at their common image;
3. solve `rho_s | rho_t(F)` exactly;
4. verify both are polynomial coordinates, not merely smooth hypersurfaces.

Passing these gates would produce a planar counterexample with no second
Jacobian search.  Interpolation is cheap; divisibility and coordinate status
are load-bearing.

For the fixed sporadic three-variable Keller map, the displayed polynomial-
graph lane now closes in **every** degree.  If `H` is the top homogeneous
form of a nonconstant `h`, the unique top tangential Jacobian row is

```text
-3x^5y^4 H(3H+xH_x),
```

and the Euler-weight operator `3+x partial_x` is injective in characteristic
zero.  For constant `h=c`, the top row is

```text
-9x^3y^4(cx+3y)(cx+2y),
```

including the nonzero hostile `-54x^3y^6` at `c=0`.  Hence the restricted
tangential Jacobian is never constant, so THM-3546 forbids any displayed
source graph `z=h(x,y)` from landing scheme-theoretically in a displayed
polynomial target graph.  This is the proved
[THM-3553](../01-canon/theorems/THM-3553-fixed-cubic-keller-map-polynomial-graph-section-no-go.md).
Graphs after nonlinear ambient coordinate changes and nongraph coordinate
hypersurfaces remain open; the plane `x=0` is the sharp hostile because it
does descend to a polynomial plane automorphism.

The categorical torus quotient is a proved hostile, not a solution:
[THM-3543](../01-canon/theorems/THM-3543-torus-quotient-ramification-square-no-go.md)
computes an exact quotient Jacobian `2(2-3v-t)^2`.  The quotient retains the
collision only by contracting a divisor and acquiring ramification.  Graph
restriction and quotient forgetting are opposite operations.

[THM-3549](../01-canon/theorems/THM-3549-torus-quotient-correction-no-go.md)
then tests the literal repair.  In collision coordinates the seed has
`Jac(R,S)=-2w^2`.  No one-sided polynomial correction works in any degree;
arbitrary simultaneous corrections of total degree at most two are
inconsistent, degree-three corrections force the constant Jacobian to zero,
and corrections affine in `w` fail for every degree in the other variable.
For a collision-preserving Keller repair, the first surviving transverse box
has sorted `w`-degrees at least `(4,5)`.  On its boundary the top coefficients
must already be

```text
a(u)=c h(u)^4,                 b(u)=d h(u)^5.
```

Thus “subtract the conductor” is not a low-order perturbation.  The open
cell requires simultaneous high-transverse corrections to both outputs,
removal of the contracted line without losing the two-point collision, and
global cancellation of every positive `w`-row.

There is a particularly sharp clean-total-degree boundary.  If
`deg h=e-1` and these `(4,5)` transverse tops dominate ordinary total degree,
then

```text
(deg P,deg Q)=(4e,5e),       H_top=u^(e-1)w.
```

The weaker gcd, Appelgate--Onishi, and Moh gates permitted `e=21`, but the
2022 sub-`125` classification excludes that `(84,105)` hostile.  For every
`e<25`, the height `5e` is below `125`, while the ratio `4:5` cannot equal
the unique surviving ratio `72:108=2:3`.  Hence the first row-dominant box
not excluded by the cited list is

```text
(deg P,deg Q)=(100,125),       e=25,       H=u^24w.
```

Repaired subleading rigidity and `(6.11)` give

```text
Q_124=(5d/(4c))H P_99,
D_3(H)=H rad(H)=u^25w^2 | P_99,
u^49w^3 | Q_124.
```

Thus `P_99` is parametrized by a degree-72 quotient: `73` homogeneous
coefficients rather than `100`.  This entire conclusion is conditional on
the displayed tops controlling ordinary total degree; a lower `w`-row with
larger `u`-degree creates a second Newton scale and must be retyped
separately.  MISTAKE-427 records why `(84,105)` is now only a hostile.

### D. Catalan self-intersection thickening and its universal cubic cover

Consider the separated formal ansatz

```text
P(v,w)=v^2+A(w),
Q(v,w)=v^3-v+vB(w),            A(0)=B(0)=0.
```

Direct coefficient comparison gives

```text
Jac(P,Q)=kappa
iff B=(3/2)A and A-(3/4)A^2=kappa*w.
```

The unique normalized branch is

```text
A=(2/3)(1-sqrt(1-3*kappa*w))
 = sum_(n>=1) C_(n-1)(3/4)^(n-1)kappa^n w^n.
```

It has the transverse double collision `(+/-1,0)->(1,0)` and constant
Jacobian as a formal, or local holomorphic, map.  Every Catalan coefficient
is nonzero in characteristic zero, so the tail never terminates and the map
is not polynomial.  This is
[THM-3545](../01-canon/theorems/THM-3545-catalan-self-intersection-keller-thickening-boundary.md),
not a counterexample.

The square root is not merely a generating-function artifact.  Adjoin

```text
r^2=1-3*kappa*w.
```

Then the Catalan map polynomializes to

```text
H(v,r)=(v^2+(2/3)(1-r), v^3-vr),
Jac_(v,r)(H)=-(2/3)r,             dr/dw=-3*kappa/(2r).
```

The two factors cancel under the chain rule, explaining the formal constant
Jacobian.  With `t=v`, `p=2r-3v^2` and an affine target change, `H` becomes

```text
G(t,p)=(p,-t^3-pt).
```

This is the universal marked-root cover of `X^3+pX+q=0`: it is finite flat
of degree three, has Jacobian `R=p+3t^2`, and its discriminant pulls back as

```text
G^*(-4p^3-27q^2)=-(p+3t^2)^2(4p+3t^2).
```

The apparent Catalan pair is part of the simple three-point fibre
`(v,r)=(1,1),(-1,1),(0,-1/2)` over `(1,0)`.  Thus the square root hides the
third marked root while also cancelling a genuine finite ramification
factor.  This is the proved
[THM-3555](../01-canon/theorems/THM-3555-catalan-thickening-universal-cubic-root-cover.md).

It gives a strong surgery gate.  Any correction fixing `R=0` pointwise has
the form

```text
G_tilde=(p+R A, -t^3-pt+R B),
Jac(G_tilde)|_(R=0)=-6t(tA+B)|_(R=0),
```

so its Jacobian still vanishes at the cusp preimage `(t,p)=(0,0)`.  A viable
polynomial deformation must move the ramification line already at order
zero while separately preserving only a selected collision fibre.  This is
more restrictive than adding a higher normal jet.

Together with the torus-quotient hostile it forms a typed pincer:

- the quotient is polynomial and colliding, but its Jacobian has a square
  ramification factor;
- the Catalan thickening is colliding and locally Keller, but nonpolynomial.

The live architecture is to add genuinely mixed corrections that reroute the
Catalan tail, move the cubic ramification line, and replace finite branching
by nonproper escape at infinity.  Truncating the series is not enough: the
first omitted coefficient becomes an explicit nonconstant Jacobian defect.

There is a concrete positive first step.  Drop separation and write

```text
P=sum_(i=0)^N p_i(v)w^i,       p_0=v^2,
Q=sum_(j=0)^N q_j(v)w^j,       q_0=v^3-v.
```

The coefficient of `w^r` in the Keller equation is the exact ladder

```text
sum_(i+j=r+1) [j p_i' q_j-i p_i q_j']
   =kappa delta_(r,0).
```

For `kappa=1`, its first Bezout row has the polynomial solution
`p_1=1`, `q_1=3v/2`.  The open question is whether later mixed rows can
terminate while satisfying the terminal equations and retaining the boundary
collision.  Thus this is not merely “try more terms”: it is a finite
degree-by-degree construction ladder with a positive first row and the
separated Catalan solution as its nonterminating hostile control.

The first mixed boxes now close exactly.  For transverse width `N=1`, the
terminal equation forces the two top coefficients to be proportional and
the Bezout row becomes a polynomial unit impossibility.  At `N=2`, a
cap-free Wronskian and top-degree split closes all three top-row cases.  At
`N=3`, exact branching and saturated Groebner checks prove the affine
coefficient variety empty for coefficient caps `D=3,4,5`.  New square/cube
common-power types first appear at `D=6`, so that is the first **internally**
unclosed recurrence cell.  This is
[THM-3557](../01-canon/theorems/THM-3557-low-width-mixed-catalan-thickening-no-go.md).

The global degree gates move the actual counterexample frontier much farther.
If every `a_j,b_j` has degree at most `D`, both components have total degree
at most `D+N`.  Since the live first reduced height is `108`, any viable
width-three instance must have `D+3>=108`, hence `D>=105`.
Therefore the `D=6` branch is valuable for discovering a renewal mechanism,
but is already globally excluded as a planar counterexample cell.

### E. Punctured Kummer collision: the missing divisor is load-bearing

The curved collision surface of the fixed three-variable Keller map is

```text
C=2-3xy-x^2z=0.
```

With `s=x^(-1)` and `v=xy`, its coordinate ring is
`k[s,s^(-1),v]`, so the surface is `G_m x A^1`, not `A^2`.  Explicit Laurent
source and polynomial target automorphisms turn the restricted map into

```text
(s,b) -> (b,4s^2),                s!=0.
```

This is a finite etale double cover, and `s -> -s` is exactly the deck
collision.  Its natural affine completion has Jacobian `-8s` and ramifies at
the missing divisor `s=0`.  More strongly, an everywhere-etale `A^2` filling
cannot preserve the same quadratic extension: etale pullback of the reduced
branch divisor has valuation one, whereas `delta=4s^2` forces every such
valuation to be even.  This is the proved
[THM-3554](../01-canon/theorems/THM-3554-punctured-kummer-collision-surface-normal-form.md).

Here the image must be typed carefully.  As a morphism onto the punctured
target `U=A^2\setminus V(beta^2-16alpha)`, the cover is finite etale and
surjective, so its topological and algebraically-closed-point images are
`U`.  Its scheme-theoretic image closure in the ambient `A^2` is all of
`A^2`, while its `k`-rational-point image can be smaller by the square class
of `beta^2-16alpha`.  None of these distinctions changes the boundary
valuation obstruction.

The Catalan, quotient, and Kummer objects now form a three-way pincer:

```text
polynomial plane + collision     -> ramification square;
unit Jacobian + collision        -> infinite Catalan tail;
finite etale + collision         -> punctured Laurent plane.
```

The sheet structure sharpens the comparison.  The Kummer slice supplies an
etale but disconnected `1+2` cover, whereas the universal cubic has connected
generic `S_3` monodromy but finite discriminant ramification.  The desired
planar architecture must combine the missing halves: connected sheets with
the finite branch locus replaced by escape to infinity.

A viable deformation must change at least one load-bearing part of the exact
quadratic extension: mix transverse terms before filling, send the missing
divisor to infinity through a nonproper modification, or pass to higher-sheet
asymptotic monodromy with no finite branch component.

### F. Cusp-square packet and nonlinear `4 -> 2` projection: immersive projection seed

Start from

```text
T=y^2-6vU,
S=y^3-9vUy,
L=v^2(8vU-y^2).                                        (7.5)
```

For every polynomial `U(v,y)`, one has

```text
S^2=T^3+27LU^2.                                         (7.6)
```

This is an exact inverse-cubic owner, not just a discriminant mnemonic.  Put

```text
E(X)=LX^3+TX+2U.
```

Then

```text
E(X)=(vX+1)
     [v(8vU-y^2)X^2+(y^2-8vU)X+2U].                   (7.6a)
```

The marked root is `X=-1/v`, or projectively `[-1:v]`, and reaches infinity
exactly at `v=0`.  The quadratic factor has square class `-L` after removing
the square `(y/v)^2`, while the full discriminant is

```text
disc_X(E)=-4LS^2.                                      (7.6b)
```

Thus the packet contains one explicit escaping root and one quadratic
Kummer pair; `L` is the odd discriminant/infinity owner and `S^2` records
finite-root collisions.  This is the proved structural content of
[THM-3556](../01-canon/theorems/THM-3556-cusp-square-packet-marked-root-kummer-owner.md).
There is also a dual cubic visible entirely in the packet observables:

```text
Y^3-3TY+2S=(Y-y)[Y^2+yY+2(9vU-y^2)],
disc_Y=-(54U)^2L.
```

Here the finite source coordinate `y` is the marked root, and the residual
quadratic again has square class `-L`.  Any projection that hides the
escaping sheet must preserve this common resolvent sidecar, not just the
scalar cusp identity.
Every natural two-coordinate output is nevertheless obstructed:

```text
Jac(T,S)=54vU(U+vU_v),
v | Jac(L,T), Jac(L,U), Jac(L,S),
y | Jac(U,S)|_(v=0),
Jac(T,U)=-2(yU_v+3UU_y).                                (7.7)
```

The last expression cannot be a nonzero constant for polynomial `U`: expand
`U=sum_i a_i(y)v^i`; descending from the highest coefficients of `U^2`
forces all `a_i` constant, after which `yU_v=constant` is impossible.

There is nevertheless a positive packet.  Take

```text
U=1+y-y^2/2-(3/2)vy(y-3).                               (7.8)
```

For `G=(L,T,U,S):A^2->A^4`, exact rational Groebner reduction shows that the
six `2x2` Jacobian minors generate the unit ideal.  Thus `G` is immersive at
every affine source point.  Exact linear algebra also shows that no constant
linear combination of those minors is one, so no constant-linear projection
is Keller.

The open architecture is a **nonlinear** projection

```text
(A(L,T,U,S),B(L,T,U,S)):A^2->A^2                       (7.9)
```

with unit Jacobian, index-one coordinate orders, and a nonperiodic/wandering
image-prime curve so that a self-cover Euler obstruction does not fire.  This
does not constitute evidence against JC(2); it is a sharply structured
search box whose full four-coordinate differential obstruction is already
removed.

The legal projection coefficients are much more constrained than the unit
minor ideal suggests.  If `Z=(L,T,U,S)` and `A_i,B_i` denote derivatives in
the four packet coordinates, then

```text
Jac(A(Z),B(Z))=sum_(i<j)(A_iB_j-A_jB_i)(Z) M_ij.       (7.10)
```

The six coefficients must descend through `k[Z]`, be integrable, and be the
Pluecker coordinates of a decomposable two-form; in particular

```text
c_12 c_34-c_13 c_24+c_14 c_23=0.                      (7.11)
```

So the exact remaining problem is a descending, integrable, decomposable
Bezout certificate for the six minors, followed by a nonproperness and
nontrivial-fibre audit.  Arbitrary source coefficients exist; constant
coefficients fail; this typed middle class is open.

### G. Other attractive architectures that close

1. **Weighted suspension -- PROVED REFUTED.**  For
   `P=xA(x^r y^s)`, the only constant-producing weight sector of a mate is
   `Q_0=yf(w)`.  The top coefficient of

   ```text
   Jac(P,Q_0)=Af+swAf'+rwA'f
   ```

   is `a_d f_N(1+sN+rd)`, nonzero in characteristic zero.
2. **Two-mode reciprocal ansatz -- PROVED REFUTED.**  For
   `P=xf(xy)+yh(xy)`, `Q=xk(xy)+yg(xy)`, disjoint support sectors force two
   Wronskians to vanish; the remaining equation is a derivative of `tfg`
   and constancy reduces the map to linear.
3. **Homogeneous Hamiltonian perturbation -- REFUTED modulo the classical
   binary-Hessian lemma.**  `P=x+H_y`, `Q=y-H_x` gives
   `Jac(P,Q)=1+det Hess(H)`.  Vanishing homogeneous binary Hessian forces
   `H` to be a power of one linear form, hence a triangular shear.
4. **Bare cubic cusp -- PROVED REFUTED.**  The simplest parametrization has
   Jacobian proportional to the escape coordinate.  Equations `(7.5)`--`(7.7)`
   show that allowing arbitrary `U` still blocks all six natural coordinate
   pairs.
5. **Gaussian--Hadamard Laurent transform -- OPEN only as a toric ansatz.**
   Equation `(5.18)` changes a Keller constant into `-2Y/X`.  A lawful search
   would have to start with a Laurent pair of bracket `-X/(2Y)` whose image is
   polynomial and retains the index-two axis/parity sidecar.

## 8. What a counterexample must now look like

Within the proved/cited scopes used here, a target-shear-minimal planar
counterexample must satisfy all of the following:

1. both reduced leading exponents `a,b` are at least two and coprime;
2. the common leading form has at least two distinct root directions;
3. every nonzero target-pencil member has at least three prime factors with
   multiplicity, hence degree at least eight;
4. in a reduced target basis the pencil has one low degree `n` and one
   generic degree `m`, with `8<=n<m`, `m>100`, `gcd(n,m)>=16`, and gcd not
   twice a prime; the cited sub-`125` list makes `(72,108)` the unique first
   pair, with the radical taxes in `(6.2e)`;
5. the subleading correction `W` vanishes, by repaired THM-3025;
6. in `b<a<2b`, every root multiplicity obeys the extra divisor tax
   `(6.10)`--`(6.12)`;
7. the first Euclidean block must survive its layer-by-layer codimensions and
   then create a genuine new face rather than assume an endless scalar
   recurrence;
8. any sparse coefficient ansatz must survive leaf peeling, every cross leak,
   fibrewise polygon closure, and the global Segre cycle binomials;
9. a magnitude-only or norm-only shadow is insufficient: determinant signs,
   coefficient phases, cycle holonomies, and global shear frustration must be
   retained;
10. exact finite even symmetry is forbidden; a near-symmetric candidate needs
   a load-bearing asymmetric defect;
11. any supplied Puiseux escape branch has escape exponent
    `m_escape<=r(d-2)` and approaches either a dark rank-one conductance
    plaquette or a channel-degenerate regime;
12. displayed polynomial-graph descent from the fixed cubic Keller map is
    impossible in every degree, so ambient descent must use a nonlinear
    coordinate change or a nongraph coordinate hypersurface;
13. a torus-quotient repair preserving the inherited collision must alter
    both outputs and begin at sorted transverse degrees at least `(4,5)`;
    its first row-dominant `4:5` box not cited away is `(100,125)`;
14. a collision thickening must evade the separated Catalan tail and, in the
    universal cubic chart, move the ramification line at order zero rather
    than fix it through higher normal jets;
15. in the mixed Catalan ansatz, widths one and two are impossible, width
    three is internally closed through coefficient degree five, and the
    global height gate forces any viable width-three cap to be at least `105`;
16. the exact quadratic Kummer collision cannot be filled across its missing
    divisor while preserving etaleness and the same function field;
17. a cusp-packet projection must solve the descending, integrable,
    decomposable minor Bezout equation, retain its escaping marked root, and
    preserve normalization, coordinate-order index, Jelonek/nonproperness
    data, and hidden discriminant components.

Items 11--17 are conditional or architecture-specific necessities, not a
complete classification of all planar counterexamples.

## 9. Connection ledger

| source | target and map | preserved predicate | destroyed information | required sidecar | cheapest decisive test |
|---|---|---|---|---|---|
| quantum hopping `H` | resistor generator, `H_xy -> |H_xy|^2` by strong dephasing | leading population semigroup | loop phase, circulation | magnetic cycle holonomies | triangle `K^3` sign check |
| differential `DF` | normalized intensity table `C/T` | asymptotic rank-one proximity | determinant Wilson phase | dark/channel-degenerate regime | THM-3548 two matrix hostiles |
| JC response matrix `A_P` | bipartite conductance graph `|A_wu|^2` | support, magnitudes, matching, leaves, polygon bounds | cancellation, rank on cycles, target range, factorization | gain phases plus Segre cycle binomials | `A_+` versus `A_-` square hostile |
| triangular row `(j,t+1)` | odd sum-row via `M` | bijection, additive row, norm up to factor two | primitive content without parity sheet | index-two/Gaussian content | `M^TM=2I`, axis completion |
| sum-row vector | scalar square norm | radius and two-adic lower data | orientation and determinant sign | ordered exponent / `P^1(F_2)` color | `(1,2),(2,1)` against `(1,4)` |
| odd-row maximum `2C_t` | filler `Q_t=2C_t+1` | collision-free row sentinel | support-point semantics | diagonal fixed point if support matters | compare `(t+1,t+1)` with `Q_t` |
| higher-dimensional collision | categorical quotient | invariant functions, collision image | transverse character | ramification divisor | THM-3543 square factor |
| ramified torus quotient | polynomial corrections `(A,B)` | polynomiality and chosen collision if imposed | low transverse repair cells | both outputs, `w`-degree `(4,5)`, common-power top | THM-3549 boundary solve |
| higher-dimensional collision | coordinate-hypersurface restriction | tangent directions and ambient unit determinant | displayed polynomial graphs are all obstructed | nonlinear coordinates or a nongraph hypersurface | THM-3553 top row, then THM-3546 four gates |
| transverse boundary collision | Catalan formal thickening, then universal cubic root cover after adjoining `r` | selected collision and the exact ramification-cancellation mechanism | polynomiality in `w`; constant Jacobian after polynomialization | mixed terms that move `R=0` and send branching to infinity | fixed-line first jet, then mixed response boxes |
| mixed Catalan polynomial ansatz | coefficient-row recurrence `E_k` | polynomiality and selected collision | low width cannot terminate the tail | width/degree state plus global degree passport | width 3 degree 6 as hostile; viable cap at least 105 |
| curved ambient collision surface | punctured Kummer cover `(s,b)->(b,4s^2)` | finite etale collision | affine-plane completeness | boundary valuation / nonconstant unit | try mixed function-field deformation |
| cusp packet `(L,T,U,S)` | marked-root/Kummer inverse cubic, then two nonlinear polynomial outputs | escaping root, odd discriminant owner, full immersion | decomposability and descent can fail in projection | Pluecker relation, integrability, order index, wandering prime | solve the descending minor-Bezout system |
| leading form `H^a,H^b` | first Euclidean coefficient tower | exact divisibility | later renewal/new faces | root multiplicities and full Newton packet | bound `(6.11)` then next layer |

## 10. Next exact experiments

1. **Root-divisor continuation.**  Insert `D_k(H)` into the next homogeneous
   Keller equation and determine whether its quotient has a new exact
   valuation law or whether a hostile cross term first breaks closure.
2. **Gain-core census.**  Enumerate small Newton supports, build the complete
   equal-sum response graph, peel singleton zero rows, minimize cycle rank
   over exact target shears, and solve the remaining gain holonomies.
3. **Parallelogram tiling.**  Start with `(7.3)`, include every cross leak,
   and search for the smallest closed multi-fibre support rather than a
   hand-picked cancellation square.
4. **Torus-correction boundary.**  Start at transverse degrees `(4,5)` with
   top coefficients `(c h^4,d h^5)`, impose the inherited two-point collision,
   and solve the remaining Jacobian rows; prioritize the clean `(100,125)`
   cell by writing `P_99=u^25w^2C_72`, only `73` homogeneous coefficients,
   and reject any retained contracted line.
5. **Mixed Catalan thickenings.**  Use width-three degree six only as the
   cheapest renewal-mechanism probe; it is globally too small.  For an actual
   counterexample search, jump to sparse coefficient cap at least `105`,
   impose the height-108 pair `(72,108)`, and solve the recurrence with support
   pruning.  In the universal cubic chart, require the deformation to move
   `R=p+3t^2=0` at order zero; reject every ansatz that merely fixes the line
   and changes its normal jet.
6. **Kummer boundary deformations.**  Add the smallest mixed transverse terms
   that change `s^2=delta/4` before filling the divisor; test whether the
   missing divisor moves to infinity or a finite ramification valuation
   survives.
7. **Cusp nonlinear projections.**  Search low invariant degree for
   `A(L,T,U,S),B(L,T,U,S)` using exact response linearization in `B`, then
   impose descent, the Pluecker relation, and differential integrability
   before testing the coordinate-order index or periodic image prime.  The
   current finite box already rejects natural and constant-linear
   projections.
8. **Coordinate-hypersurface descent.**  Stop searching displayed graphs,
   which THM-3553 excludes in every degree.  Instead enumerate low-degree
   nonlinear ambient coordinate changes and nongraph coordinate
   hypersurfaces, then test exact divisibility and coordinate status, using
   the plane `x=0` as the positive hostile control.
9. **Finite-dephasing diagnostics.**  Use the first nonzero loop correction,
   not the resistor limit, as a numerical detector for the coefficient
   holonomies most likely to control cancellation.

## 11. Scope boundary

- The strong-dephasing theorem is a finite-dimensional singular limit, not a
  reduction of JC to electrical networks.
- THM-3548 supplies necessary intensity, polygon, Segre, shear, and conditional
  escape filters.  None is sufficient, and the Puiseux bound assumes the
  displayed branch expansion.
- The arithmetic row maps are exact, but norms and residues are not invariant
  under `GL_2(C)`.  They can schedule an integral support search; they cannot
  impose a complex planar-JC obstruction without an integral normalization
  and coefficient sidecar.
- The `(6,4)` proof is an elementary internal certificate for a class already
  excluded by stronger cited results.
- The new divisor `(6.11)` concerns only the first post-subleading layer and
  only `b<a<2b`; it is not an infinite tower.
- The cited 2022 sub-`125` list is imported rather than reproved.  It excludes
  the weaker height-`105` screen and leaves `(72,108)` as the current first
  passport; it does not assert that such a counterexample exists.
- The immersive cusp packet has a proved marked-root/Kummer inverse-cubic
  factorization, but its nonlinear decomposable projection and graph-descent
  criterion are construction programs, not counterexamples.
- The Catalan thickening is a formal/local holomorphic map in `(v,w)`, not a
  polynomial endomorphism.  Its polynomialized `(v,r)` model is a ramified
  finite cubic cover, not Keller; THM-3555 excludes only corrections fixing
  the whole ramification line, while mixed surgeries moving it remain open.
- THM-3557 closes mixed Catalan widths one and two in all coefficient degrees
  and width three only through cap five.  Its cap-six cell is internally
  open but globally excluded by the cited degree gates; no mixed thickening
  is constructed.
- The punctured Kummer map is a genuine finite etale collision, but its source
  is `G_m x A^1`; THM-3554 excludes only affine-plane fillings with the same
  quadratic function field.
- THM-3553 excludes polynomial graphs only in the displayed coordinates of
  one fixed ambient cubic map.  It does not exclude graphs after nonlinear
  ambient changes or nongraph coordinate hypersurfaces.
- THM-3549 closes one-sided, low-total-degree, and transverse-affine repairs
  of one fixed quotient seed; the mixed `(4,5)`-and-above correction box is
  still open.  Its `(100,125)` transfer assumes the transverse top rows also
  dominate ordinary degree.  THM-3550 is independently hostile-audited.
- `JC(2)` remains open.

The most productive conceptual compression is:

```text
equal-sum exponent fibres
    -> complex Jacobian gain graph
    -> resistor shadow under dephasing
    + holonomy sidecar for cancellation;

triangular square rows
    -> Gaussian--Hadamard sum rows
    -> exact parity/index loss;

leading common form
    -> Euclidean divisibility
    -> multiplicity-sensitive codimension;

collision near-counterexamples
    -> ramified polynomial quotient / nonpolynomial Catalan thickening /
       ramified connected cubic root cover / punctured etale Kummer cover
    -> move finite branching to nonproper escape through a mixed-boundary
       deformation;

cusp-square packet
    -> marked root escaping at `v=0` + quadratic Kummer pair
    -> descending integrable decomposable projection problem.
```

Those are genuine maps with explicit preserved and destroyed data.  They
leave a narrower counterexample search without pretending that an analogy is
a proof.
