# Planar Jacobian counterexample search through dephasing gain graphs and square-row completion

**Research synthesis, 2026-08-18.**  Status labels are local to each claim.
`JC(2)` remains **OPEN**.  No counterexample is claimed here.  The session
produced one repaired canonical proof, a new root-multiplicity-sensitive
constraint on a hypothetical counterexample, an exact low-degree jet
closure, several refuted construction architectures, and two concrete open
counterexample programs.

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
4. **PROVED strong-dephasing limit.**  A Hermitian hopping matrix with
   `|H_xy|^2=c_xy`, under uniform site dephasing of coherence-decay rate
   `lambda`, has population generator

   ```text
   (2/lambda)L_c + O(lambda^-2).
   ```

   Thus `p(lambda*tau/2) -> exp(tau L_c)p(0)`.  Kirchhoff theory is exactly
   the leading slow limit.  Phases first return through magnetic loop
   corrections: triangles contribute at the next order.
5. **PROVED exact bridge, with a decisive loss ledger.**  The coefficient
   equation `Jac(P,Q)=1`, for fixed `P`, is a complex bipartite gain graph
   whose rows are grouped by equal exponent-vector sums.  Strong dephasing
   keeps squared edge magnitudes but discards the cycle holonomies that
   control cancellation, rank, and range.  A resistor graph is therefore a
   useful pruning shadow, not a lossless JC carrier.
6. **PROVED arithmetic merge.**  The user's triangular square rows and odd
   equal-sum rows are related by the Gaussian--Hadamard map

   ```text
   M(x,y)=(y-x,y+x),  det M=-2,  M^T M=2I.
   ```

   Hence odd-row norms are exactly twice the corresponding triangular-row
   norms.  The filler `3,11,27,51,...` is one more than the transformed row
   maximum.  It is a collision-free scalar sentinel, while the diagonal
   exponent is the functorial support filler.
7. **FINITE-EXACT positive seed; OPEN projection.**  An explicit
   four-coordinate cusp-square packet `A^2 -> A^4` is everywhere immersive,
   although none of its six natural two-coordinate projections, nor any
   constant-linear projection, is Keller.  Nonlinear polynomial projection
   is the strongest self-contained counterexample search seed found here.

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
- Least-used sidecar: the normalization/order-index distinction in the fixed
  higher-dimensional Keller packets.  It is not a planar theorem, but it is
  the correct hostile test for proposed inverse-cover constructions.

### Five live concepts

1. leading-form Euclidean divisibility and root multiplicity;
2. equal-sum coefficient fibres and magnetic holonomy;
3. Gaussian--Hadamard row transport and its index-two parity sidecar;
4. cusp-square discriminant packets and nonlinear projection;
5. coordinate-hypersurface descent from an ambient Keller collision.

The session changed the board as follows.  The row transport became an exact
map, but failed the ordinary Keller predicate by a Laurent Jacobian factor.
The resistor analogy became exact at leading dephasing order, but the gain
holonomy became the mandatory sidecar.  The cusp construction lost every
obvious two-coordinate projection, but gained an everywhere-immersive
four-coordinate packet.  The Euclidean lane gained an explicit
root-multiplicity divisor rather than another untyped recurrence guess.

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
of `H`, hence `g>=2`.  The corrected Appelgate--Onishi theorem says a
counterexample cannot have either component degree equal to a product of at
most two primes; its proof history includes the later Nowicki--Nakai repairs
([original](https://doi.org/10.1016/0022-4049(85)90099-4),
[lemma repair](https://doi.org/10.1016/0022-4049(88)90069-2)).  Thus each
component degree must have at least three prime factors counted with
multiplicity.  Moh's degree computation shows that both degrees cannot be at
most 100 ([Moh's summary](https://www.math.purdue.edu/~ttm/jacobian.html)).

These are cited global restrictions, not results of this companion.

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
variables sends one polynomial coordinate graph into another, its restriction
is Keller in `n` variables; any collision on the graph descends.

For an ambient three-variable collision, the search becomes:

1. find a source coordinate polynomial `rho_s` vanishing on two colliding
   preimages;
2. find a target coordinate `rho_t` vanishing at their common image;
3. solve `rho_s | rho_t(F)` exactly;
4. verify both are polynomial coordinates, not merely smooth hypersurfaces.

Passing these gates would produce a planar counterexample with no second
Jacobian search.  Interpolation is cheap; divisibility and coordinate status
are load-bearing.

The categorical torus quotient is a proved hostile, not a solution:
[THM-3543](../01-canon/theorems/THM-3543-torus-quotient-ramification-square-no-go.md)
computes an exact quotient Jacobian `2(2-3v-t)^2`.  The quotient retains the
collision only by contracting a divisor and acquiring ramification.  Graph
restriction and quotient forgetting are opposite operations.

### D. Cusp-square packet and nonlinear `4 -> 2` projection: strongest explicit seed

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

This packages a cubic inverse discriminant as one genuine component times a
square.  Every natural two-coordinate output is obstructed:

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

### E. Other attractive architectures that close

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
3. both component degrees have at least three prime factors with
   multiplicity, and at least one degree exceeds 100;
4. the subleading correction `W` vanishes, by repaired THM-3025;
5. in `b<a<2b`, every root multiplicity obeys the extra divisor tax
   `(6.10)`--`(6.12)`;
6. the first Euclidean block must survive its layer-by-layer codimensions and
   then create a genuine new face rather than assume an endless scalar
   recurrence;
7. any sparse coefficient ansatz must survive leaf peeling and all cross
   leaks of its complete product support;
8. a magnitude-only or norm-only shadow is insufficient: determinant signs,
   coefficient phases, and cycle holonomies must be retained;
9. exact finite even symmetry is forbidden; a near-symmetric candidate needs
   a load-bearing asymmetric defect;
10. a cover/cusp construction must retain normalization, coordinate-order
    index, Jelonek/nonproperness data, and hidden discriminant components.

The last four are architectural necessities within the stated programs, not
a complete classification of all planar counterexamples.

## 9. Connection ledger

| source | target and map | preserved predicate | destroyed information | required sidecar | cheapest decisive test |
|---|---|---|---|---|---|
| quantum hopping `H` | resistor generator, `H_xy -> |H_xy|^2` by strong dephasing | leading population semigroup | loop phase, circulation | magnetic cycle holonomies | triangle `K^3` sign check |
| JC response matrix `A_P` | bipartite conductance graph `|A_wu|^2` | support, magnitudes, matching, leaves | cancellation, rank on cycles, target range | gain phases on a cycle basis | `A_+` versus `A_-` square hostile |
| triangular row `(j,t+1)` | odd sum-row via `M` | bijection, additive row, norm up to factor two | primitive content without parity sheet | index-two/Gaussian content | `M^TM=2I`, axis completion |
| sum-row vector | scalar square norm | radius and two-adic lower data | orientation and determinant sign | ordered exponent / `P^1(F_2)` color | `(1,2),(2,1)` against `(1,4)` |
| odd-row maximum `2C_t` | filler `Q_t=2C_t+1` | collision-free row sentinel | support-point semantics | diagonal fixed point if support matters | compare `(t+1,t+1)` with `Q_t` |
| higher-dimensional collision | categorical quotient | invariant functions, collision image | transverse character | ramification divisor | THM-3543 square factor |
| higher-dimensional collision | coordinate graph restriction | tangent directions, ambient unit determinant | nothing if coordinates verified | coordinate/divisibility certificate | THM-3546 four gates |
| cusp packet `(L,T,U,S)` | two polynomial outputs | cusp-square discriminant if chosen invariantly | immersion can be lost in projection | all six minors, order index, wandering prime | nonlinear degree-box response solve |
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
4. **Cusp nonlinear projections.**  Search low invariant degree for
   `A(L,T,U,S),B(L,T,U,S)` using exact response linearization in `B`, then
   reject by Jacobian factors, coordinate-order index, or periodic image
   prime.  The current finite box already rejects natural and
   constant-linear projections.
5. **Coordinate-graph descent.**  For an ambient Keller collision, solve the
   divisibility `rho_s | rho_t(F)` at increasing coordinate degree, with
   positive controls from triangular automorphisms and hostile smooth
   noncoordinate hypersurfaces.
6. **Finite-dephasing diagnostics.**  Use the first nonzero loop correction,
   not the resistor limit, as a numerical detector for the coefficient
   holonomies most likely to control cancellation.

## 11. Scope boundary

- The strong-dephasing theorem is a finite-dimensional singular limit, not a
  reduction of JC to electrical networks.
- The arithmetic row maps are exact, but norms and residues are not invariant
  under `GL_2(C)`.  They can schedule an integral support search; they cannot
  impose a complex planar-JC obstruction without an integral normalization
  and coefficient sidecar.
- The `(6,4)` proof is an elementary internal certificate for a class already
  excluded by stronger cited results.
- The new divisor `(6.11)` concerns only the first post-subleading layer and
  only `b<a<2b`; it is not an infinite tower.
- The immersive cusp packet and graph-descent criterion are construction
  programs, not counterexamples.
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
    -> multiplicity-sensitive codimension.
```

Those are genuine maps with explicit preserved and destroyed data.  They
leave a narrower counterexample search without pretending that an analogy is
a proof.
