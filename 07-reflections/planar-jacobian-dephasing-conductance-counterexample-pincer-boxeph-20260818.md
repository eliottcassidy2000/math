# Planar Jacobian counterexamples through dephasing, conductance, and hidden phase

**Status: PROVED synthesis with new canon THM-3544, THM-3545,
THM-3548--THM-3550, THM-3586, and THM-3587, integrating
THM-3551--THM-3557;
VERIFIED-EXACT quotient, dephasing, Catalan, and packet experiments; OPEN
counterexample constructions.**  No polynomial planar
counterexample and no proof of `JC(2)` is claimed.  The main outcome is a
smaller counterexample passport, a three-way near-counterexample pincer, and
three sharply typed construction cells.

The motivating quantum-to-classical statement is true with normalization
and graph-channel qualifications.  Its most useful contribution to the
Jacobian problem is also its warning: dephasing deletes exactly the cycle
phases that distinguish determinant one from determinant zero.  The right
JC search strategy is therefore not to replace amplitudes by conductances,
but to use conductances as a magnitude filter while retaining Wilson phases,
coefficient-polygon phases, and global rank-one sidecars.

## 1. Inheritance pass and live concept board

The closest proved planar mechanisms are the one-, two-, and three-degree
source-fibre closures THM-2063, THM-2071, and THM-2118; the all-positive-
weight proper-power gate THM-2102; and the exact target-response quotient
THM-2230.  The canonical hostile object is the explicit three-dimensional
Keller collision in THM-1300.  It cannot simply be projected: THM-3543 shows
that its categorical torus quotient keeps the collision only by acquiring a
square ramification divisor.  The corrected near miss is MISTAKE-421: an
etale full fibre does not make the derivative of a chosen coordinate
eliminant nonzero.  The least-used sidecars brought into this session are
Dirichlet energy, cycle holonomy, the normal multiplier of a hypersurface,
and the global Segre equations behind coefficient cancellation.

The live board was:

| lane | object | representation | invariant / obstruction | status |
|---|---|---|---|---|
| Anchor | hypothetical planar Keller pair | target pencil and Newton faces | internal `>=6/8`; cited `Omega>=3`, height `>=108` | PROVED + CITED |
| Anchor | THM-1300 quotient | `(v,w)->(R,S)` | collision plus conductor `-2w^2` | PROVED + EXACT |
| Niche | dephased quantum walk | edge amplitudes and conductances | phases disappear at leading order, return as cycle flux | PROVED finite-dimensional expansion |
| Niche | Jacobian coefficients | fibre polygons inside a Segre matrix | local closure plus global even-cycle binomials | PROVED |
| Niche | Keller derivative plaquette | consecutive fibres `T=bc`, `T+kappa=ad` | factor allocation, curls, and coefficient-span observer | PROVED + EXACT |
| Wildcard | self-intersecting boundary curve | formal transverse thickening | exact constant Jacobian, Catalan nontermination | PROVED + EXACT |
| Wildcard | invariant graph in the 3D map | coordinate hypersurface | displayed graphs fail; nonlinear coordinate charts remain | PROVED / OPEN boundary |
| Wildcard | curved THM-1300 collision surface | punctured Kummer cover | finite etale collision, nonconstant source unit | PROVED + EXACT |
| Wildcard | inverse-cubic owner packet | noninjective differential-rank-two `A^2->A^4` packet | normal multiplier misses the singular-minor ideal in all degrees | PROVED + VERIFIED-EXACT CLOSED |
| Wildcard | nodal cylinder | cap-38 width-`(4,6)` two-scale scaffold | all-line passport and period pass; residual Keller PDE | PROVED gates + FINITE-EXACT OPEN |

The board changed repeatedly.  First, the quotient, formal boundary, and curved
collision surface became a three-way pincer: polynomial-but-ramified,
Keller-but-nonpolynomial, and etale-but-punctured.  Then the dephasing analogy
became a phase-aware search architecture, classical bounds moved the first
degree cell to `(72,108)`, displayed graphs closed in all degrees, and the
differential-rank-two four-coordinate packet exposed decomposability as another
lost phase sidecar.  Its exact image ideal then closed that packet in every
target degree.  Finally, the nodal cylinder passed the missing normal-ideal
gate and moved the first surviving construction to a two-scale `(4,6)` width
cell whose collision, every-line, and period tests can all be paid exactly.
In parallel, retaining the coherent differential plaquette before taking
magnitudes turned the determinant into two consecutive reducible fibres and
closed every unbalanced four-factor cell, every positive toric carrier, and
every coefficient-span cell of dimension at most two.

## 2. The strong-dephasing theorem, with the normalization exposed

Let `G=(V,E)` be a finite simple undirected graph.  Give edge `{i,j}` a
positive conductance `c_ij` and choose a Hermitian hopping Hamiltonian

```text
H_ij=sqrt(c_ij) exp(i phi_ij),       H_ji=conj(H_ij).      (1)
```

Arbitrary real onsite energies are allowed.  Put
`Pi_i=|i><i|` and consider

```text
dot rho=-i[H,rho]
 + eta sum_i kappa_i(2Pi_i rho Pi_i-Pi_i rho-rho Pi_i),
kappa_i>0.                                               (2)
```

On the matrix unit `E_ij`, `i!=j`, the dephaser in `(2)` acts by

```text
-eta(kappa_i+kappa_j)E_ij.                               (3)
```

If the physical evolution is observed at time `eta tau`, the diagonal
populations converge uniformly on compact slow-time intervals to

```text
dot p_i=sum_j [2c_ij/(kappa_i+kappa_j)](p_j-p_i).         (4)
```

For uniform `kappa_i=1`, `(4)` is exactly the variable-speed resistor
generator

```text
dot p_i=sum_j c_ij(p_j-p_i).                             (5)
```

### Mechanism

Let `Pcal` project matrices onto their diagonal, `Qcal=1-Pcal`, and put
`K=-i[H,-]`.  Then `Pcal K Pcal=0`.  On the slow space the leading Schur
complement is

```text
-Pcal K (D_kappa|_Qcal)^(-1) K Pcal,                    (6)
```

where `D_kappa(E_ij)=-(kappa_i+kappa_j)E_ij`.  Direct
substitution in `(6)` gives the off-diagonal rate in `(4)`.  Only
`|H_ij|^2` survives.

For symmetric rates, the Dirichlet form is

```text
-<f,Lf>=sum_{ {i,j} } c_ij |f_i-f_j|^2.                 (7)
```

Kirchhoff harmonicity, hitting probabilities, and effective resistance now
follow from the usual finite Laplacian identities.  The time normalization
matters: for generator `(5)` on a connected `N`-vertex graph, the commute
time is `N R_eff`.  The familiar `2(sum_e c_e)R_eff` belongs to the
constant-speed discrete-time normalization.

### Exact boundary conditions

The slogan needs four qualifications.

1. Nonuniform dephasing produces rates `2c_ij/(kappa_i+kappa_j)`.  One global
   rescaling recovers all original conductances exactly only when the edge
   sums `kappa_i+kappa_j` are constant.  On a connected nonbipartite graph
   this forces uniform rates; on a connected bipartite graph the two parts
   may carry complementary constants.
2. If dephasing does not separate every pair of vertices, coherent blocks
   survive and the slow state need not be a nodewise Markov chain.
3. Parallel coherent edges add as amplitudes before squaring.  The original
   parallel conductances add only if the edge channels are orthogonal or are
   separately dephased.
4. Disconnected components retain their masses.  Resistance across
   components is infinite, and the dictionary is componentwise.

These are normalization and channel issues, not failures of the singular
limit.

## 3. The first phases to return are cycle fluxes

For uniform off-diagonal decay `Gamma`, write the physical generator as

```text
K-Gamma Qcal.
```

After the initial layer, its population reduction begins

```text
G_eff=Gamma^(-1) A_2+Gamma^(-2) A_3+O(Gamma^(-3)),       (8)
```

with

```text
(A_2 p)_i=2 sum_j |H_ij|^2(p_j-p_i),                    (9)
(A_3)_ij=-6 Im(H_ij(H^2)_ji)
         =-6 sum_k Im(H_ij H_jk H_ki),       i!=j.      (10)
```

The sign in `(10)` uses `K=-i[H,-]`; replacing `H` by `-H` reverses the odd
term.  On the leading-normalized slow scale, physical time
`s=(Gamma/2)tau`, the rate is

```text
q_ij=c_ij-(3/Gamma)sum_k Im(H_ij H_jk H_ki)
          +O(Gamma^(-2)).                                (11)
```

For `(2)`, `Gamma=2eta`.  Thus the first phase-sensitive coordinate is an
oriented triangle Wilson flux.  Triangle-free graphs have no correction at
this order; on a tree all phases are removable by vertex gauge at every
dephasing strength.  More generally, path counting says that a cycle of
length `g` first permits gauge-sensitive data at slow order
`Gamma^{-(g-2)}`.  In particular the four-edge plaquette needed for a
`2 x 2` determinant is invisible to the classical limit and to the triangle
correction, but appears at fourth commutator order.

That last order statement has an exact finite control.  On the unit cycle
`0-1-2-3-0`, put the entire Wilson variable `z=exp(i Phi)` on the edge
`3->0`.  With `B=Pcal K Qcal`, `C=Qcal K Pcal`, and
`D=Qcal K Qcal`, exact matrix-unit computation gives

```text
B D C=0,
B D^2 C(z)-B D^2 C(1)
  =(z+z^(-1)-2)
    [ 2 -4  6 -4
     -4  2 -4  6
      6 -4  2 -4
     -4  6 -4  2 ].                                  (11a)
```

Every row of the displayed matrix sums to zero, and
`z+z^(-1)-2=2(cos Phi-1)`.  Hence the `C_4` flux first appears in the raw
fourth Schur-walk coefficient and at order `Gamma^(-2)` after leading time
normalization, exactly as the girth rule predicts.  This is `FINITE-EXACT`:
the [script](../04-computation/dephasing_c4_flux_order_control.py) and
[stored output](../05-knowledge/results/dephasing_c4_flux_order_control.out)
claim the Wilson factor and its first order, not convention-dependent
phase-blind memory terms in a fully normalized effective generator.

The exact control extends to every girth.  With
`A_m=BD^(m-2)C`, every raw word below graph girth `g` is phase blind, while
the first Wilson monomials occur in `A_g`.  On the unit `g`-cycle, if `S` is
the cyclic shift and `R_g=(I-S)^g`, then

```text
Delta A_g=i(-1)^((g-1)/2)(z-z^-1)R_g,       g odd,
Delta A_g=  (-1)^(g/2)(z+z^-1-2)R_g,        g even.    (11b)
```

Thus odd cycles return as sine/skew/oriented signals and even cycles as
cosine/symmetric signals; overlapping girth cycles add their embedded copies.
The statement is for literal population coordinates and the raw Schur series
`B(Gamma-D)^(-1)C` (convergent for `Gamma>||D||`, or formal in
`Gamma^-1`).  Arbitrary phase-dependent coordinate changes can manufacture
earlier-looking terms, and nonconstant onsite potentials remove some
zero-onsite comparison symmetries without changing the girth threshold.

On the weighted `K_2,2` built from
`A=[a b;c d]`, put `W=ad conjugate(bc)` and `R=|abcd|`.  The same fourth word
gives

```text
Delta A_4=(W+conjugate(W)-2R)(I-S)^4,                 (11c)
|det A|^2=|ad|^2+|bc|^2-2R-(Delta A_4)_(00)/2.        (11d)
```

So the resistor shadow plus its first plaquette correction reconstructs the
determinant modulus exactly.  For a Keller matrix, the scalar coefficient is
`q=(|ad|-|bc|)^2-|kappa|^2 in [-|kappa|^2,0]`; after `kappa=1`, its operator
norm is at most `16`.  When `|ad||bc|` grows, the Wilson phase is forced toward
constructive alignment, while THM-3548's channel-degenerate branch remains a
separate unsqueezed regime.  The complete signs, graph-wide formula, sharp
boundaries, and hostile controls are in the
[cycle-flux note](dephasing-cycle-flux-girth-law-boxeph-20260820.md),
[companion](../04-computation/dephasing_cycle_flux_girth_control.py), and
[output](../05-knowledge/results/dephasing_cycle_flux_girth_control.out).

This is the first decisive bridge lesson.  Kirchhoff theory is obtained for
free only after the information carrying determinant interference has been
discarded.  Finite-dephasing spectroscopy, rather than the infinite limit,
is the natural analogue of a Jacobian phase sidecar.

## 4. Six typed connection contracts

| source | target and map | predicate preserved | information destroyed | required sidecar | cheapest decisive test |
|---|---|---|---|---|---|
| quantum hopping | resistor edge, `h_ij -> |h_ij|^2` | graph and leading jump rate | all gauge holonomy | Wilson phases on a cycle basis | compare a tree, triangle, and `K_2,2` |
| Jacobian coefficient pair | channel energy, `z_ab -> |z_ab|^2` | channel support and cancellation capacity | complex closing phase | fibre polygon phases | polygon inequality and exact coefficient sum |
| coefficient channels | independent fibre networks | each local zero sum | common factorization `p_aq_b` | Segre even-cycle binomials | rank-one minors before solving fibres |
| Keller differential | intensity matrix `DF -> (|DF_ij|^2)` | entry growth and near-rank-one profile | determinant plaquette phase | `arg(ad overline(bc))` | equal-intensity determinant-zero/one pair |
| Keller differential | matching carrier `[[a,b],[c,d]] -> (T=bc,T+kappa=ad)` | determinant and two factor partitions | edge allocation and both row curls | UFD allocation plus integrability | toric `SL_2` hostile and balanced-factor census |
| THM-1300 collision | torus quotient `(x,y,z)->(xy,x^2z)` | polynomiality and a collision | character/normal direction | a coordinate hypersurface or normal multiplier | test divisibility in THM-3546 |

The table explains why a bare Laplacian resemblance cannot imply `JC(2)`.
It also explains how the resemblance remains useful: every destroyed column
becomes an explicit variable that a counterexample search must retain.

## 5. Exact conductance-shadow gates for a Keller pair

[THM-3548](../01-canon/theorems/THM-3548-planar-keller-conductance-shadow-gates.md)
packages four new necessary filters.

### 5.1 Differential intensities become rank one at infinity

For a pointwise differential

```text
A=DF=[a b;c d],       C=(|A_ij|^2),       T=sum C_ij,
r=|ad|,               s=|bc|,
```

the fixed determinant gives

```text
|r-s|<=|kappa|<=r+s,
|det C|<=|kappa|T/2,
dist_F(C/T,{rank<=1})<=|kappa|/T.                       (12)
```

Hence every high-gradient counterexample branch has an asymptotically
rank-one intensity table.  There are only two regimes after subsequencing:

- **dark plaquette:** `r+s->infinity`, `r/s->1`, and the two large matching
  phases align so their difference remains `kappa`;
- **channel-degenerate shadow:** `r+s` stays bounded while one entry grows
  and its determinant partner collapses.

The same four unit intensities can have determinant zero or one, so `(12)`
is a necessary magnitude theorem, not a phase-free characterization.

### 5.2 Coefficient fibres are polygons coupled by a Segre variety

Write

```text
P=sum_a p_a X^a,        Q=sum_b q_b X^b,
z_ab=det(a,b)p_aq_b.                                    (13)
```

The coefficient of `X^(s-(1,1))` in `Jac(P,Q)` is

```text
sum_(a+b=s) z_ab.                                       (14)
```

Every nonconstant fibre therefore closes a complex polygon.  If it has `k`
nonzero channels, conductances `c_e=|z_e|^2`, and total energy `E`, then

```text
max sqrt(c_e)<=sum_rest sqrt(c_e),
c_max<=(k-1)E/k.                                        (15)
```

A singleton fibre is impossible; a two-channel fibre has equal magnitudes
and opposite phases.  But local polygons are not enough.  The matrix

```text
W_ab=p_aq_b                                             (16)
```

has rank one globally, so all its even-cycle alternating products agree.
This is the coefficient version of Wilson holonomy: independently choosing
a closing phase in each fibre usually destroys the possibility of one
global pair `(P,Q)`.

### 5.3 Legal target shears have global frustration

For `u=grad P`, `v=grad Q`, `alpha=||u||^2`, and
`beta=<v,u>/alpha`,

```text
||v+h u||^2=|kappa|^2/alpha+alpha|h+beta|^2.            (17)
```

The pointwise Schur minimizer is `h=-beta`.  A polynomial target shear can
use it only if `beta=B(P)` globally on all `P`-fibres.  If that happens, the
sheared gradients become Hermitian-orthogonal, their squared norms multiply
to a unit, and both gradients are constant.  The map is a target shear of an
affine automorphism.  A counterexample must therefore have

```text
beta notin C[P].                                        (18)
```

THM-2230 shows that the legal shears tested here exhaust all mates with the
same first component and constant response.  The obstruction is genuine
global shear frustration, not failure to guess the best pointwise scalar.

### 5.4 Puiseux contact cannot be arbitrarily flat

If a degree-`d` Keller map has a specified Laurent--Puiseux escape

```text
x(t)=b t^(-r)+...,
F(x(t))=a+q t^m+...,              r,m>0,                (19)
```

then

```text
||DF(x(t))||_F^2 >= const |t|^(-2(r+m)),
||DF(x(t))||_F^2 <= const |t|^(-2r(d-1)),
m<=r(d-2).                                              (20)
```

The first inequality is the exact two-dimensional identity
`||DF^{-1}||_F=||DF||_F/|kappa|` applied to the differentiated branch.  The
second is polynomial degree.  The normalized intensity matrix approaches
rank one at least as fast as `O(|t|^(2(r+m)))`.  This does not import the
curve-selection theorem needed to produce `(19)` from arbitrary
nonproperness; it is conditional on the displayed branch.

### 5.5 The coherent plaquette is a consecutive-factor pencil

The phase sidecar can be retained algebraically before taking magnitudes.
Write

```text
a=P_x,       b=P_y,       c=Q_x,       d=Q_y,
T=bc,                    T+kappa=ad.                   (20a)
```

The determinant gives the two consecutive fibres in `(20a)`, but their
factor allocation must also satisfy

```text
a_y=b_x,                 c_y=d_x.                      (20b)
```

The all-degree UFD/curl analysis in
[THM-3587](../01-canon/theorems/THM-3587-consecutive-keller-fibre-factor-toric-and-coefficient-span-gates.md)
and its
[`reflection`](keller-consecutive-factorization-integrability-plaquette-boxeph-20260820.md)
proves the following additional passport.  In a hypothetical nonautomorphism,
all four derivative entries are nonunits,

```text
Omega(T),Omega(T+kappa)>=2,
omega(T)+omega(T+kappa)>=4,                            (20c)
```

and equality in the distinct-factor bound is forced to be balanced `(2,2)`.
A carrier `T=h(x^p y^q)` is impossible in every degree.  More generally,
`T=h(R)` is tame whenever every distinguished level `R-r`, for a root of
`h(h+kappa)`, is irreducible; a scalar-carrier construction must spend a
reducible selected fibre.

There is an independent coefficient observer.  Let `S_F` be the span of the
nonconstant coefficient matrices of `DF`.  If `dim S_F<=2`, its Frobenius
annihilator contains a rank-one matrix and exposes a constant directional
derivative, so the map is tame.  At rank three the only surviving annihilator
is invertible and a constant source gauge gives

```text
F=(lambda x+H_y,lambda y-H_x),
T=-(H_xy^2+mu),          T+kappa=lambda^2-H_xy^2.      (20d)
```

In this trace gauge the distinct-factor floor strengthens from four to five.
Thus a minimal balanced `(2,2)` matching table must have full coefficient
span four; the first open rank-three cell must already have a reducible
distinguished `H_xy` level.

The sharp positive control is a three-shear automorphism with balanced
fibres, gradient-coprime `T`, and full span four.  Hence none of these fields
is sufficient.  A useful generator that pays determinant and factor balance
up front is

```text
E_+(V)E_-(U)E_+(W)
 =[1+UV, V+W(1+UV); U,1+UW],                           (20e)
```

leaving only two explicit curl PDEs, quotient gauges, collision, every-line,
and nonproperness.  The nodal boundary frame factors into the opposite
three-shear orientation plus one final right shear; its unipotent branch
holonomy is therefore exactly the coherent datum omitted by the matching
products.

## 6. The global degree passport is much smaller

[THM-3544](../01-canon/theorems/THM-3544-planar-keller-target-pencil-total-degree-six-floor.md)
and
[THM-3550](../01-canon/theorems/THM-3550-prime-degree-exclusion-and-pencil-height-eight-floor.md)
give the strongest **internal** numerical squeeze from this session.

For a nonautomorphic complex planar Keller pair:

```text
every nonzero target-pencil member has composite total degree >=6;
the target-pencil height is >=8;
every member has degree >=4 along every linear source-fibre direction. (21)
```

The degree-six floor first combines the linear/quadratic/cubic fibre
closures with a projective root of the top binary form.  Degree five needs a
second weight: ordinary weight forces `lambda L^5`, the fibre floor supplies
`cy^4`, and weight `(4,5)` exposes the power-free binomial
`lambda x^5+cy^4`.

The prime exclusion repeats that idea without fixing the prime.  If
`R_p=lambda x^p`, increase the positive weight ratio until the first lower
`y`-term ties `x^p`.  If this face were an `m`th power, its vertex `(p,0)`
would give `m|p`, hence `m=p`; its other vertex would then have positive
`y`-exponent divisible by `p` despite having total degree at most `p-1`.

After a target change, the pencil has a reduced basis

```text
n=deg R < m=deg S,                                     (22)
```

with exactly two degree values: `[R]` has degree `n`, and every other
direction has degree `m`.  Its leading forms have the typed common base

```text
R_n=cH^a,       S_m=dH^b,
n=a deg H,      m=b deg H,      a,b>=2,                (23)
```

where `H` is power-free.  Both `n,m` are composite, `n>=6`, and `m>=8`.

Three classical `CITED` gates make the viable box much smaller than this
elementary boundary.  Nagata's repaired Appelgate--Onishi theorem, together
with Magnus, says that a component whose degree has at most two prime
factors counted with multiplicity belongs to an automorphic pair.  Completing
each nonzero pencil member to a target basis therefore gives

```text
Omega(deg R)>=3                  for every pencil member. (21a)
```

Guccione--Guccione--Horruitiner--Valqui list all hypothetical
counterexample degree pairs of height below `125`, leaving only `(72,108)`
and its transpose.  Applied to `(22)`, this gives

```text
m>=108;       if m<125, then (n,m)=(72,108).            (21b)
```

Finally, Shastri's one-point-at-infinity criterion says the common leading
base has at least two distinct projective roots: a top form `lambda L^n`
would make the pair automorphic.  These are imported results, not new proofs
of this session.

Consequently the internal height-eight equality is a sharp hostile boundary,
not a viable counterexample cell.  If that internal floor were sharp, the
degree pattern would be `(6,8)` and only two leading architectures would
survive:

```text
(R_6,S_8)=(cL^6,dL^8),
(R_6,S_8)=(c(xy)^3,d(xy)^4).                            (24)
```

The cited degree theorem excludes both, and the cited infinity theorem
separately excludes the one-root row.  The two-root architecture remains a
useful negative control for the internal gates and also enters the repaired
subleading rigidity of THM-3025.  More precisely, in the mirrored
`g=2,(a,b)=(3,4)` normalization,
THM-3016 and repaired THM-3025 force the degree-seven/degree-five rows into
`Q_7=lambda H P_5`, or force both to vanish; a generic translation places
the empty branch in the same relation.  Those two theorems still await an
independent hostile audit after MISTAKE-422, so this is a secondary gate, not
part of the fully audited degree chain.

## 7. The polynomial/algebraic/punctured collision pincer

Three exact objects now approach a counterexample from different sides.

### 7.1 Polynomial and colliding, but ramified

In `w=2-3v-t` coordinates the THM-1300 torus quotient is

```text
R=w(4v+6-3(v+1)^2w),
S=w^2((v+1)(v+2)-(v+1)^3w),
Jac(R,S)=-2w^2.                                         (25)
```

It has generic degree three and the exact fibre

```text
(R,S)^(-1)(0,0)=V(w) union {(0,2)}.                     (26)
```

Thus it is a polynomial, noninjective plane map inherited from a genuine
Keller collision.  Its sole fatal defect is not small: the whole line `w=0`
is contracted and supplies the Jacobian square.  THM-3543 proves that every
polynomial target-invariant postcompression keeps this factor.

### 7.2 Keller and colliding, but algebraic rather than polynomial

On the self-intersecting boundary curve

```text
v -> (v^2,v^3-v),          (+1),(-1) -> (1,0),          (27)
```

take the separated transverse ansatz

```text
P=v^2+A(w),
Q=v^3-v+vB(w).                                          (28)
```

THM-3545 proves

```text
Jac(P,Q)=kappa
iff B=3A/2 and A-3A^2/4=kappa w.                        (29)
```

The unique normalized solution is

```text
A=(2/3)(1-sqrt(1-3kappa w)),                            (30)
```

whose Taylor coefficients are positive Catalan multiples and never
terminate.  This is an exact formal/holomorphic/algebraic local Keller map
with a transverse collision, but not a polynomial map on `A^2`.

Equations `(25)` and `(30)` are dual conductors.  Flattening `w^2` in `(25)`
requires a pole or fractional power; polynomiality fails on the quotient
side.  Inverting the quadratic conductor in `(29)` produces the square root
and infinite Catalan tail; polynomiality fails on the boundary-thickening
side.  A real counterexample must couple variables strongly enough to avoid
both one-variable inversion obstructions.

### 7.3 Etale and colliding, but punctured

[THM-3554](../01-canon/theorems/THM-3554-punctured-kummer-collision-surface-normal-form.md)
restricts the same three-dimensional map to its curved collision component

```text
C=2-3xy-x^2z=0.                                        (30a)
```

Because `x(3y+xz)=2` there, `x` is a unit and the surface is not `A^2` but

```text
G_m x A^1 = Spec C[s,s^(-1),v],       s=x^(-1).        (30b)
```

Explicit Laurent source and polynomial target automorphisms put the
restricted map in the exact normal form

```text
(s,b) -> (b,4s^2): G_m x A^1 -> A^1 x G_m.            (30c)
```

It is a finite etale Kummer cover of degree two and its deck involution
contains the known collision.  Filling `s=0` gives the affine map
`(s,b)->(b,4s^2)` with Jacobian `-8s`.  More strongly, an everywhere-etale
`A^2` filling with the same quadratic function-field extension is
impossible: the pullback of the reduced divisor `delta=0` would have odd
valuation by etaleness and even valuation because `delta=4s^2`.

The three failures are sharply complementary:

```text
quotient:     polynomial plane + collision, but a contracted divisor;
thickening:   local Keller + collision, but an infinite algebraic tail;
restriction:  finite etale + collision, but a punctured Laurent plane.
```

A genuine counterexample has to retain the finite collision, remove the
puncture without restoring a finite branch divisor, and terminate
polynomially.  Live escapes must change the Kummer function field through
mixed corrections, send the missing divisor to infinity by a nonproper
affine modification, or use higher-sheet noninjectivity carried entirely by
the asymptotic set.

## 8. The quotient repair box begins much later than expected

[THM-3549](../01-canon/theorems/THM-3549-torus-quotient-correction-no-go.md)
sets

```text
P=R+A(v,w),              Q=S+B(v,w)                    (31)
```

for the seed `(25)` and proves:

1. keeping either output is impossible in every correction degree;
2. arbitrary corrections of total degree at most two have empty
   constant-Jacobian coefficient variety;
3. through total degree three the constant is forced to zero, with or
   without preserving the collision `(0,2),(-3/2,0)`;
4. corrections affine in `w`, with arbitrary polynomial dependence on
   `v`, fail identically; and
5. a collision-preserving Keller correction must have sorted final
   `w`-degrees at least `(4,5)`.

At the first surviving boundary, the transverse tops are already forced to

```text
P=c h(u)^4 w^4+lower rows,
Q=d h(u)^5 w^5+lower rows,           u=v+1.             (32)
```

This follows from the top coefficient

```text
5a'b-4ab'=0.                                            (33)
```

The first live correction problem is therefore a resonant common-power
network whose lower rows must cancel the inherited `-2w^2`, not a generic
low-degree perturbation.

## 9. Ranked counterexample architectures

### A. The cited first degree cell: a `(72,108)` coherent-cycle ansatz

If the reduced height is below `125`, `(21b)` leaves one degree cell.  In a
maximal common-base normalization its ordinary leading forms are

```text
R_72=c K^2,              S_108=d K^3,       deg K=36,  (34)
```

and Shastri forces `K` to have at least two distinct projective roots.  A
cheap sparse hostile seed is `K=x^17 y^19`; a dense control is a product of
36 distinct linear forms.  Neither leading pair is Keller--the top gradients
share the expected base--but they expose the exact rows that lower terms must
repair.  For every root direction of `K`, attach a lower Newton face which is
itself a proper power, raises the linear-fibre degree to at least four, and
shares no global gradient factor.

Search those lower coefficients in this order:

1. reject singleton coefficient fibres and polygon-energy violations;
2. impose all Segre rank-one cycle binomials;
3. solve the exact Jacobian fibre sums;
4. impose a normalized two-point collision; and
5. only then run nonproperness/Puiseux tests.

The key creative change is to require at least one genuine three-or-more-
channel fibre cycle.  Two-channel cancellations are phase-rigid and tend to
collapse into binomial/common-factor geometry; a three-channel polygon is
the first place holonomy can move while magnitudes stay admissible.

The low-degree scaffold

```text
(x^3+alpha y^2)^2,       (x^4+beta y^2)^2
```

remains a sharp negative control for the new internal height-eight argument,
but the classical degree and one-point-at-infinity gates exclude it from the
live counterexample search.

### B. High-transverse repair of the inherited quotient

Use `(31)` with the exact collision constraints and top form `(32)`.  Split
the first search into:

- `h` constant, where the top transverse rows have no additional root
  divisor; and
- `h` nonconstant, where every root of `h` must be checked for fibre-degree
  drop and a critical gradient before solving lower rows.

The lower `w`-rows should be organized as coefficient polygons, not as a
flat Groebner box.  The fixed source term `-2w^2` is an external amplitude in
the degree-two fibre; all other nonconstant fibres must close to zero.  The
global Segre cycles then decide whether the separately chosen closures come
from one pair `A,B`.

There is also a concrete ordinary-degree threshold.  Put `N=deg h+1`.  If
the displayed top-`w` rows in `(32)` dominate total degree, then the pair has
degrees exactly

```text
(4N,5N).                                                (34a)
```

The cited sub-`125` list excludes every such pair with `N<25`; none has ratio
`72:108`.  The first row-dominant box not classically excluded is therefore

```text
deg h=24,                 (deg P,deg Q)=(100,125).      (34b)
```

Below that threshold, a candidate must create a lower-`w`, higher-`u` row
which supplies a genuine second Newton scale.  If the resulting height stays
below `125`, that second scale must land exactly in the `(72,108)` cell.

### C. Mixed finite completion of the Catalan thickening

Keep the boundary collision `(27)` but drop separation:

```text
P=sum_(i=0)^N p_i(v)w^i,
Q=sum_(j=0)^N q_j(v)w^j,
p_0=v^2,                 q_0=v^3-v.                    (35)
```

The exact coefficient ladder is

```text
sum_(i+j=r+1) [j p_i' q_j-i p_i q_j']
   =kappa delta_(r,0).                                  (36)
```

The first Bezout row has the simple solution

```text
p_1=1,                    q_1=3v/2                     (37)
```

for `kappa=1`.  Separation turns `(36)` into the infinite Catalan tail.
THM-3557 proves that mixed width one and two fail in every coefficient
degree, and that width three fails through coefficient degree five.  Its
first internally open cell is therefore `(N,D)=(3,6)`, where `D` bounds all
coefficient degrees.

That is not yet a globally admissible cell.  The ordinary degrees in `(35)`
are at most `max(3,D+N)`, while every solution retains the boundary
collision and would be a counterexample.  The cited degree floor therefore
forces

```text
D+N>=108.                                              (37a)
```

At width three one needs `D>=105`, and at width four `D>=104`; below height
`125`, a target-reduced basis must have degree values `(72,108)` (the
displayed components may both have degree `108`).  Thus low width remains
structurally attractive, but only with very high longitudinal coefficients
or a deliberately nonuniform row profile.

THM-3555 identifies what this ladder is trying to desingularize.  After
adjoining the Catalan square root, the separated map is affinely equivalent
to the universal marked-root cubic cover

```text
(t,p) -> (p,-t^3-pt).
```

It is a connected generic `S_3` cover, but ramifies on `p+3t^2=0`.  Any
polynomial correction which fixes that entire line pointwise still has zero
Jacobian at its cusp.  Thus a viable mixed construction must move the branch
curve at order zero, retain only selected collision fibres, and export every
remaining sheet merger to infinity.  This is a much narrower surgery problem
than arbitrary termination of `(36)`.

### D. Coordinate-hypersurface descent from the 3D counterexample

THM-3546 gives a literal route that the torus quotient does not test.  Find
source and target coordinate polynomials `rho_s,rho_t` such that

```text
rho_t(F)=rho_s A                                         (38)
```

and the source coordinate hypersurface contains two colliding preimages.
Block-triangular differentiation then factors the ambient Jacobian unit into
the tangential Jacobian and normal multiplier; both are constants.  A
collision descends to a planar Keller counterexample automatically.

Interpolation through the collision is cheap.  The load-bearing tests are
global divisibility `(38)` and proof that both `rho`'s are polynomial
coordinates.  The torus quotient fails because it discards the normal
character and contracts a divisor; an invariant graph must retain that
normal direction until the unit factorization is made.

THM-3554 closes the most obvious curved surface: `C=0` contains the collision
and its restriction is etale, but its coordinate ring has the nonconstant
unit `x`, so it is `G_m x A^1`, not a coordinate plane.  A concurrent
`FINITE-EXACT` calculation first found no invariant same-graph section of
degree at most three, but THM-3553's top-form argument closes the displayed
graph lane in every degree.  For

```text
f=(F_1(x,y,h),F_2(x,y,h)),             h in C[x,y],
```

if `D=deg h>=1`, the top forms are

```text
A=x^3y^3 h_D,             B=3x^3y^2 h_D,
A=(y/3)B,                 Jac(A,B)=-(1/3)B B_x !=0.   (38a)
```

The last nonvanishing is robust: writing `B=3x^3y^2h_D`, the operator
`3+x partial_x` has positive eigenvalues `3,...,D+3` on the degree-`D`
monomials.  This is the unique top Jacobian row, of degree `2D+9`.  If `h=c`
is constant, including `c=0`, the same identity holds with

```text
A=x^2y^3(cx+3y),          B=3x^2y^2(cx+3y).           (38b)
```

Thus `Jac(f)` is never constant.  THM-3546 would force it to be a nonzero
constant if the displayed source graph mapped scheme-theoretically into
*any target graph in the displayed target coordinates*.  No such graph pair
exists.  What remains open is a nongraph coordinate hypersurface, a graph
after nonlinear ambient coordinate changes, or a different three-dimensional
Keller map.  This scope is essential: the plane `x=0` does restrict to a
Keller map and maps into `F_3=0`.

### E. Sparse resonant additive layers are a negative control

The tempting two-layer family

```text
F(z)=z+H_d(z)+K_(2d-1)(z)                               (39)
```

with homogeneous vector layers does not provide the required complexity.  A
degree split of `det JF=1` gives

```text
tr JH=0,
tr JK+det JH=0,
cross(JH,JK)=0,
det JK=0.                                               (40)
```

The last equation makes the components of `K` proportional.  After linear
conjugacy, `K=(0,T)`; the cross equation is `Jac(H_1,T)=0`.  Coprimality of
`d` and `2d-1` forces common powers of one linear form, and substitution in
the first two rows collapses the family to

```text
(x,y)->(x,y+alpha x^d+beta x^(2d-1)).                   (41)
```

Thus a two-resonance additive shear network is tame.  This pushes a viable
counterexample toward three or more coupled layers/cycles, consistent with
the coefficient-polygon diagnosis.

### F. Nonlinear projection of the differential-rank-two inverse-cubic packet

THM-3556 supplies a different positive object.  For a polynomial `U(v,y)`,
put

```text
T=y^2-6vU,       S=y^3-9vUy,       L=v^2(8vU-y^2).
```

Then `S^2=T^3+27LU^2`, and `LX^3+TX+2U` splits into one marked root which
escapes at `v=0` and a quadratic Kummer pair.  For one explicit `U_*`, the
map

```text
Z=(L,T,U_*,S):A^2 -> A^4                              (41a)
```

has differential rank two everywhere: its six `2 x 2` minors `M_ij` generate
the unit ideal.  It is not injective.  Indeed the two distinct conjugate
points

```text
p_+/-=(8/27 +/- 4sqrt(5)/9, 3/2 +/- 3sqrt(5)/5)
```

have the common packet value

```text
Z(p_+)=Z(p_-)=(-6724/3645,57/20,27/40,27/40).        (41a')
```

The radical-free certificate is
`27v-20y+22=y^2-3y+9/20=0`.  Thus the packet is a
differential-rank-two polynomial surface parameterization with a genuine
geometric double fibre, not an algebraic embedding.  Every polynomial target
projection retains the collision.

The two tangent sheets are locally compatible with unit area: for

```text
B_0=(9627/304384)T-(205687/1369728)U_*,
```

one has `Jac(L,B_0)=1` at both points.  Hence the obstruction below is global,
not a first-jet contradiction at the self-intersection.  Nevertheless no
natural two-coordinate projection is Keller.

There is a dual visible cubic `Y^3-3TY+2S` with marked root `Y=y`; its
discriminant has the same odd Kummer owner `-L` as the escaping-root cubic.
A projection which hides the marked sheet must therefore retain a common
resolvent square class as well as the scalar cusp equation.

The would-be construction asks for nonlinear target polynomials `A,B` with

```text
Jac(A(Z),B(Z))
 =sum_(i<j)(A_iB_j-A_jB_i)(Z) M_ij=1.                 (41b)
```

This is a particularly literal dephasing bridge.  An arbitrary Bezout
combination of the six minors exists, just as the classical intensity shadow
may satisfy local conductance constraints.  A legal projection retains the
discarded coherent sidecar: its six coefficients must descend through
`C[L,T,U_*,S]`, lie on the Pluecker quadric

```text
c_12c_34-c_13c_24+c_14c_23=0,                         (41c)
```

and arise simultaneously from ambient polynomial potentials as
`dA wedge dB`.  These conditions must hold before pullback; relations only
modulo the image ideal need not lift.  Even a closed, de Rham-exact,
decomposable polynomial two-form need not equal `dA wedge dB` for polynomial
`A,B`, so “exactness” alone is not the missing integrability gate.
An exact hostile in four variables is

```text
omega=(dy-y dx) wedge (dz+z dx)=d(y dz+yz dx).
```

Its kernel contains `partial_t` and
`partial_x+y partial_y-z partial_z`; their common polynomial invariants are
`C[yz]`, so two polynomial potentials annihilated by the kernel have zero
wedge.  Thus this closed, exact, decomposable `omega` is not `dA wedge dB`.

Two finite exact relaxations first closed the bottom of this search.  For all target
polynomials of total degree at most two, the `139 x 91` bracket-coefficient
matrix has full/deleted-constant ranks `67,67`.  For total degree at most
three, all `34` nominal target functions and `561` pair brackets give a
`336 x 561` matrix with ranks `187,187`.  In each case the constant row is
already forced by the nonconstant rows even when bivector coefficients vary
independently, before Pluecker or potential-realizability constraints.
Therefore target degree at most three is impossible.

The packet image algebra first appeared as a finite shadow.  In integral
coordinates `(L,T,W,R)=(L,T,2U_*,2S)`, besides

```text
F_3=R^2-4T^3-27LW^2,
```

there is an exact `28`-term quartic relation `G_4`.  Pullback dimensions
through target degree seven agree exactly with the Hilbert function of a
`(3,4)` complete intersection: the kernel dimensions are

```text
0,0,1,6,20,50,104.                                   (41c')
```

At the packet double value, `grad G_4=(-21/5)grad F_3`, the singularity
expected from the two immersed sheets.  A subsequent exact elimination now
promotes this shadow: `ker Q[L,T,W,R]=(F_3,G_4)` is an absolutely prime
complete intersection.  Adjoining the marked root `y` gives its smooth
normalization, although the original plane parameterization is birational
and nonfinite because two vertical `W=0` lines are reached only from source
infinity.  The linear-subresultant coefficient `alpha` is simultaneously
the rational-inverse denominator, the main conductor owner, and the exact
normal multiplier

```text
det d(F_3,G_4,A,B)(Z)=9 alpha(Z) Jac(A(Z),B(Z)).
```

The generic unordered-pair conductor has a genus-two model; it contains the
known quadratic double fibre and has one explicitly recorded triple fibre.
The original
[relation probe](../04-computation/cusp_square_packet_image_relation_probe.py)
is retained as an independent low-degree control, while the
[full image companion](../04-computation/cusp_square_packet_image_ideal_thm3556.py)
and [output](../05-knowledge/results/cusp_square_packet_image_ideal_thm3556.out)
carry the exact proof gates.

The full image theorem turns the projection equation into an all-degree
obstruction.  If `N_ij=det(dF_3,dG_4,e_i,e_j)`, then

```text
N_ij(Z)=9 alpha(Z)M_ij                                  (41c'')
```

for all six oriented pairs.  Exact rational grevlex reduction gives

```text
alpha notin (F_3,G_4,N_ij),
Phi(alpha)=81/50,                                      (41c''')
```

where `Phi` extracts the `LR^2` coefficient of the reduced normal form and
annihilates the ideal.  Therefore even arbitrary descended coefficients
`C_ij` cannot satisfy `sum C_ij(Z)M_ij=1`: multiplying such an identity by
`9alpha(Z)` and using the full kernel would contradict `(41c''')`.  Actual
potential pairs are a strict subclass.  The
[all-degree companion](../04-computation/cusp_square_all_degree_descended_bivector_no_go_thm3556.py)
and [output](../05-knowledge/results/cusp_square_all_degree_descended_bivector_no_go_thm3556.out)
make all six Hodge signs, the thirteen-element basis, the nonzero remainder,
and three prime-field controls reproducible.

Before this closure, the exact collision activated the cited planar degree
floor.  Since the
source degrees of `(L,T,U_*,S)` are `(6,4,3,5)`, a target cap `E` gives source
height at most `6E`; the GGHV bound forces `E>=18`.  At the first cap, target
reduction must have source degrees `(72,108)`, and the top packet ledger is
forced to

```text
A_top~L_top^12=(L_top^6)^2,
B_top~L_top^18=(L_top^6)^3.                           (41c'''')
```

This cap-`18` `2:3` cusp cell remains an independent degree/face hostile, but
`(41c''')` proves that no packet cap survives.  The lesson is architectural:
a replacement collision surface must make its descended normal multiplier
belong to the image singular-minor ideal before Pluecker or potential
integrability is even worth testing.

### F'. A nodal cylinder passes the normal-ideal gate

The normalization

```text
(u,t)->(u,X=t^2-1,Y=t(t^2-1)),
H=Y^2-X^2(X+1)                                        (41N1)
```

is a differential immersion and identifies `t=+1` with `t=-1`.  Its Hodge
multiplier is `X`, but now

```text
X=-((3X+1)/2)H_X-(9/4)YH_Y               modulo H,    (41N2)
```

so the descended normal-ideal gate passes.  An explicit closed ambient
two-form also contracts to one; a diagonal-invariant argument proves that
particular certificate is not `dA wedge dB`, leaving the correct global
integrability sidecar exposed.

Classical degree and face gates put the first target cap at `38`, with reduced
source tops

```text
[t^35(t+lambda u)]^2,       [t^35(t+lambda u)]^3.     (41N3)
```

The inhomogeneous nodal relation is load-bearing: it reduces the literal cube
of `Y^12+lambda uXY^11` from ambient cap `39` to cap `38`.  But the resulting
width-`(2,3)` scaffold is dead.  Gwozdziewicz's injectivity-on-one-line
theorem demands noninjectivity on every affine source line, while his quoted
Abhyankar--Lang Newton theorem gives, for the reduced `(72,108)` pair,

```text
N_B=(3/2)N_A,                 width_u(B)=(3/2)width_u(A). (41N4)
```

The highest-`u` Jacobian row also forces
`a_m=alpha h^(m/d), b_n=beta h^(n/d)`, `d=gcd(m,n)`.  A root of `h` kills
width `(2,3)` to an immersed degree-at-most-two line map, hence to an
injective line.  Thus `(2,3)` is impossible and `(41N4)` makes `(4,6)` the
first surviving width.  At target cap `38`, its shared lower base satisfies
`deg_t h<=32`, genuinely below the degree-`35` infinity base.

There is an exact positive scaffold at that boundary.  Put

```text
h=XY^10,   C=Y^12+lambda uXY^11,   D=C+mu u^2h.       (41N5)
```

Then `D^2` has cap/source degree/width `(26,72,4)`, while the nodal relation
reduces `D^3` to `(38,108,6)` without changing its pullback.  With
`eta=3X+1`, the low rows

```text
A=D^2+u^2-1+eta Y/2,
B_0=D^3+u^3-u+3eta uY/4                              (41N6)
```

give the same immersed nodal boundary with unit jet on all three degree-drop
lines `t=-1,0,1` and retain the four-point conductor fibre.  On every other
affine line the restriction degrees remain in ratios
`(72,108),(70,105),(68,102)`, or `(4,6)`; Abhyankar--Moh divisibility then
rules out an injective restriction.  Thus `(41N6)` pays the full every-line
necessary gate, not merely a sampled set.

The remaining global sidecar is visible but not fatal.  For
`A=a+Yc,B=b+Yd`, the nodal action

```text
Phi_F(u)=integral_(-1)^1 A B_t dt,
Jac(A,B)=1  ==>  Phi_F=2u+constant.                   (41N7)
```

has an exact beta-integral coefficient formula.  No pair with the normalized
nodal boundary and only first-conductor-layer terms can satisfy the Keller
PDE; a genuine `(X,Y)^2` correction is forced.  Conversely, a three-moment
dual built from `X^(j+2)(X+1)`, `j=0,1,2`, repairs `(41N7)` for `(41N6)` by
target-cap-at-most-ten corrections, preserving widths, all root-line jets,
the four-point fibre, and every-line degree passport.  Its full Jacobian is
still nonconstant.  This is therefore a sharply typed counterexample seed,
not a counterexample: degree, Newton width, collision, normal ideal, all
affine lines, first jets, the first conductor layer, and the global period
have been separated from the residual Keller PDE.  Exact formulas are in the
[nodal-cylinder scaffold](nodal-cylinder-keller-projection-scaffold-boxeph-20260820.md)
and its width/period companion.

### G. Two nonparallel invariants with cancelling response classes

THM-3551 closes three all-degree one-invariant families, and THM-3552 closes
a broad two-channel Kummer family even when its first component is a
polynomial submersion of high Newton area and high generic genus.  The
obstruction is not degree: on the normalized generic fibre the forced mate
differential is nonzero holomorphic and hence not exact.

The first honest escape is therefore

```text
T=x^a y^b,       S=x^c y^d,       ad-bc=+/-1,
P=x Phi(T,S)+Psi(T,S),                                (41d)
```

with genuinely nonlinear dependence on both `T` and `S`.  The unimodular
exponent determinant minimizes new toric ramification, but is not enough by
itself: an exact hostile which is only affine-linear in `S` still carries a
nonzero holomorphic Newton-adjoint class and has no rational mate.  A positive
signal is instead two or more response differentials of the second kind,
with zero residues and opposite cohomology classes, whose cancellation is
structural before any coefficient cap is raised.

This is the fibre-cohomology version of the dephasing lesson.  One channel
cannot close; two are rigidly opposite; three are the first number that can
carry variable phase while retaining zero total response.  The cheap search
works in the two-dimensional charge lattice, computes pole divisors and
cohomology classes first, and solves for a polynomial mate only after a
class cancellation is visible.

### H. Alternating-factor defect invoices

A second speculative ansatz starts with several smooth nonparallel factors
`f_i` and common-power leading forms

```text
H=product_i f_i^(e_i),          P_0=H^2,       Q_0=H^3. (41e)
```

Rather than asking one approximate root to absorb every Jacobian defect,
choose the multiplicities so consecutive correction rows invoice different
factors.  The desired finite-termination signal is a periodic owner sequence
whose obstruction operator is nilpotent.  The first two-factor cusp probe
`f=y^2-x^3`, `g=y^2-2x^3`, `H=f^2g^3` failed at the next defect (rank `15`
versus augmented rank `16`), so it is a hostile, not evidence.  The live
experiment uses at least three factors or a different multiplicity vector
and records owner, residual divisibility, solution dimension and target-shear
freedom at every row.  This lane is `FINITE-EXACT exploratory`, not canon.

## 10. The consolidated counterexample passport

A hypothetical planar counterexample now has to satisfy all of the following
simultaneously.

| coordinate | necessary profile | source |
|---|---|---|
| total target degrees | internally `6<=n<m`, composite, `m>=8`; cited `Omega>=3`, `m>=108`, and sub-`125` pair only `(72,108)` | THM-3544/3550; Nagata; GGHV |
| leading forms | `R_n=cH^a,S_m=dH^b`, `a,b>=2`, with at least two projective roots | top Jacobian + THM-2102; Shastri |
| linear source fibres | every nonzero pencil member has degree at least four in every direction | THM-2063/2071/2118 |
| positive-weight faces | every face of every pencil member is a proper power | THM-2102/2740 |
| coefficient fibres | polygon closure and concentration bound | THM-3548 |
| global coefficients | all Segre/even-cycle binomials | THM-3548 |
| target response | nonremovable shear frustration | THM-2230/3548 |
| escape differential | dark or channel-degenerate rank-one intensity shadow | THM-3548 |
| Puiseux contact | `m_contact<=r(d-2)` for a branch in the displayed form | THM-3548 |
| derivative matching carrier | two reducible consecutive fibres; distinct-factor total at least four, balanced if equality; no toric ray; coefficient span at least three | THM-3587 |
| rank-three derivative span | in constant-trace gauge the matching distinct-factor total is at least five | THM-3587 |
| sparse fibres | avoid the classified `f(x)+g(x)z^d` cells | THM-3418 |
| inherited quotient repair | both outputs changed; transverse degrees at least `(4,5)`; row-dominant first live box `(100,125)` | THM-3549 + cited degree list |
| inherited finite collision | exact quadratic Kummer filling of the puncture is impossible | THM-3554 |
| displayed graph descent | impossible in every polynomial degree; nonlinear coordinate hypersurfaces remain | THM-3546/3553 |
| mixed Catalan width | `D+N>=108`; width three requires coefficient cap at least `105` | THM-3557 + cited degree list |
| cusp-square packet projection | closed in every target degree: its descended normal multiplier is outside the image singular-minor ideal | THM-3556 |
| nodal-cylinder projection | cap `38`; widths at least `(4,6)` with lower base degree `<=32`; genuine conductor-square term; every-line, holonomy, and period gates | THM-3586 + Gwozdziewicz/Lang |
| response channels | avoid all one-invariant rays and the two-channel cyclic Kummer cell; cancel fibre cohomology | THM-3551/3552 |

At the cited exceptional height `108`, the Puiseux inequality becomes
`m_contact<=106r`.  This is not sharp enough to prove properness, but it ties
the actual first degree cell to a quantitative infinity profile.

## 11. Underexplored operations and cheap probes

Several ideas remain worthwhile precisely because the connection contracts
are now explicit.

### Cycle-rank filtration of coefficient supports — OPEN

Treat the nonzero `W_ab` as edges of a bipartite graph.  On a tree, rank-one
phases can be gauged from vertex phases and there is no independent holonomy,
mirroring the exact gauge triviality of quantum hopping on a tree.  This
suggests the testable hypothesis:

> A finite Keller support whose active coefficient graph has no cycle, or
> whose nonconstant sum fibres all have at most two channels, is forced into
> a tame/binomial class.

This is not proved.  The cheap test is an exhaustive support enumeration at
small degree, applying singleton, polygon, Segre, and fibre-degree gates
before coefficient elimination.  THM-3551/3552 supply all-degree evidence
for the one-/two-channel boundary, but do not prove the general hypothesis.

### Cancellation slack as a search score — HEURISTIC

For a coefficient fibre define

```text
delta_s=sum_rest sqrt(c_e)-sqrt(c_max).                 (42)
```

Negative slack rejects a support/magnitude assignment.  Zero slack forces
all smaller amplitudes to align and is phase-rigid.  Positive slack measures
the room available for a closing polygon.  This is coordinate-dependent and
not a Keller invariant, but it is a useful branch-ordering score for exact
search.

### Finite-dephasing tomography — PROVED mechanism, OPEN use

The leading classical generator reveals edge magnitudes; the next orders
reveal triangle and longer-cycle Wilson data.  Applied metaphorically to
coefficient networks, this suggests reconstructing a candidate in layers:

```text
magnitudes -> short-cycle phases -> full Segre holonomy -> integrability.
```

The destroyed data are restored in the same order that commutator walks
discover cycles.  No physical experiment implies a Jacobian theorem, but the
ordering supplies a principled exact-solver schedule.

### Electrical resistance on the support graph — HEURISTIC

Effective resistance can rank which coefficient-cycle phases are most
globally constrained: an edge of high leverage belongs to few alternate
routes and should tolerate less phase adjustment.  The proposed observable
is leverage `c_e R_eff(e)` on the coefficient bipartite graph.  The cheapest
probe is to compare it with Groebner infeasibility across small support
atlases.  No coordinate invariance or proof implication is claimed.

### Normal multiplier as the missing resistor terminal — PROVED gate, OPEN generator

In THM-3546 the normal derivative is the second factor whose product with the
tangential Jacobian is the ambient unit.  The torus quotient deletes this
terminal and obtains `w^2`.  A graph/hypersurface search should therefore
record the normal multiplier as an extra boundary node rather than eliminate
it.  More generally, for a local-complete-intersection image with Hodge
identities `N_I(Z)=mu(Z)M_I`, a descended unit-minor certificate forces

```text
mu in (image relations, normal minors N_I).            (42a)
```

The cusp-square packet fails `(42a)` by the exact quotient dual
`Phi(mu)!=0`; the nodal cylinder passes it by the explicit Bezout identity
`(41N2)`.  This is a reusable first gate, not a sufficient criterion:
closedness, Pluecker decomposability, and global polynomial potentials remain
separate.  For a hypersurface `(42a)` specializes to the normal divisibility
test `(38)`, not similarity of Laplacians.

## 12. What changed, and what did not

The session produced seven durable shifts.

1. The motivating dephasing statement is exact for a simple graph with
   uniform node dephasing, with explicit corrections for nonuniform rates,
   multiedges, disconnected graphs, and time normalization.
2. The classical resistor shadow is provably insufficient for the Jacobian
   condition; the missing datum is a Wilson phase already visible in a
   `2 x 2` plaquette.
3. The new internal proof gives composite pencil degrees at least six and
   height at least eight; classical cited results sharpen the live search to
   `Omega>=3`, height at least `108`, and the unique sub-`125` cell
   `(72,108)`.
4. The fixed 3D collision yields a three-way quotient/thickening/puncture
   pincer, a sharply delayed quotient-repair cell, and an all-degree no-go for
   displayed polynomial graph descent.
5. The inverse-cubic packet is now a differential-rank-two *noninjective*
   polynomial parameterization with a prime complete-intersection image, but
   its conductor multiplier fails the singular-minor ideal test.  This closes
   every nonlinear target degree, not only the former cap-three census.
6. Before magnitudes are taken, the differential plaquette is a pencil of two
   consecutive reducible fibres.  UFD allocation plus the two curls close the
   unbalanced and toric cells, while coefficient observers route the last
   rank-three boundary to a five-factor Hessian gauge.
7. The most credible searches are no longer generic coefficient boxes.  They
   are coherent-cycle constructions in the `(72,108)` Newton cages, the
   `(100,125)` row-dominant quotient repair (or a genuine second scale), cubic
   branch-curve surgery at `D+N>=108`, and the cap-`38` nodal-cylinder cell
   whose multiplier gate passes and whose marked boundary already has a
   four-point fibre.  Nonparallel
   invariant channels and alternating factor owners are the two higher-risk
   mechanisms most clearly outside the proved one-/two-channel walls.

What did not change is the theorem status of the central problem:

```text
JC(2) remains OPEN.                                      (43)
```

The conductance picture supplies new exact filters and candidate generators,
not a proof by analogy.  Its real value is that it identifies, with equations,
which phases and global compatibility conditions a counterexample must keep.
