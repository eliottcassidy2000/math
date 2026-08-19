# Planar Jacobian counterexamples through dephasing, conductance, and hidden phase

**Status: PROVED synthesis with new canon THM-3544, THM-3545,
THM-3548, THM-3549, and THM-3550, integrating incoming THM-3554;
VERIFIED-EXACT quotient and Catalan experiments; OPEN counterexample
constructions.**  No polynomial planar
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
| Anchor | hypothetical planar Keller pair | target pencil and Newton faces | every member composite, degree at least six; height at least eight | PROVED |
| Anchor | THM-1300 quotient | `(v,w)->(R,S)` | collision plus conductor `-2w^2` | PROVED + EXACT |
| Niche | dephased quantum walk | edge amplitudes and conductances | phases disappear at leading order, return as cycle flux | PROVED finite-dimensional expansion |
| Niche | Jacobian coefficients | fibre polygons inside a Segre matrix | local closure plus global even-cycle binomials | PROVED |
| Wildcard | self-intersecting boundary curve | formal transverse thickening | exact constant Jacobian, Catalan nontermination | PROVED + EXACT |
| Wildcard | invariant graph in the 3D map | coordinate hypersurface | normal/tangent unit factorization | PROVED criterion, OPEN search |
| Wildcard | curved THM-1300 collision surface | punctured Kummer cover | finite etale collision, nonconstant source unit | PROVED + EXACT |

The board changed twice.  First, the quotient, formal boundary, and curved
collision surface became a three-way pincer: polynomial-but-ramified,
Keller-but-nonpolynomial, and etale-but-punctured.  Second, the dephasing analogy
stopped being an attempted proof bridge and became a phase-aware search
architecture.

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

This is the first decisive bridge lesson.  Kirchhoff theory is obtained for
free only after the information carrying determinant interference has been
discarded.  Finite-dephasing spectroscopy, rather than the infinite limit,
is the natural analogue of a Jacobian phase sidecar.

## 4. Five typed connection contracts

| source | target and map | predicate preserved | information destroyed | required sidecar | cheapest decisive test |
|---|---|---|---|---|---|
| quantum hopping | resistor edge, `h_ij -> |h_ij|^2` | graph and leading jump rate | all gauge holonomy | Wilson phases on a cycle basis | compare a tree, triangle, and `K_2,2` |
| Jacobian coefficient pair | channel energy, `z_ab -> |z_ab|^2` | channel support and cancellation capacity | complex closing phase | fibre polygon phases | polygon inequality and exact coefficient sum |
| coefficient channels | independent fibre networks | each local zero sum | common factorization `p_aq_b` | Segre even-cycle binomials | rank-one minors before solving fibres |
| Keller differential | intensity matrix `DF -> (|DF_ij|^2)` | entry growth and near-rank-one profile | determinant plaquette phase | `arg(ad overline(bc))` | equal-intensity determinant-zero/one pair |
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
The open question is whether polynomial `v`-dependence in higher rows can
make the ladder terminate while the terminal equations and collision remain
valid.  This is a finite, degree-by-degree construction problem and has a
positive exact first row plus a hostile separated control.

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
| sparse fibres | avoid the classified `f(x)+g(x)z^d` cells | THM-3418 |
| inherited quotient repair | both outputs changed; transverse degrees at least `(4,5)`; row-dominant first live box `(100,125)` | THM-3549 + cited degree list |
| inherited finite collision | exact quadratic Kummer filling of the puncture is impossible | THM-3554 |
| displayed graph descent | impossible in every polynomial degree; nonlinear coordinate hypersurfaces remain | THM-3546/3553 |

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
before coefficient elimination.

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

### Normal multiplier as the missing resistor terminal — OPEN

In THM-3546 the normal derivative is the second factor whose product with the
tangential Jacobian is the ambient unit.  The torus quotient deletes this
terminal and obtains `w^2`.  A graph/hypersurface search should therefore
record the normal multiplier as an extra boundary node rather than eliminate
it.  The exact test is divisibility `(38)`, not similarity of Laplacians.

## 12. What changed, and what did not

The session produced five durable shifts.

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
5. The most credible searches are no longer generic coefficient boxes.  They
   are coherent-cycle constructions in the `(72,108)` Newton cages, the
   `(100,125)` row-dominant quotient repair (or a genuine second scale), and
   cubic branch-curve surgery which moves ramification while exporting it to
   infinity.

What did not change is the theorem status of the central problem:

```text
JC(2) remains OPEN.                                      (43)
```

The conductance picture supplies new exact filters and candidate generators,
not a proof by analogy.  Its real value is that it identifies, with equations,
which phases and global compatibility conditions a counterexample must keep.
