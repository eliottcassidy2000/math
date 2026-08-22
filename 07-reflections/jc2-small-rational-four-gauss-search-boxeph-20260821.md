# Small rational data in the planar problem: Wright's elementary-Jacobian boundary and the non-elementary frontier

**Status: CITED + INDEPENDENTLY REPROVED for Wright's all-length
elementary-Jacobian criterion; PROVED all-degree for the displayed
three-through-six-factor reductions; FINITE-EXACT small-rational four-, five-,
and six-word hostile atlases; OPEN only outside the elementary subgroup and
outside the stated finite atlases.**  No
counterexample to the planar Jacobian conjecture and no proof of `JC(2)` is
claimed.  Wright's 1978 weak Jacobian theorem is the decisive endpoint: a
planar Keller map is invertible exactly when its normalized Jacobian matrix is
a finite product of elementary matrices.  The compact reduced-word proof
developed below independently rederives the direction needed here.  Every
polynomially integrable three-element Gauss word is an explicit
composition of triangular automorphisms, and the first curl collapses every
four-element word to an explicit tame boundary.  A second `x`-degree split
closes every five-element word, while a constant-factor Bruhat reduction and
a lexicographic leading-term split close every six-element word.  The exact
scouts are therefore reproducible replays and positive-control censuses inside
proved strata.  The four-word scout solves both curls over complete rational
coefficient spaces and finds no
`delta!=0` nonunit-eligible first-curl survivor in its precisely stated
`20,400`-row universe.  At length five, two staged linear correction
equations emerge, and the finite signed-monomial base is empty in every
correction degree.  These are controls inside the proved all-degree theorem,
not evidence for a remaining length-five cell.  The first six-word atlas then
eliminates `10,000` signed-monomial base quadruples against the complete
degree-six middle-correction space, again as a finite control.  There is no
live elementary depth: the first honest generator is a certified
non-elementary `SL_2/E_2` core, such as Cohn's matrix.

The exact companion is
[`jc2_small_rational_four_gauss_search.py`](../04-computation/jc2_small_rational_four_gauss_search.py),
with frozen
[`output`](../05-knowledge/results/jc2_small_rational_four_gauss_search.out).
The first six-word companion is
[`jc2_small_rational_six_gauss_search.py`](../04-computation/jc2_small_rational_six_gauss_search.py),
with frozen
[`output`](../05-knowledge/results/jc2_small_rational_six_gauss_search.out).

## 1. Why “small rationals” should mean a small structural language

The verified dimension-three counterexample in
[`THM-1300-jacobian-counterexample-dixmier-A3-explicit.md`](../01-canon/theorems/THM-1300-jacobian-counterexample-dixmier-A3-explicit.md)
has small integral coefficients and a collision at small rational points.  It
is tempting to enumerate planar maps with similarly small expanded
coefficients.  That is the wrong first box: the cited planar degree gates are
already enormous, while a sparse composition can acquire large expanded
degree and denominators very quickly.

The search therefore enumerates a **derivative word**, not a dense pair of
expanded coordinate polynomials.  This implements the connection contract

```text
small rational coefficient idea -> derivative matrix class
map: coefficients become shear parameters around a certified core
preserved: exact determinant, rational height, elementary/non-elementary class
lost: row integrability and global properness
sidecar: both row curls, factor passport, coefficient span, collision test
decisive test: reject E_2 first; solve curls over QQ only in a non-elementary class
```

The closest proved mechanism is
[`THM-3587-consecutive-keller-fibre-factor-toric-and-coefficient-span-gates.md`](../01-canon/theorems/THM-3587-consecutive-keller-fibre-factor-toric-and-coefficient-span-gates.md):
a smallest derivative-table survivor may have balanced matching fibres
`(2,2)` and coefficient span four.  Its displayed three-factor Gauss chart was
left open.  The canonical hostile is still the determinant-one polynomial
matrix whose curls fail; determinant tables are never counted as candidates.
The all-length endpoint is classical and is cited to
[Wright, *The amalgamated free product structure of GL_2... and the weak
Jacobian theorem for two variables*](https://www.sciencedirect.com/science/article/pii/0022404987900041).

The final concept board was:

| lane | object | live predicate | cheapest decisive test | outcome |
|---|---|---|---|---|
| Anchor | four-element word | both curls, balanced factors, span four | `x`-degree of `UV` plus exact replay | entire chart tame |
| Niche | three-element word | polynomial integrability | `x`-degree of the lower curl | all-degree tame theorem |
| Wildcard | involution-free cell | escape known equivariant closures | five exact signed-permutation tests | all `168` controls evade the tests but are tame |
| Next box | Cohn core with elementary decorations | retain nonzero `SL_2/E_2` class while repairing both curls | transport cokernel and exact class certificate | open |
| Boundary | constant `V=v,Z=-1/v` | defeat the naive six-word leading term | exact Gauss cancellation | full word collapses to an affine-triangular map |
| Hostile | `U=W=x,V=0` | determinant without integrability | first curl | defect `-1` |

## 2. The all-degree three-Gauss theorem

Let `k` be a characteristic-zero field and put

```text
E_+(f)=[1 f;0 1],                    E_-(f)=[1 0;f 1].
```

For arbitrary `U,V,W in k[x,y]`, the three-element word is

```text
M_3=E_+(V)E_-(U)E_+(W)
   =[1+UV, V+W(1+UV); U, 1+UW].                       (1)
```

### Theorem

If both rows of `(1)` are closed polynomial one-forms, so that `M_3=DF`
for a polynomial map `F`, then `F` is a polynomial automorphism.  More
precisely, unless `U=0`, there are univariate polynomials `A,Phi,Psi` with

```text
W=A'(y),
s=x+A(y),
U=Phi'(s),
q=y+Phi(s),
V=Psi'(q),
F=(s+Psi(q),q)                                             (2)
```

up to additive target constants.  Its inverse is

```text
s=P-Psi(Q),        y=Q-Phi(s),        x=s-A(y).           (3)
```

Thus three elementary factors cannot support a counterexample; the later
theorems raise this lower bound to six.

### 2.1 The lower curl forces `W_x=0`

The lower row of `(1)` has curl equation

```text
U_y=(UW)_x.                                               (4)
```

Suppose first that `U` is nonzero.  If `W=0` or `deg_x W=0`, conclusion
`(6)` already holds.  Otherwise write

```text
m=deg_x U,          n=deg_x W,
U=u_m(y)x^m+...,    W=w_n(y)x^n+....
```

If `n>=2`, the left side of `(4)` has `x`-degree at most `m`, whereas
`(UW)_x` has degree `m+n-1>m` and nonzero leading coefficient
`(m+n)u_mw_n`.  If `n=1`, comparison of the `x^m` coefficients gives

```text
u_m'=(m+1)w_1u_m.                                        (5)
```

For nonzero polynomials this is impossible: the right side has degree at
least `deg u_m`, while the derivative has degree `deg u_m-1`; if `u_m` is
constant, `(5)` is zero equals a nonzero polynomial.  Therefore

```text
W=W(y).                                                   (6)
```

Choose `A'=W` and set `s=x+A(y)`.  Equation `(4)` becomes

```text
(partial_y-W partial_x)U=0.
```

In the polynomial coordinates `(s,y)`, this derivation is just
`partial_y`, so its kernel is `k[s]`.  Hence `U=H(s)`; choose `Phi'=H` and
put `q=y+Phi(s)`.  The lower row of `(1)` is exactly `dq`.

### 2.2 The upper curl forces `V=V(q)`

After `(4)` and `(6)`, the upper curl reduces without division to

```text
U V_y-(1+UW)V_x=0.                                      (7)
```

Let

```text
L=U partial_y-(1+UW)partial_x.
```

In the polynomial coordinates `(s,q)`, direct calculation gives

```text
L(s)=-1,                 L(q)=0.                         (8)
```

Thus `L=-partial_s`, and `(7)` says `V=K(q)`.  Choosing `Psi'=K` gives
the first row `d(s+Psi(q))`.  Formulas `(2)--(3)` follow.

### 2.3 Boundary audit

- If `U=0`, then `M_3=[1,V+W;0,1]`; the remaining curl says
  `(V+W)_x=0`, and the map is a single triangular shear.
- Constants and zero values of `V` or `W` cause no division and are already
  included in `(2)`.
- Both coordinate changes used above are polynomial shears.  No generic-fibre,
  irreducibility, or analytic assumption is hidden.
- Characteristic zero is load-bearing only for antiderivatives and the
  nonvanishing leading coefficient in the degree comparison.

This proves the whole three-element chart, not merely a rational box.

## 3. The four-element equations

Now put

```text
M_4=E_-(Z)E_+(V)E_-(U)E_+(W).
```

Write

```text
A=1+UV,                 B=V+W(1+UV),
delta=U_y-(UW)_x.                                           (9)
```

Then

```text
M_4=[A,B; ZA+U,ZB+1+UW].                                  (10)
```

Its determinant is identically one.  Expanding the two row curls and using
the first in the second gives exactly

```text
U V_y-(1+UW)V_x+delta V = W_x,                            (11)
A Z_y-B Z_x = -delta.                                     (12)
```

For fixed `U,W`, both `(11)` and then `(12)` are linear systems in the
coefficients of `V` and `Z`.  This is the key search order: determinant,
first curl, second curl, factor/span passport, and only then collision.

The defect `delta` is the symmetry-breaking coordinate inherited from the
three-factor theorem.  A `delta=0` solution is already a three-shear prefix;
it is therefore the first boundary to remove.  Section 3.3 then closes the
apparently harder `delta!=0` branch as well.

### 3.1 All-degree closure of the `delta=0` four-word stratum

In fact, the entire `delta=0` stratum is tame, without degree or coefficient
bounds.  When `delta=0` and both curls hold, the right three-factor prefix has
both rows closed.  By the theorem of Section 2 it is the Jacobian of an
automorphism

```text
F_0=(p,q),                dp=A dx+B dy.
```

Equation `(12)` becomes

```text
J(p,Z)=A Z_y-B Z_x=0.
```

In the polynomial coordinates `(p,q)`, where `J(p,q)=1`, this is simply
`partial Z/partial q=0`.  Hence `Z=K(p)`.  Choosing `Omega'=K`, the full
four-word map is

```text
(p,q+Omega(p)),
```

a final triangular target shear after `F_0`.  Therefore any counterexample
in this four-word orientation must satisfy

```text
delta != 0.                                                (12a)
```

### 3.2 Curvature repair is a Hamiltonian-cokernel problem

After solving the first curl, write `dp=A dx+B dy` and let

```text
eta=U dx+(1+UW)dy
```

be the lower row of the three-factor prefix.  Its curvature coefficient is
`-delta`.  Rotate `eta` and `dp` into the polynomial vector fields

```text
D=-(1+UW) partial_x+U partial_y,
X_p=-B partial_x+A partial_y.
```

The first-row identity and the determinant identity give, respectively,

```text
div(X_p)=0,             D(p)=-(1+UW)A+UB=-1,             (12b)
div(D)=delta.                                               (12c)
```

Thus `(A,B)=grad p` is unimodular: `(1+UW)A-UB=1`.  This is a
Bezout certificate that `p` has no affine critical point; it does **not** by
itself say that `p` is a coordinate.  If the closed final lower row is `dq`,
then the last elementary factor replaces `eta` by `eta+Z dp`, or equivalently

```text
X_q=D+Z X_p.
```

Consequently

```text
div(X_q)=delta+X_p(Z)=delta+J(p,Z),
```

and the second curl is exactly the divergence-killing equation

```text
J(p,Z)=-delta.                                             (12d)
```

Thus the true obstruction is the class of `delta` in the algebraic
Hamiltonian cokernel

```text
k[x,y] / Im(J(p,·)).                                      (12e)
```

The map `Z -> J(p,Z)` is differentiation along the Hamiltonian trajectories
inside the fibres of `p`.  Solving for `Z` repairs the divergence locally,
but primitives can fail globally around a fibre; periods or residues are the
natural restoration sidecar.  Section 3.3 shows that the four-word equations
always resolve into a tame boundary; for the staged operators at length five
and beyond, the cheapest test remains a cokernel/period obstruction before
increasing a correction degree.

This has the same operation type as the nodal conductor-action obstruction:
a local correction can kill curvature only when its global fibre-period
class vanishes.  This is a typed mechanism analogy, not an identification of
the nodal and Gauss objects.

### 3.3 All-degree closure of the entire four-word chart

The first curl actually closes every remaining value of `delta`; no bound on
any parameter is needed.

**Theorem.**  Let `k` have characteristic zero.  For arbitrary
`U,V,W,Z in k[x,y]`, if both rows of `M_4` are closed, the integrated map is
a polynomial automorphism.  More precisely, it is an explicit composition
of affine maps and triangular shears.

Put

```text
R=UV,                         K=1+R.
```

The first curl has the scalar divergence form

```text
R_y=V_x+(KW)_x.                                             (12f)
```

**The `U=0` boundary.**  Equation `(12f)` gives

```text
V+W=C(y).
```

Thus `A=1,B=C(y)` and `p=x+Lambda(y)`, where `Lambda'=C`.  Here
`delta=0`, so `(12)` is `J(p,Z)=0`; in the coordinates `(p,y)` this says
`Z=S(p)`.  The full map is the pair of shears

```text
(p,y+Psi(p)),                Psi'=S.                       (12g)
```

If `U!=0`, `V=0` is possible only when `W_x=0`, and is included in the
univariate branch below.

**The `W_x!=0` branch.**  Assume `U,V!=0`.  Let

```text
ell=deg_x W>=1,             D=deg_x R,
R=r_D(y)x^D+... .
```

If `D>=1` and `ell>=2`, the term `(KW)_x` in `(12f)` has the unique top
`x`-degree `ell+D-1>D`, an immediate contradiction.  If `D>=1` and
`ell=1`, comparison of the `x^D` coefficients gives

```text
r_D'=(D+1)w_ell r_D,                                      (12h)
```

where `w_ell(y)` is the nonzero leading `x`-coefficient of `W`.  No nonzero
polynomial satisfies `(12h)`: its derivative has smaller degree, while the
right side has at least the original degree.  Hence `D=0`.

Because `k[y][x]` is a domain, `D=0` forces `U,V in k[y]`.  Equation `(12f)`
becomes

```text
K'=K W_x.                                                (12i)
```

If `ell>=2`, the right side of `(12i)` has positive `x`-degree unless
`K=0`.  If `ell=1`, `(12i)` is the same impossible polynomial logarithmic-
derivative equation unless `K=0`.  Therefore

```text
UV=-1.
```

The two polynomials `U,V in k[y]` are consequently constants `u,v` with
`uv=-1`.  Now `A=0,B=v`, and the second curl gives

```text
vZ_x+uW_x=0,              Z=u^2W+L(y).                   (12j)
```

The `W` terms cancel from the final lower row.  Up to target constants the
map is

```text
(p,q)=(vy,ux+G(y)),        -uv=1,                         (12k)
```

which is triangular after an affine coordinate exchange.

**The `W_x=0` branch.**  Write `W=W(y)`, choose `A_0'=W`, and set

```text
s=x+A_0(y).
```

In the polynomial coordinates `(s,y)`, equation `(11)` becomes

```text
V_s=(UV)_y.                                               (12l)
```

If `V=0`, then `p=s`; equation `(12)` is `Z_y=-U_y`, with both derivatives
taken at fixed `s`.  Hence `Z=-U+C(s)`, the full lower row is
`dy+C(s)ds`, and the map is `(s,y+Psi(s))`.

Suppose `V!=0`.  The first closed row has a polynomial potential `H`:

```text
dp=(1+UV)ds+Vdy=ds+dH,
H_s=UV,                    H_y=V.                        (12m)
```

Let `n=deg_y H>=1`.  The equation `H_s=U H_y` first forces
`deg_y U<=1`.  Write

```text
U=u_1(s)y+u_0(s),          H=h_n(s)y^n+... .
```

The coefficient of `y^n` in `(12m)` is

```text
h_n'=n u_1 h_n.                                          (12n)
```

The same polynomial degree comparison forces `u_1=0`.  Thus `U=U(s)` and
`delta=0`; Section 2 makes the three-factor prefix an automorphism, and
Section 3.1 makes the final factor a target shear.  This includes `W=0`,
constant `W`, `U=0`, and every zero-parameter boundary, completing the
proof.

The opposite alternating word is not a new cell.  With
`S=[[0,1],[1,0]]`, one has `S E_+(f) S=E_-(f)`, and the affine change
`F(z) -> S F(Sz)` preserves polynomial invertibility.  Adjacent elementary
factors of the same sign combine.  Therefore:

```text
every integrable elementary Gauss word of length at most four is tame.    (12o)
```

A counterexample expressed in this derivative language needs elementary
length at least five at this stage; Section 3.7 raises the bound to six.

### 3.4 The five-word correction system

The structural next box is therefore

```text
M_5=E_+(R_5)E_-(Z)E_+(V)E_-(U)E_+(W).                  (12p)
```

Let the right three-word prefix have rows

```text
(A,B),                    (U,1+UW),
epsilon=A_y-B_x,          delta=U_y-(UW)_x.
```

After the fourth factor, write the lower row as `(C,D)`:

```text
C=ZA+U,                   D=ZB+1+UW.
```

Its closure is exactly the first linear correction equation

```text
A Z_y-B Z_x+epsilon Z=-delta.                           (12q)
```

Once `(12q)` holds, left multiplication by `E_+(R_5)` leaves `(C,D)`
fixed and changes the top row to `(A+R_5C,B+R_5D)`.  Its closure is exactly
the second linear correction equation

```text
C (R_5)_y-D (R_5)_x=-epsilon.                           (12r)
```

Thus a bounded five-word search should enumerate only `U,V,W`, solve the
whole chosen coefficient space for `Z` in `(12q)`, and only then solve the
whole chosen space for `R_5` in `(12r)`.  Determinant, factor passport, and
collision tests come afterward.

There is another all-degree rejection inside this system.  If `epsilon=0`
and `(12q)` is solvable, the four-word prefix has both rows closed and is an
automorphism `(p,q)` by Section 3.3.  Equation `(12r)` becomes
`J(q,R_5)=0`; in the coordinates `(p,q)` it forces `R_5` to depend only on
`q`.  The last factor is the target shear `(p+Phi(q),q)`.  Consequently the
still-unclosed branch at this intermediate stage must satisfy

```text
epsilon != 0.                                             (12s)
```

Section 3.7 closes this branch too.

Equations `(12q)--(12r)` are the precise two-stage residual-correction
analogue of the Hamiltonian discussion in Section 3.2: the fourth factor
repairs the lower curvature, and the fifth repairs the surviving upper
curvature.  Fibre cokernels, periods, or residues should be computed at each
stage before increasing either correction degree.

### 3.5 First all-degree closures inside the five-word chart

**Affine-rightmost theorem.**  If `W` is affine-linear, every integrable
five-word `(12p)` is tame, with no degree bounds on `U,V,Z,R_5`.

If `W_x=0`, the rightmost factor is the Jacobian of a source shear; canceling
it leaves a four-word map.  Suppose `W_x!=0`.  Zero values of `U` or `V`
merge adjacent factors and again shorten the word.  For `U,V!=0`, let
`m=deg U`, `n=deg V`, and let `Z_N` be the top homogeneous part of a nonzero
`Z`.  Suppose first that `V` is nonconstant.  For every nonzero `Z`, the
unique top-degree part of `(12q)` is

```text
-partial_x(U_m V_n W_1 Z_N)=0.                           (12t)
```

The `x`-bearing linear form `W_1` cannot divide a pure power of `y`, so this
is impossible.  The boundary `Z=0` would require `delta=0`, which Section
2.1 excludes for `W_x!=0`.  Hence `V=v` is constant.

If `U` is nonconstant and `Z` is nonconstant, the same top identity `(12t)`
again gives a contradiction.  Thus a hypothetical solution in this branch
has `Z=z` constant, and equation `(12q)` reduces to

```text
(zv+1)delta-zW_x=0.                                     (12u)
```

If `U` is nonconstant, the top part of `delta` forces `zv+1=0`, after which
`zW_x=0` is impossible.  Thus `U=u` is constant as well.  Put `K=1+uv`.
At this constant-`U,V` boundary `Z` is again arbitrary; the cancellation
`K=0` is precisely why it could not be included in the preceding top-degree
argument.  For `K!=0`, the shift `H=Z+u/K` satisfies

```text
K H_y-vH_x-K(WH)_x=0.                                  (12v)
```

An `x`-degree comparison forces `H=0`: for `deg_x H>=1` it gives the same
logarithmic-derivative contradiction as `(12h)`, and the constant-in-`x`
case does too.  Thus `Z=-u/K`, the closed lower row is `(0,1/K)`, and final
closure forces `R_5=-K^2W+S(y)`, an affine-triangular map.  For `K=0`, one
instead gets `Z=u^2W+L(y)` and `epsilon=0`, so Section 3.4 supplies the final
target shear.  This proves the theorem, including every constant boundary.

Inside that tame stratum, the smallest apparent coefficient cell satisfies
the sharper statement that it dies before the last factor.

**Theorem.**  Suppose `U,V,W` are nonzero homogeneous linear forms.  Then
equation `(12q)` has no polynomial solution `Z`, of any degree.  Hence no
integrable five-word `(12p)` has such a base triple.  If one of `U,V,W`
vanishes, adjacent factors merge and Section 3.3 applies.

Let `Z_N` be the top total-homogeneous part of a hypothetical nonzero `Z`.
The unique top-degree part of the left side of `(12q)` is

```text
-partial_x(U V W Z_N)=0.                                 (12w)
```

Indeed, its degree is `N+2`, while `-delta` has degree at most one.  Equation
`(12w)` says that the homogeneous product `UVWZ_N` is a pure power of `y`.
Unique factorization forces

```text
U=uy,                     V=vy,                 W=wy     (12x)
```

with `u,v,w` nonzero.  In this remaining case

```text
A=1+uvy^2,
B=(v+w)y+uvwy^3,
epsilon=A_y,
delta=u.                                                   (12y)
```

Expand `Z=sum_(j=0)^N z_j(y)x^j` by `x`-degree.  If `N>=1`, the coefficient
of `x^N` in `(12q)` is

```text
(A z_N)'=0,
```

which is impossible because `A` is nonconstant and `z_N` is a nonzero
polynomial.  The `N=0` boundary is different and must be checked separately:

```text
(A z_0)'=-u,              A z_0=-uy+C.                   (12z)
```

The quadratic `A=1+uvy^2` cannot divide the nonzero affine numerator.  The
case `Z=0` also fails because `delta=u!=0`.  This exhausts the cell.

Thus the homogeneous-linear five-word scheme cannot be rescued by merely
raising the degree allowed for `Z`.

There is a complementary all-degree **monomial-height strengthening**.
Suppose `U,V,W` are arbitrary nonconstant single monomials, with arbitrary
positive exponents and nonzero scalar coefficients.  The same top-degree
identity `(12w)` holds.  If any of the three monomials has positive
`x`-exponent, the product `UVWZ_N` cannot be a pure `y`-power, so no nonzero
`Z` exists.  Nor can `Z=0`: the two possible terms in
`delta=U_y-(UW)_x` have different total degrees, and they cannot both vanish.

If all three monomials are pure `y`-powers, the `x`-degree argument above
again removes every `x`-dependent `Z`.  For `Z=z_0(y)`, equation `(12q)` is

```text
((1+UV)z_0)'=-U_y,
(1+UV)z_0=-U+C.                                         (12aa)
```

But `deg(1+UV)=deg U+deg V>deg(C-U)`, so divisibility is impossible.  Thus
no triple of nonconstant single monomials admits a polynomial first
correction, at any exponent or correction degree.

Allowing a single monomial to have exponent zero can make the first
correction solvable, but still cannot make the five-word map noninvertible.
If `W in k[y]`, its rightmost factor is the Jacobian of a source shear and
cancels, leaving a four-word map.  Suppose therefore that the monomial `W`
has positive `x`-exponent.  The same top-degree comparison eliminates every
case except `V=v` and `Z=z` constant, or the separate case in which both
`U=u,V=v` are constant.

In the first case, lower-row closure reduces exactly to

```text
(zv+1)delta-zW_x=0.                                    (12ab)
```

If the monomial `U` is nonconstant, its top product with `W` forces
`zv+1=0`, after which `(12ab)` contradicts `zW_x!=0`.  Thus `U=u` is also
constant.  At this separate constant-`U,V` boundary, `Z` is arbitrary.  Put
`K=1+uv`.  When `K!=0`, set `H=Z+u/K`; the first correction equation is

```text
K H_y-v H_x-K(WH)_x=0.                                 (12ac)
```

An `x`-degree comparison, identical to `(12h)`, forces `H=0`.  Hence
`Z=-u/K`, the closed lower row is `(0,1/K)`, and final closure forces

```text
R_5=-K^2W+S(y).
```

The resulting map is affine-triangular.  When `K=0`, instead
`Z=u^2W+L(y)` and `epsilon=0`, so Section 3.4 makes the fifth factor a target
shear after a four-word automorphism.  Zero monomial factors simply shorten
the word.  Therefore every integrable five-word whose base `U,V,W` consists
of single monomials, constants included, is tame.  A live base must contain
a genuinely multi-term parameter and cannot merely be a homogeneous-linear
triple.

### 3.6 A separated-variable nonlinear five-word cell is empty

The first natural way to add lower-order `x`-dependence also closes in all
correction degrees.  Let

```text
U=alpha(x+a(y)),          V=v(y),
W=c(y)x+b(y),             alpha in k*,
```

where `c` is nonzero and `v` is nonconstant.  Then `(12q)` has no polynomial
solution `Z`, for arbitrary `a,b,c,v`.

Write `Z=z_N(y)x^N+...`.  If `N>=1`, the coefficient of `x^(N+1)` in
`(12q)` is

```text
(v z_N)'=(N+2)c(v z_N).                                 (12ad)
```

The nonzero polynomial `v z_N` cannot satisfy this logarithmic-derivative
equation.  Hence `Z=z(y)`.  The coefficient of `x` in `(12q)` is then

```text
(vz)'-2c(vz)=2c.                                        (12ae)
```

Setting `r=vz+1` turns `(12ae)` into `r'=2cr`, so polynomial degree forces
`r=0`.  But `vz=-1` is impossible for nonconstant `v`.  This proves the
claim, including `Z=0`.

This separated cell was a useful hostile stepping stone.  The apparent
escape directions listed by its proof do not survive the complete argument
below.

### 3.7 All-degree closure of every five-word chart

**Theorem.**  Every integrable elementary Gauss word of length at most five
over a characteristic-zero field is tame.

It suffices to close `(12p)`.  If `W_x=0`, its rightmost factor is a source
shear and cancels, leaving a four-word map.  Suppose

```text
e=deg_x W>=1,
m=deg_x U,              n=deg_x V,
```

and first take nonzero `U,V`.  Put `d=m+n`.  When `d>=1`, let
`u_m,v_n,w_e` be the leading `x`-coefficients and let a hypothetical nonzero
`Z` have `x`-degree `N` and leading coefficient `z_N`.

If `n+N>0`, the coefficient above the degree of `-delta` is, for `e=1`,

```text
(u_m v_n z_N)'-(d+N+1)w_1u_mv_nz_N,                    (12af)
```

and, for `e>=2`,

```text
-(e+d+N)w_eu_mv_nz_N.                                  (12ag)
```

The first is an impossible polynomial logarithmic-derivative equation; the
second is plainly nonzero.  The boundary `Z=0` would force `delta=0`, which
Section 2.1 excludes when `W_x!=0`.  Hence any solution in the `d>=1`
branch must have `n=N=0`.

Now `V=v(y)` and `Z=z(y)`.  At the top `x`-degree, the `e=1` equation is

```text
[u_m(vz+1)]'=(m+1)w_1u_m(vz+1),                         (12ah)
```

whereas for `e>=2` it is simply `vz=-1`.  Equation `(12ah)` gives the same
conclusion.  Thus `v,z` are nonzero constants.  But the exact first
correction equation for constant `V=v,Z=z` is

```text
(zv+1)delta-zW_x=0,                                    (12ai)
```

which reduces to `zW_x=0`, a contradiction.  Therefore `d>=1` has no
solution.

It remains that `U=u(y),V=v(y)`.  Set `K=1+uv`.  If `K` is nonzero, the same
top comparison eliminates every `x`-dependent `Z`; for `Z=z(y)` it gives

```text
Kz=-u.                                                   (12aj)
```

For `e=1`, `(12aj)` first appears as
`(Kz+u)'=w_1(Kz+u)`; for `e>=2` it is the leading algebraic coefficient.
Since `gcd(K,u)=1`, equation `(12aj)` forces `K` to be a unit.  Nonzero
`u,v` are therefore constants.  Writing them again as `u,v`, if
`K=1+uv!=0`, the complete solution and final closure are

```text
Z=-u/K,                  R_5=-K^2W+S(y).                 (12ak)
```

The final rows are `(K,v+S/K)` and `(0,1/K)`, hence triangular.  The
load-bearing cancellation boundary `K=0` must be split off before the top
argument: here `uv=-1`,

```text
Z=u^2W+L(y),             epsilon=0,
```

and Section 3.4 makes the fifth factor a target shear after a four-word
automorphism.  Zero `U` or `V` merges adjacent factors.  Finally, the
opposite alternating orientation is obtained by simultaneous source/target
coordinate exchange, and adjacent equal signs combine.  This proves the
theorem.

### 3.8 The six-word correction system and its all-degree closure

The next alternating word is

```text
M_6=E_-(S_6)E_+(R_5)E_-(Z)E_+(V)E_-(U)E_+(W).          (12al)
```

Keep `A,B,epsilon,delta` from Section 3.4 and define

```text
C=ZA+U,                  D=ZB+1+UW,
kappa=C_y-D_x
     =delta+A Z_y-B Z_x+epsilon Z.                      (12am)
```

Unlike at length five, `kappa` need not vanish: the sixth factor can repair
the lower row after the fifth factor has repaired the top.  For fixed
`U,V,W,Z`, put

```text
E=A+R_5C,                F=B+R_5D.
```

The top row closes exactly when the first linear correction equation holds:

```text
C (R_5)_y-D (R_5)_x+kappa R_5=-epsilon.                 (12an)
```

Once `(12an)` holds, the final lower row `(C+S_6E,D+S_6F)` closes exactly
when

```text
E (S_6)_y-F (S_6)_x=-kappa.                             (12ao)
```

These are linear PDEs in `R_5` and then `S_6`, after the four base
parameters `U,V,W,Z` have been fixed.  If `kappa=0`, the five-word prefix has
both rows closed and is an automorphism by Section 3.7; equation `(12ao)`
then makes `S_6` a final target shear.  Thus a counterexample at length six
would have to satisfy

```text
kappa != 0.                                               (12ap)
```

The apparent branch `(12ap)` also closes without degree or coefficient
bounds.

**Theorem.**  Every integrable elementary Gauss word of length at most six
over a characteristic-zero field is tame.

It suffices to treat `(12al)`.  Three reductions first remove the source and
constant-factor boundaries.

**Source-shear boundary.**  If `W in k[y]`, choose `H'=W` and let
`T(x,y)=(x+H(y),y)`, so `DT=E_+(W)`.  Composing the integrated map with
`T^{-1}` cancels the rightmost factor; the other five parameters are merely
pulled back through a polynomial coordinate change.  Section 3.7 applies.

**Constant middle factor.**  If `Z=0`, the adjacent factors `E_+(R_5)` and
`E_+(V)` merge, leaving four factors.  If `Z=z in k*`, put

```text
w_z=[0,-1/z;z,0].
```

Direct matrix multiplication gives the exact Bruhat identity

```text
E_+(R_5)E_-(z)E_+(V)
 =E_+(R_5+1/z) w_z E_+(V+1/z),                           (12ap1)
```

and conjugation gives

```text
w_z E_+(f) w_z^-1=E_-(-z^2 f),
w_z E_-(f) w_z^-1=E_+(-f/z^2).                          (12ap2)
```

Pushing `w_z` through the three factors on its right rewrites the full word
as

```text
E_-(S_6) E_+(R_5+1/z) E_-(-z^2(V+1/z))
E_+(-U/z^2) E_-(-z^2W) w_z.                             (12ap3)
```

The five displayed variable factors alternate.  Absorb the constant right
factor by a linear source change with derivative `w_z^-1`; after the
corresponding argument pullback, one has an integrable five-word map, hence
an automorphism by Section 3.7.

It remains to assume `Z` is nonconstant and `W_x!=0`.  A zero value of `U`
or `V` merges adjacent factors, so take both nonzero and set `A=1+UV` as
above.

**The `A!=0` branch.**  Order monomials lexicographically with `x>y`.  Since
`W_x!=0`, write `LM(W)=x^e y^f` with `e>=1`.  If `A` is nonconstant then
`LM(A)=LM(UV)`; if `A` is a nonzero constant, then `U,V` are constants.
In either case there is no leading cancellation in

```text
LM(B)=LM(WA),
LM(C)=LM(ZA)>LM(U),
LM(D)=LM(ZWA)>LM(UW),1.                                 (12ap4)
```

Equation `(12an)` is equivalently the divergence equation

```text
(C R_5)_y-(D R_5)_x=-epsilon.                           (12ap5)
```

The leading term of `epsilon=A_y-B_x` is `-partial_x LM(B)`: compared with
an `A_y` term, the `B_x` term gains the lexicographically positive ratio
`x^(e-1)y^(f+1)`.  If `R_5` is nonzero, the same comparison makes the leading
term of the left side of `(12ap5)` equal to
`-partial_x LM(DR_5)`.  It strictly exceeds the leading term on the right by
the factor `LM(ZR_5)>1`.  This is impossible.  The boundary `R_5=0` is also
impossible because `epsilon!=0`.

**The cancellation `A=0`.**  Here `UV=-1`, so `U=u,V=v` are constants with
`uv=-1`.  The rows and curvatures simplify to

```text
A=0,  B=v,  C=u,  D=vZ+1+uW,
epsilon=0,             kappa=-D_x.                      (12ap6)
```

If `R_5=0`, the top row is `(0,v)`.  Final lower closure is
`-v(S_6)_x=D_x`, hence

```text
S_6=uD+T(y),
```

and the final rows are `(0,v)` and `(u,vT(y))`, an affine-triangular map.
If `R_5!=0`, top closure is

```text
u(R_5)_y=(D R_5)_x,
(R_5)_y=((D/u)R_5)_x.                                  (12ap7)
```

This is exactly the lower-curl equation from Section 2.1 with parameters
`R_5` and `D/u`.  Since `R_5` is nonzero, that lemma forces `(D/u)_x=0`.
Thus `kappa=0`; the five-word prefix has both rows closed and is tame by
Section 3.7, while final closure makes `S_6` a target shear in its polynomial
coordinates.

The opposite alternating orientation follows by simultaneous source/target
coordinate exchange, and adjacent factors of the same sign combine.  Hence

```text
every integrable elementary Gauss word of length at most six is tame.    (12aq)
```

The `10,000`-system six-word computation in Section 5.2 is therefore a finite
hostile replay inside this theorem, not evidence for an unclosed
`kappa!=0` cell.

### 3.9 All finite elementary depth closes: Wright's criterion

The six-factor proof exposes an induction that removes the apparent
length-seven frontier entirely.

**Cited theorem; independently reproduced direction.**  Let `k` be a
characteristic-zero field and let `F in k[x,y]^2` have `det DF=1`.  If
`DF in E_2(k[x,y])`, then `F` is a tame polynomial automorphism.  Wright's
1978 weak Jacobian theorem proves the stronger equivalence; the following
reduced-word proof independently audits the direction used by this search.

Merge adjacent elementary factors of the same sign.  Suppose first that all
remaining parameters are nonconstant and that the rightmost factor is
`E_+(f_1)`.  The two rows of an alternating product are consecutive vectors
in the continuant recurrence

```text
v_0=(0,1),  v_1=(1,f_1),  v_j=v_(j-2)+f_j v_(j-1).     (12aq1)
```

Under lexicographic order `x>y`, the product term in `(12aq1)` strictly
dominates the skipped term, so for `j>=2`

```text
LM(v_j,1)=product_(i=2)^j LM(f_i),
LM(v_j,2)=product_(i=1)^j LM(f_i).                      (12aq2)
```

The same ratio holds for `v_1`.  If `(f_1)_x!=0`, write
`LM(f_1)=x^e y^r`, `e>=1`.  In every nonconstant row `(P,Q)` from
`(12aq1)`, `LM(Q)=LM(P)LM(f_1)`.  The leading monomial of `Q_x` strictly
exceeds every term of `P_y`: if `e>1` its `x` exponent is larger; if `e=1`
the `x` exponents tie and its `y` exponent is larger.  This also covers a
leading `P` independent of `y`, whose derivative simply disappears.  Thus
the row cannot be closed.

Consequently `f_1 in k[y]`.  Choose `A'=f_1` and
`T(x,y)=(x+A(y),y)`.  If `DF=N E_+(f_1)`, the chain rule says

```text
D(F o T^-1)(z)=N(T^-1 z),                              (12aq3)
```

an integrable elementary word one factor shorter.  Polynomial coordinate
substitution preserves nonconstancy.  A rightmost lower factor is treated by
exchanging `x,y` and using lexicographic order `y>x`.

Constants also shorten the word.  Zero parameters delete and endpoint
constants are affine source or target gauges.  For an internal `c in k*`,

```text
E_+(a)E_-(c)E_+(b)
 =E_+(a+c^-1) w_c E_+(b+c^-1),
w_c=[0 -c^-1;c 0],                                     (12aq4)
```

while

```text
w_c E_+(h) w_c^-1=E_-(-c^2h),
w_c E_-(h) w_c^-1=E_+(-h/c^2).                         (12aq5)
```

Push `w_c` to the right, absorb it by a linear source change, and obtain one
fewer elementary factor after polynomial parameter substitution.  Induction
ends at a constant Jacobian.  Every operation is affine or triangular, so the
integrated map is tame.

Thus no rational-height or word-length enlargement inside `E_2` can reach a
counterexample.  The live invariant is the class in
`SL_2(k[x,y])/E_2(k[x,y])`, not elementary depth.  Cohn's classical matrix

```text
[1+xy x^2;-y^2 1-xy]                                  (12aq6)
```

is the smallest canonical hostile: determinant one but not elementary.

### 3.10 The seven-word staged equations as a historical replay

The next alternating word is

```text
M_7=E_+(T_7)E_-(S_6)E_+(R_5)E_-(Z)E_+(V)E_-(U)E_+(W).  (12ar)
```

For fixed `U,V,W,Z,R_5`, retain `C,D,kappa` and put

```text
E=A+R_5C,             F=B+R_5D,
rho=E_y-F_x
   =epsilon+C(R_5)_y-D(R_5)_x+kappa R_5.                (12as)
```

Algebraically, `rho` need not vanish because the last factor can repair the
top row.  After the sixth factor define

```text
G=C+S_6E,             H=D+S_6F,
sigma=G_y-H_x
     =kappa+E(S_6)_y-F(S_6)_x+rho S_6.                  (12at)
```

The bottom row must already close before the final upper factor, so the first
linear correction equation is `sigma=0`, or

```text
E(S_6)_y-F(S_6)_x+rho S_6=-kappa.
```

Once it holds, the final top row `(E+T_7G,F+T_7H)` closes exactly when

```text
G(T_7)_y-H(T_7)_x=-rho.                                 (12au)
```

The superseded bounded search plan would fix a base quintuple
`(U,V,W,Z,R_5)`, solve the whole chosen `S_6` space in `(12at)`, and then
solve the whole `T_7` space in `(12au)`.  If `rho=0`, a solvable first stage
makes the six-word prefix integrable and tame by `(12aq)`; the last factor is
a target shear.  Consequently the first genuinely live branch is

```text
rho != 0.                                                (12av)
```

This is no longer a live counterexample branch.  Section 3.9 applies to the
entire seven-factor matrix, so every solution of `(12at)--(12au)` is tame
regardless of `rho`.  The staged equations are retained only as exact
transport-operator controls and as a record of the finite-depth route that
led to Wright's all-length boundary.

## 4. Exact finite universe and gauges

Everything in this section is over `QQ`.  Sections 3.3, 3.7, and 3.8 already
prove the four-, five-, and six-word conclusions without finite parameter
bounds.  The censuses are retained as reproducible exact profiles of the curl
equations and as banks of balanced, full-span, symmetry-breaking positive
controls; they are not the basis for the all-degree conclusions.

1. The Jacobian is normalized to one.  Source/target affine gauges may set
   `F(0)=0` and `DF(0)=I`.  In this atlas the stronger conditions
   `U(0)=V(0)=W(0)=Z(0)=0` are imposed explicitly, and they indeed make
   `DF(0)=I`.  They are **atlas restrictions**, not consequences of a proved
   zero-constant Gauss normal form: constant elementary factors can cancel
   or refactor.
2. The small nonzero coefficient bank is

   ```text
   B={-2,-1,-1/2,1/2,1,2}.
   ```

3. The rightmost direction is one of

   ```text
   W in {y,x,x+y,x-y,x+2y,x-2y,2x+y,2x-y}.
   ```

   These are exactly the primitive projective integer directions of height at
   most two, up to overall sign.  Projectivizing the overall nonzero scale is
   legitimate for the geometric counterexample predicate.  For the scalar
   homothety

   ```text
   F_lambda(z)=lambda^-1 F(lambda z),
   DF_lambda(z)=DF(lambda z),
   ```

   a linear parameter `W=c ell` becomes `c lambda ell`; take
   `lambda=c^-1`.  The other parameter coefficients move under this pullback,
   so this is only a geometric gauge slice and not an equality of raw
   rational-height boxes.
4. `U` has support one, two, or three in

   ```text
   {x,y,x^2,xy,y^2},
   ```

   and every active coefficient belongs to `B`.  There are exactly
   `5*6+binom(5,2)*6^2+binom(5,3)*6^3=2,550` choices, hence
   `20,400` pairs `(U,W)`.
5. For every pair, equation `(11)` is solved over the **entire** rational
   vector space `deg V<=4`, with the explicit atlas condition `V(0)=0`.
   Coefficients of `V` are not sampled.
6. Every surviving equation `(12)` is solved over the **entire** rational
   vector space `deg Z<=8`, with `Z(0)=0`.  The constant in `Z` is also a
   left constant lower shear, hence a target `GL_2` redundancy.
7. Only after classifying these solution spaces are their free scalars sampled
   in `B` (and `B union {0}` for `Z`) to make the small-rational point bank.
8. As a separate length-five hostile control, `U,V,W` are independently
   selected from

   ```text
   {+/-x,+/-y,+/-x^2,+/-xy,+/-y^2},
   ```

   giving `10^3=1,000` base triples.  Equation `(12q)` is solved over the
   **entire** rational space `deg Z<=6`, this time including the constant
   coefficient.  The fifth parameter `R_5` is not sampled because every
   system dies at this first correction stage.  Section 3.5 proves the
   stronger monomial theorem, while Section 3.7 closes the entire five-word
   chart.
9. As a separate length-six hostile control, `U,V,W,Z` are independently
   selected from the same ten signed nonconstant monomials, giving `10^4`
   base quadruples.  Equation `(12an)` is solved over the **entire** rational
   space `deg R_5<=6`, constant included, of dimension `28`.  Every system
   has coefficient/augmented rank `(28,29)` modulo `1009`, so all are
   inconsistent over characteristic zero.  The `S_6` stage is not reached.
   Section 3.8 proves the stronger all-degree six-word theorem.

The finite four-word histograms below make no claim about a four- or
five-term `U`, `deg V>=5`, or their raw rational-height boxes.  The finite
five- and six-word controls make no bounded claim about multi-term
parameters, nonzero intermediate constants, or another correction space.
The proved all-degree theorems above retain their separately stated
quantifiers.

## 5. Complete first-curl result

The exact hybrid rank computation gives

| `W` | inconsistent | only `V=0` | nonzero one-dimensional kernel |
|---|---:|---:|---:|
| `y` | 0 | 2,546 | 4 |
| `x` | 2,550 | 0 | 0 |
| `x+y` | 2,550 | 0 | 0 |
| `x-y` | 2,550 | 0 | 0 |
| `x+2y` | 2,550 | 0 | 0 |
| `x-2y` | 2,550 | 0 | 0 |
| `2x+y` | 2,550 | 0 | 0 |
| `2x-y` | 2,550 | 0 | 0 |

The `2,546` zero-only rows have `A=1` and therefore fail the nonunit edge and
balanced-factor passport.  Every `x`-bearing direction dies before the
second curl.  The four nonzero families are exactly

```text
alpha in {-2,-1,1,2},
s=x+y^2/2,
U=alpha s,
q=y+alpha s^2/2,
V=t q.                                                    (13)
```

They all have `delta=0`.  Thus, after the proved unit-edge rejection of the
zero-only rows, the primary statistic is

```text
delta-nonzero nonunit-eligible first-row survivors = 0.  (14)
```

This is not inferred from a random sample.  Full coefficient-column rank
modulo `1009` soundly certifies all `17,850` inconsistent systems and all
`2,546` zero kernels: full modular rank forces full rational rank, and an
augmented rank of `16` against coefficient rank `15` forces inconsistency.
Every modular rank drop is replayed by exact python-FLINT rational row
reduction; exactly the four displayed kernels survive.

### 5.1 Length-five signed-monomial hostile replay

For the separate universe in item 8 of Section 4, all `1,000` first
correction systems `(12q)` are inconsistent in the complete space
`deg Z<=6`, constant included.  Thus no `R_5` system and no collision test is
reached.  This bounded replay is an implementation control for the
all-exponent monomial proof in Section 3.5, not the source of that theorem.

### 5.2 Length-six signed-monomial hostile replay

For the separate universe in item 9 of Section 4, all `10,000` top-closure
systems `(12an)` are inconsistent in the complete 28-dimensional space
`deg R_5<=6`, constant included.  Modulo `1009`, every coefficient/augmented
rank pair is `(28,29)`; maximal coefficient rank makes this a sound
characteristic-zero certificate.  Thus no `S_6` system is reached.  The
two-scale extension takes `U,V,Z` from the same bank,
`W=y^2+/-x`, and the entire 45-dimensional space `deg R_5<=8`; all `2,000`
systems have ranks `(45,46)` and are likewise inconsistent.  The companion
also checks a nonzero-`kappa`, all-factor identity control and a height-two
nonconstant control integrating to an explicit tame automorphism.  These
bounded replays are implementation controls for Sections 3.8--3.9, not the
source of the all-length theorem.

## 6. Complete second-curl result and collision consequence

For each of the four values of `alpha` and each nonzero `t in B`, let

```text
p=s+t q^2/2.                                               (15)
```

The complete `deg Z<=8`, `Z(0)=0` solution space of `(12)` is

```text
Z=z p,                    z in QQ.                        (16)
```

The `24` second-row systems all have this one-dimensional kernel.  Sampling
`z in B union {0}` produces `4*6*7=168` exact small-rational maps

```text
F_(alpha,t,z)=(p,r),       r=q+z p^2/2.                   (17)
```

Every one has coefficient-matrix span four.  More importantly, every one
has the global polynomial inverse

```text
p=P,
q=R-zP^2/2,
s=P-tq^2/2,
y=q-alpha s^2/2,
x=s-y^2/2.                                                 (18)
```

Thus the collision consequence is exact and global:

```text
counterexample survivors = 0,
global inverse certificates = 168.                       (19)
```

As a redundant hostile, the exact `49`-point rational bank
`(B union {0})^2` has no image collision for any of the `168` maps.  The
bounded point test is not the reason for `(19)`; formula `(18)` is.

## 7. The balanced span-four controls are genuine

When `z=0`, the `24` maps in `(17)` occupy the exact balanced distinct-factor
cell of THM-3587 and have coefficient span four.  Work in the polynomial
coordinates `(s,q)`, where `y=q-alpha s^2/2`.  Their derivative edges are

```text
c=alpha s,
a=1+alpha t s q,
d=1+alpha s q-alpha^2 s^3/2,
b=alpha t s q^2+(1+t-alpha^2 t s^3/2)q-alpha s^2/2.     (20)
```

The first three are absolutely irreducible: `c` is linear, while `a` and
`d` are primitive linear polynomials in `q`.  They are pairwise nonassociate
on their matching sides.

If `t!=-1`, `b` is a primitive quadratic in `q`.  Its discriminant is

```text
Delta(s)=(1+t-alpha^2 t s^3/2)^2+2 alpha^2 t s^3.        (21)
```

As a quadratic in `z=s^3`, the discriminant of `(21)` is

```text
-4 alpha^4 t^3 !=0.                                     (22)
```

Its two nonzero `z`-roots are distinct, so `Delta(s)` has simple roots and
is not a square over `C[s]`.  Hence `b` is absolutely irreducible.

At `t=-1`, instead

```text
b=alpha s[-q^2+(alpha/2)s^2q-s/2].                      (23)
```

The bracket is irreducible because its `q`-discriminant
`alpha^2s^4/4-2s` has a simple zero at `s=0`.  Consequently, even here
`T=bc` has exactly the two distinct factors `s` and the bracket.  In every
case `T+1=ad` also has exactly two distinct factors.  Therefore

```text
(omega(T),omega(T+1))=(2,2)                              (24)
```

for all `24` controls.  Exact coefficient matrices have rank four.  These
are useful positive controls: balanced factors, full span, small rationals,
and no obvious linear involution still coexist with explicit tameness.

No absolute factor count is asserted here for the `z!=0` maps.

## 8. Involution guardrail and what “symmetry-breaking” means here

There is a clean reduction from reflection-intertwining cells to a cited
closure theorem.
Suppose a complex Keller map satisfies

```text
F sigma=tau F,                                             (25)
```

where `sigma,tau` are affine reflections with fixed affine lines `L,M`.
Then `F(L)` lies in `M`.  Differentiate `(25)` along the `+1` tangent of
`L`.  The invertible matrix `DF` sends it to a nonzero `+1` tangent of
`M`, so the one-variable polynomial restriction `F|_L:L->M` has
nowhere-vanishing derivative.  Over `C` that derivative is constant;
therefore the restriction is affine and injective.  The **CITED**
[Gwozdziewicz injective-line theorem](https://arxiv.org/abs/alg-geom/9305008)
then makes `F` an automorphism.  This argument does
not identify all involutions and does not apply to the central involution,
whose fixed locus is a point rather than a line.

The script asks, for each map and each of five literal source involutions

```text
-I, diag(-1,1), diag(1,-1), exchange, negative exchange,
```

whether there is a rational linear target involution intertwining the map.
All `168` answers are negative.  Thus these maps escape the cheapest linear
equivariant closures but remain tame by `(18)`: symmetry breaking is
necessary as a search priority, not sufficient as a counterexample signal.

For the separate central commuting case, Miyanishi's
[*Equivariant Jacobian Conjecture in dimension two*](https://arxiv.org/abs/2110.06709)
is **CITED**: it proves the commuting effective finite-group case of even
order over `C`.  That is a commuting-action theorem, not permission to
identify arbitrary source and target intertwiners, and it is not a dependency
of the finite atlas because no named intertwiner survives.

## 9. Hostiles, failure boundary, and next exact box

The cheapest hostile is

```text
U=x,       W=x,       V=0.
```

The four-word determinant is still one, but the first curl defect is `-1`.
This demonstrates why determinant and displayed matching factors are only
the first filter.

The load-bearing hostile for the six-word leading-term proof is the constant
Bruhat cancellation.  For arbitrary `U,W` and `v in k*`, take

```text
V=v,  Z=-1/v,  R_5=v(1+vU),  S_6=W/v^2+T(y).
```

The final rows are `(0,v)` and `(-1/v,vT(y))`.  Thus the terms nominally led
by `ZUVW` can cancel, but the result is affine-triangular rather than a
counterexample.  The smallest exact instance

```text
U=W=x,  V=1,  Z=-1,  R_5=1+x,  S_6=x
```

has constant final matrix `[0,1;-1,0]`.  This is why Section 3.8 separates
constant `Z` before applying its lexicographic argument.

The present atlas does **not** say that rational coefficients cannot produce
a planar counterexample.  It says something more surgical:

1. By Wright's criterion, every integrable finite elementary Gauss word is
   tame, with no word-length, degree, or coefficient-height restriction.
2. In the normalized `20,400`-row four-word replay, every nonunit-eligible
   `delta!=0` attempt dies at the first curl.  All remaining sampled points are
   four sequential shears and hence
   have a printed inverse, even though they are balanced, span four, and
   evade all five literal involution tests.
3. The `1,000` signed-monomial five-word systems are all inconsistent at the
   first correction stage, a finite control inside the theorem.
4. In the exact six-word atlas, all `10,000` signed-monomial choices of
   `(U,V,W,Z)` make `(12an)` inconsistent for the complete 28-dimensional
   space `deg R_5<=6`: modulo `1009`, every coefficient/augmented rank pair is
   `(28,29)`.  Maximal coefficient rank makes this a characteristic-zero
   certificate.  The `S_6` stage is therefore never reached.  This is a
   finite control inside the all-degree theorem.  The added `2,000` two-scale
   bases against the full degree-eight correction space similarly all have
   rank pair `(45,46)`.
5. The seven-factor staged system is also tame for every `rho`; it is a
   transport-equation replay, not a counterexample frontier.

No enlargement inside `E_2(k[x,y])` can escape the theorem.  The next
decisive atlas should preserve a certified nonzero class in
`SL_2(k[x,y])/E_2(k[x,y])`.  The smallest canonical seed is Cohn's matrix
`[1+xy,x^2;-y^2,1-xy]`.  The exact next box is therefore:

- multiply the Cohn core on both sides by small elementary decorations,
  explicitly retaining its non-elementary double-coset class;
- build the monomial-incidence graph of both row curls and reject exposed
  factorial ladders before raising coefficient degree;
- retain only support cycles whose exact multiplier holonomy satisfies the
  cycle-product compatibility, with reciprocal small rational gains tested
  first;
- solve both complete correction spaces, not a sampled coefficient list; and
- retain positive tame controls and the literal involution filters, but use
  neither bounded point injectivity nor symmetry breaking as a conclusion.

The first object worth recording is not a determinant table but a
non-elementary matrix whose two rows are both closed.  Only then should
integration, factorization, full coefficient span, nonproperness, and exact
collision elimination consume the larger computational budget.  The detailed
router is in
[`jc2-cohn-parity-cycle-repair-codex-20260821.md`](jc2-cohn-parity-cycle-repair-codex-20260821.md).

## 10. Reproduction and scope ledger

Run

```bash
python3 04-computation/jc2_small_rational_four_gauss_search.py
python3 -O 04-computation/jc2_small_rational_four_gauss_search.py
python3 04-computation/jc2_small_rational_six_gauss_search.py
python3 -O 04-computation/jc2_small_rational_six_gauss_search.py
diff -u \
  05-knowledge/results/jc2_small_rational_four_gauss_search.out \
  <(python3 04-computation/jc2_small_rational_four_gauss_search.py)
diff -u \
  05-knowledge/results/jc2_small_rational_six_gauss_search.out \
  <(python3 04-computation/jc2_small_rational_six_gauss_search.py)
```

The program uses exact `Fraction` polynomial arithmetic, a sound full-rank
certificate modulo `1009`, and python-FLINT rational row reduction for every
rank drop.  It checks `1,331` load-bearing identities and consequences,
including both four-word curls, determinant one, all inverse certificates,
coefficient spans, literal involution tests, the rational collision bank, and
the `1,000` five-word first-correction systems.  `python3` and `python3 -O`
have identical logic.  The six-word companion independently checks `10,000`
base quadruples, all 28 middle-correction columns, `2,000` two-scale bases
against all 45 degree-eight columns, both maximal modular-rank certificates,
and two exact tame controls with nonzero staged defects.

Scope labels:

- Wright's all-length elementary-Jacobian criterion: **CITED + INDEPENDENTLY
  REPROVED**; the displayed length-three-through-six reductions and
  seven-word staged identities: **PROVED**;
- four-Gauss universe, histograms, kernels, spans, and collision bank, plus
  the `1,000` signed-monomial five-word control, `10,000` signed-monomial
  six-word degree-six systems, and `2,000` two-scale degree-eight systems:
  **FINITE-EXACT**;
- exchange/even-group equivariant filters: **CITED**, not dependencies;
- non-elementary Cohn-core correction cells, any larger rational-height box,
  and any `JC(2)` conclusion: **OPEN**.
