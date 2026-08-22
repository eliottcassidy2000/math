---
id: THM-3561
title: "Rational Keller triple collision and Danielewski polynomial completion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  middle-coordinate zero leaf of the fixed
  three-dimensional Keller map is the punctured plane
  A2 minus V(1-x^2 w).  In a symplectic punctured-plane chart its restriction
  is the rational Keller pair (q/D^2,xD(D+2)), D=1+x^2q, with determinant -3
  and an exact triple collision.  Its maximal subalgebra of polynomial target
  observables that extend across D=0 is the smooth Danielewski ring
  C[b,c,e]/(c^2e-b(b+4)).  The resulting polynomial map A2 to this surface is
  etale of generic degree three, retains the triple collision, and misses one
  point.  The target has Euler characteristic two and Picard group Z, so it is
  not A2.  Its global symplectic form is nevertheless exact and H^1_dR=0,
  so topology alone does not obstruct a Darboux pair.  The complete
  homogeneous-weight cell, all one-generator Hamiltonian slices, and the
  complete two-weights-by-two-weights cell are excluded.  A genuine planar
  counterexample through this completion would require at least two weights
  in each output and at least three in one of them, hence at least five
  nonzero homogeneous pieces in total.  No such pair is claimed.
source: root/dimension3_descent_hostiles, 2026-08-18
audit: >
  Independent hostile audits rederived the punctured-plane chart, maximal
  polynomial intersection, smooth completion, exact image, Picard class,
  symplectic primitive and de Rham complex.  A separate current-byte audit
  checked the homogeneous and complete two-by-two weight obstructions,
  including every zero-weight and R=T boundary.  Ordinary and optimized
  companion replays are byte-identical; active-gate, hash, documentation,
  and diff checks pass.
depends_on:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
related:
  - THM-3404-factorized-danielewski-principal-parts-and-finite-cover-obstruction
  - THM-3554-punctured-kummer-collision-surface-normal-form
  - THM-3559-affine-target-coordinate-pullback-no-go
script: 04-computation/jc2_rational_keller_danielewski_completion_thm3561.py
output: 05-knowledge/results/jc2_rational_keller_danielewski_completion_thm3561.out
script_sha256: 2ee110af0f66a024bc4cb74c39ad4b81cba9380a2d5c968f5eb4532d6603b746
output_sha256: f6153a9921c10f964031180f0312bbc02e247ac09352bdfd1af75737f96a5bae
hash_basis: LF-normalized bytes
---

# THM-3561 -- rational Keller collision and Danielewski completion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
does **not** prove `JC(2)` false.  It turns one descent attempt into a precise
near-counterexample and an exact residual coefficient system.

All rings and varieties below are over `C`.  Compactly supported Euler
characteristic is used in Section 5.

## 1. The middle-coordinate leaf is a punctured plane

Retain the fixed map `F=(F1,F2,F3)` of THM-1300 and put `u=1+xy`.  Its middle
coordinate has the two exact forms

```text
F2=(9u^2-6u-2)y+3xu^2z,                                (1)
xF2=2+u(3x^2zu+9u^2-15u+4).                           (2)
```

Therefore `u` is a unit on the smooth surface `S=V(F2)`, with explicit
inverse

```text
u^(-1)=-(3x^2zu+9u^2-15u+4)/2.                        (3)
```

Smoothness also follows immediately from the ambient Keller condition: the
row `dF2` of the invertible matrix `JF` never vanishes.

In the localization at `u`, set

```text
w=3z+2(4u+1)y^2/u^2.                                  (4)
```

Direct expansion gives

```text
F2=u^2(y+xw).                                          (5)
```

Thus on `S`,

```text
y=-xw,             u=1-x^2w,
z=(w-2(4u+1)x^2w^2/u^2)/3,                             (6)
```

and these formulas are mutually inverse ring maps:

```text
C[S] ~= C[x,w,(1-x^2w)^(-1)].                          (7)
```

The omitted curve `1-x^2w=0` is smooth and isomorphic to `G_m`.  Since the
polynomial is irreducible, localization of the UFD `C[x,w]` also gives the
complete unit group

```text
C[S]^* = C^* u^Z.                                      (8)
```

In particular `S` is not `A2`.  More strongly, every morphism `A2->S` pulls
`u` back to a constant unit.  Hence no such morphism is dominant, and no
polynomial `A2` chart can contain both the fixed collision point, where
`u=1`, and either orbit point, where `u=-1/2`.

## 2. A rational planar Keller pair with a triple collision

On `S`, direct substitution into the other two target coordinates gives

```text
F1=uw/3,              F3=2x(1+2u)/(3u^2).              (9)
```

Put

```text
q=w/u,                D=1+x^2q=u^(-1),
a=3F1,                c=3F3/2.                         (10)
```

Then `S` is the principal open set `D!=0` in the `(x,q)`-plane, and the
restricted map is

```text
a=q/D^2,              c=xD(D+2).                       (11)
```

It has the exact standard-plane determinant

```text
Jac_(x,q)(a,c)=-3.                                     (12)
```

The three original collision points become

```text
r0=(0,-3/4),          r+=(1,-3),          r-=(-1,-3),  (13)
```

with `D`-values `1,-2,-2`, respectively, and

```text
(a,c)(r0)=(a,c)(r+)=(a,c)(r-)=(-3/4,0).                (14)
```

Thus `(11)` is a literal constant-Jacobian, noninjective rational plane map.
It is regular on a punctured plane, not a polynomial endomorphism of `A2`:
`a` has a double pole on `D=0`.

## 3. The exact polynomial-observable intersection

View `C[a,c]` as a subring of `C(x,q)` through `(11)`.  The largest part of
this **polynomial target ring** whose pullback extends across `D=0` is

```text
C[a,c] intersect C[x,q]
  = C[b,c,e],                                           (15)

b=ac^2=(D-1)(D+2)^2,
e=a(b+4)=q(D+3),                                       (16)
```

and its only defining relation is

```text
c^2e=b(b+4).                                           (17)
```

The qualifier in `(15)` is important: this is an intersection with the
polynomial ring `C[a,c]`, not a claim about every rational target function.

### Proof of `(15)`

Give `x,q` weights `1,-2`.  Then `D,b` have weight zero, while `a,c,e` have
weights `-2,1,-2`.  Both rings in `(15)` are graded, so membership can be
tested one homogeneous weight at a time.

For weight `r>=0`, every homogeneous element of `C[a,c]` has the unique form

```text
c^r P(b),              P in C[b],                      (18)
```

which is already in `C[b,c]`.  For weight `-s<0`, it has the form

```text
c^(-s)P(b),                                             (19)
```

where polynomiality in `a,c` first forces

```text
b^m divides P(b),             m=ceil(s/2).              (20)
```

Indeed a monomial of weight `-s` is `a^i c^(2i-s)=b^i
c^(-s)`, and `2i-s>=0` is exactly `i>=m`.

Substitute `(16)` into `(19)`.  After using
`D-1=x^2q`, condition `(20)` cancels every possible negative power of `x`
and `D+2`; the only remaining denominator is `D^s`.  But

```text
b+4=D^2(D+3).                                          (21)
```

Since `D` is prime and `D+3` is a unit at its generic point, cancellation of
`D^s` is equivalent to

```text
(b+4)^m divides P(b).                                  (22)
```

Together `(20)` and `(22)` make `(19)` equal to

```text
e^(s/2) P0(b),                 s even,
c e^((s+1)/2) P0(b),           s odd,                  (23)
```

for some `P0 in C[b]`.  These belong to `C[b,c,e]`.  Conversely all three
generators in `(16)` are visibly polynomial in `(x,q)`, proving `(15)`.

Equation `(17)` follows directly.  The functions `b,c` are algebraically
independent because

```text
Jac_(x,q)(b,c)=-3c^2                                  (24)
```

is generically nonzero.  Hence `(17)` is the full kernel, not merely one
relation among additional hidden ones.  This proves the intersection
statement.

## 4. Polynomial completion on a smooth Danielewski surface

Let

```text
Y=Spec C[b,c,e]/(c^2e-b(b+4)).                         (25)
```

The gradient of its defining equation is

```text
(-(2b+4), 2ce, c^2).                                   (26)
```

If it vanished, then `c=0,b=-2`, contradicting `(25)`.  Thus `Y` is smooth.

Equations `(16)` define a polynomial morphism

```text
Phi:A2_(x,q) -> Y,
Phi=( (D-1)(D+2)^2, xD(D+2), q(D+3) ).                 (27)
```

It extends the rational Keller pair through the missing divisor `D=0`.

### 4.1 `Phi` is everywhere etale

On `c!=0`, the pair `(b,c)` is a coordinate chart on `Y`, and

```text
Jac_(x,q)(b,c)=-3c^2!=0.                               (28)
```

On `c=0`, equation `(25)` gives `b=0` or `b=-4`; there `(c,e)` is a target
chart.  The source equation

```text
c=xD(D+2)=0                                            (29)
```

has exactly the three cases `x=0` (hence `D=1`), `D=0`, and `D=-2`.  Exact
differentiation gives

```text
Jac_(x,q)(c,e)=6(D+1)(D^2+2D-2),                       (30)
```

whose values in those three cases are, respectively,

```text
12,                  -12,                  12.         (31)
```

The differential is therefore an isomorphism everywhere, proving that
`Phi` is etale.

### 4.2 Generic degree and image

The zero-weight parameter `D` satisfies

```text
b=(D-1)(D+2)^2,
D^3+3D^2-(b+4)=0.                                     (32)
```

The polynomial in `(32)` is irreducible over `C(b,c)`: it is primitive and
the quotient of `C[b,c,D]` by it is the domain `C[c,D]`, after eliminating
`b`.  Generically,

```text
x=c/(D(D+2)),                 q=(D-1)/x^2,             (33)
```

so `Phi` has generic degree exactly three.

Its geometric image is also exact:

```text
Phi(A2)=Y minus {(-4,0,0)}.                             (34)
```

First, a preimage of `(-4,0,0)` would satisfy
`b+4=D^2(D+3)=0`.  If `D=0`, then `e=3q` and
`x^2q=-1`, so `e!=0`.  If `D=-3`, then `c=3x`, so `c=0`
would force `x=0` and hence `D=1`.  Both are contradictions.

Conversely, if `b+4!=0`, put `a=e/(b+4)`.  The relation on `Y` says
`b=ac^2`.  If `b!=0`, any root in `(32)` avoids `0,-2`, and `(33)`
reconstructs a preimage.  If `b=0,c!=0`, choose
`D=1,x=c/3,q=0`; if `b=c=0`, choose `D=1,x=0,q=a=e/4`.
If `b=-4,c!=0`, choose `D=-3`,
`x=c/3`, and `q=-36/c^2`.  Finally, if `b=-4,c=0,e!=0`, choose
`D=0`, `q=e/3`, and either square root `x^2=-3/e`.  This covers every point
except the one in `(34)`.

### 4.3 The triple collision survives

At all three points in `(13)`, direct substitution in `(27)` gives

```text
Phi(r0)=Phi(r+)=Phi(r-)=(0,0,-3).                       (35)
```

Thus `Phi` is an explicit noninjective polynomial etale map from the actual
affine plane to a smooth affine surface.  It is not yet a planar Jacobian
counterexample because that target surface is not an affine plane.

This is precisely a nonfinite escape of the type left open by THM-3404.
Introduce the ordinary two-arm Danielewski model

```text
Y0: cd=b(b+4).                                          (35a)
```

The morphism `Y->Y0` is induced by `d=ce`; its ring inclusion is not finite.
THM-3404 proves that the free arm class of `Y0` cannot die under a finite
dominant normal factorial cover.  The present construction uses the typed
escape that theorem explicitly leaves open: first the nonfinite affine
modification `Y->Y0`, and then the generically finite but nonfinite `Phi`,
whose image misses the point in `(34)`.  Thus there is no conflict with the
finite-cover obstruction; `Phi` realizes its missing-point/nonproper
boundary case.

## 5. The target obstruction is topological and divisor-theoretic

On `c!=0`, equation `(25)` eliminates `e`, giving

```text
Y minus V(c) ~= G_m x A1.                              (36)
```

On `c=0`, the target is the disjoint union of the two affine lines

```text
L0=V(c,b),               L4=V(c,b+4).                  (37)
```

Therefore

```text
chi_c(Y)=0+1+1=2,             chi_c(A2)=1.             (38)
```

In particular `Y` is not isomorphic to `A2`.

There is also an exact Picard obstruction.  Localizing at `c` gives the UFD

```text
C[Y][c^(-1)] ~= C[c,c^(-1),b].                         (39)
```

Nagata localization says that the class group is generated by `L0,L4`.
The localized unit group is `C^*c^Z`, and

```text
div(c)=L0+L4.                                          (40)
```

This is the only relation, so

```text
Cl(Y)=Z^2/<(1,1)> ~= Z.                                (41)
```

Because `Y` is smooth, `Pic(Y)=Cl(Y)=Z`, while `Pic(A2)=0`.  The completion
has moved the obstruction from the nonconstant unit on the source leaf to a
nontrivial divisor class on the target.

## 6. Exact residual coefficient system for `JC(2)`

The bracket inherited from the standard form `da wedge dc` extends over
`Y`.  On the generators it is

```text
{b,c}=c^2,             {c,e}=-2(b+2),
{b,e}=-2ce.                                             (42)
```

Consequently, for `P,Q in C[Y]`,

```text
{P,Q}=
 c^2(P_bQ_c-P_cQ_b)
-2ce(P_bQ_e-P_eQ_b)
-2(b+2)(P_cQ_e-P_eQ_c),                                (43)
```

computed modulo `(17)`.  Every element of `C[Y]` is, by `(15)`, both a
polynomial in the target variables `(a,c)` and a polynomial after pullback
to `(x,q)`.  Hence

```text
Jac_(x,q)(P(Phi),Q(Phi))=-3{P,Q}(Phi).                  (44)
```

This produces a sharp construction target:

```text
find P,Q in C[b,c,e]/(c^2e-b(b+4)) and kappa in C^*
such that {P,Q}=kappa.                                  (45)
```

Any solution of `(45)` would make `(P(Phi),Q(Phi))` a polynomial planar
Keller pair with determinant `-3kappa`; equation `(35)` would make it
noninjective.  It would therefore be a genuine counterexample to `JC(2)`.

No solution of `(45)` is known or claimed here.  The obvious compression
`(b,c)` fails by `(24)`, and adding `e` removes the universal `c^2` factor
only at the price of the nonconstant bracket `-2(b+2)`.  Thus the
Danielewski completion is a concrete residual coefficient system, not a
disguised proof of the planar conjecture.

The affine-linear cell is already empty.  If `P,Q` are affine-linear in
`b,c,e`, their bracket is a linear combination of `c^2`, `ce`, and `b+2`.
Killing the `b` coefficient also kills the only constant term, so a constant
bracket in this cell is necessarily zero.  This is a bounded exact hostile,
not an extrapolation to nonlinear `P,Q`.

### 6.1 Exact symplectic form: a failed topological obstruction

The Poisson tensor in `(42)` is everywhere nondegenerate.  On `c!=0`, its
inverse symplectic form is

```text
omega=db wedge dc/c^2.                                  (46)
```

It is tempting to use `Pic(Y)!=0` or `H^2_(dR)(Y)!=0` to declare `(46)`
nonexact and thereby exclude `(45)`.  That inference is false.  Put

```text
alpha=-(ce/4)db + (e(b+2)/4)dc + (c(b+2)/8)de.          (47)
```

Differentiating `b(b+4)-c^2e=0` and using
`(b+2)^2=c^2e+4` gives, in the differential module of `Y`,

```text
c alpha=db,             hence alpha=db/c on c!=0,
d alpha=omega.                                                (48)
```

Thus the physical symplectic form passes the global exactness gate.  The
passing mechanism is the two-arm cancellation in `(47)`, not the deletion of
the divisor `c=0`.

The complete algebraic de Rham calculation makes the boundary precise.  Give

```text
wt(b)=0,                  wt(c)=1,                  wt(e)=-2.   (49)
```

The Euler field contracts every nonzero-weight subcomplex.  In weight zero,
write

```text
S=b(b+4),         eta=ce dc,         theta=c omega=db wedge dc/c.
```

Then the de Rham complex reduces to

```text
Omega^0_0=C[b],
Omega^1_0={f(b)db+g(b)eta},
Omega^2_0=C[b]theta,

d(f db+g eta)=(Sg)' theta.                             (50)
```

Consequently

```text
H^1_dR(Y)=0,
H^2_dR(Y)=C[b]/{(Sg)':g in C[b]} ~= C.                 (51)
```

For example, if `H'=h`, the functional `h |-> H(0)-H(-4)` kills every
`(Sg)'` and sends `1` to `4`.  Hence `[omega]=0` while `[c omega]!=0`.

Since `H^1_dR(Y)=0`, equation `(45)` with `kappa=1` is equivalent to the
more pointed primitive equation

```text
alpha=P dQ+d phi,                 P,Q,phi in C[Y].      (52)
```

There is already a rational solution,

```text
{b,-1/c}=1,
alpha=b d(-1/c)+d(b/c).                                 (53)
```

The remaining debt is exactly polynomial clearing of the pole `1/c`, not
symplectic exactness.

### 6.2 All homogeneous Darboux pairs fail

The bracket raises the weight `(49)` by one.  If homogeneous `P,Q` solve
`{P,Q}=1`, their weights are `r` and `-r-1`; after swapping, take `r>=0`.
The all-weight intersection proof in Section 3 gives, on `c!=0`,

```text
P=c^r f(b),
Q=c^(-r-1) S(b)^m g(b),       m=ceil((r+1)/2).          (54)
```

Writing `h=S^m g`, direct use of `{b,c}=c^2` gives

```text
{P,Q}=-(r+1)f'h-r f h'.                                 (55)
```

For `r=0`, `(55)` is divisible by `S`.  For `r>=2`, it is divisible by
`S^(m-1)`.  For `r=1`, it is

```text
-2f'Sg-f(Sg)',                                          (56)
```

whose leading coefficient, for nonzero `f,g`, is

```text
-(2 deg(f)+deg(g)+2) lc(f)lc(g),                        (57)
```

at positive degree.  No case is a nonzero constant.  If only one of two
candidate outputs were homogeneous, its bracket with the other's weight
decomposition would isolate one complementary homogeneous component at
weight zero, giving the same contradiction.  Therefore every solution of
`(45)` must mix at least two weights in **both** outputs, and its constant
term must be a cancellation among at least two complementary weight pairs.

The three generator slices fail for sharper reasons.  The image of
`{b,-}` lies in `(c)`.  The locally nilpotent derivation `{c,-}` has kernel
`C[c]`, but a slice would make `C[Y]=C[c,Q]`, contradicting the reducible
special fibre

```text
C[Y]/(c)=C[b,e]/(b(b+4)).                               (58)
```

Finally, on `e=1`, the derivation `{e,-}` becomes `2r d/dr` on a Laurent
coordinate ring `C[r,r^(-1)]`, whose image has zero constant Laurent
coefficient.  It cannot contain `1`.  These are structural hostiles, not a
general no-go for mixed nonlinear pairs.

### 6.3 The complete two-by-two weight cell is empty

Put `S=b(b+4)`.  The homogeneous piece of weight `t` is

```text
A_t={c^t f(b): f in C[b],
     S^ceil(-t/2) divides f when t<0},                  (59)
```

and, on `c!=0`,

```text
{c^r f,c^s g}=c^(r+s+1)(s f'g-r f g').                (60)
```

Suppose that each of `P,Q` has exactly two nonzero weight pieces and that
`{P,Q}=1`.  Subtract constants from `P,Q` first.  If this removes a
weight-zero piece, one output becomes homogeneous and Section 6.2 applies;
thus every retained weight-zero coefficient is nonconstant.  Write the
supports as

```text
r0<r1,                 s0<s1.                          (61)
```

A single complementary pair `r_i+s_j=-1` would itself be a forbidden
homogeneous Darboux pair.  Hence at least two pairs must contribute to weight
zero.  Ordering the four sums forces the unique cross matching

```text
r0+s1=r1+s0=-1,
delta:=r1-r0=s1-s0>0.                                  (62)
```

The two extreme brackets, of output weights `-delta` and `delta`, must
vanish.  A useful root-valuation lemma fixes their signs.  If `p<0<q`, then
at either root of `S` the leading Wronskian coefficient of a putatively
vanishing bracket is

```text
q ord(f)-p ord(g) != 0.                                (63)
```

If `p<0=q`, formula `(60)` reduces to `-p f g'`; it can vanish only when the
retained weight-zero coefficient `g` is scalar, already removed above.  The
case `p=0>q` is symmetric.  Therefore the low extreme has both weights
strictly negative and the high extreme has both weights nonnegative.  There
are integers `R,T>=1` and nonzero polynomials `f,g,F,G` such that

```text
P=c^(-R)f+c^(T-1)F,
Q=c^(-T)g+c^(R-1)G,                                    (64)
```

where `S^ceil(R/2)|f` and `S^ceil(T/2)|g`.  The two extreme equations and
the weight-zero coefficient are

```text
Rfg'-Tf'g=0,
(R-1)F'G-(T-1)FG'=0,                                   (65)

H=(R-1)f'G+RfG'-TF'g-(T-1)Fg'=1.                       (66)
```

At either root `beta` of `S`, the middle two terms of `(66)` vanish.  The
first can survive only when `R=2`, and the last only when `T=2`: order one
is necessary, while the coefficient kills the boundary value `R=1` or
`T=1`.  Hence `R=2` or `T=2`.

Assume `R=2`.  The first equation in `(65)`, together with `H(beta)=1` on
both arms, forces `T=2k` even and forces `f` to have simple zeros there.  The
case `T=2` is included: the same Wronskian makes `f` and `g` simultaneously
simple even if the last term of `(66)` supplies the nonzero arm value.  UFD
comparison in the two vanishing Wronskians now gives

```text
f=A h,             g=B h^k,             S|h,
F=L K^(2k-1),      G=M K,                                (67)
```

for nonzero constants `A,B,L,M` and nonzero polynomials `h,K`.  Substitution
factors the entire remaining coefficient:

```text
H=(K h'+2h K')
  (A M-k(2k-1)L B K^(2k-2)h^(k-1)).                    (68)
```

If `H=1`, the first factor must be a unit.  This is impossible because
`deg h>=2` and

```text
deg(Kh'+2hK')=deg K+deg h-1>=1,                        (69)
```

with nonzero leading coefficient
`(deg h+2 deg K)lc(h)lc(K)`.  The case `T=2` and `R` even is symmetric.
Thus no nonzero constant bracket exists when both outputs have at most two
weights.

Consequently any solution of `(45)` must have support sizes

```text
(2,>=3),             (>=3,2),             or (>=3,>=3), (70)
```

after constants are removed.  In particular it needs at least five nonzero
homogeneous pieces.  This is an all-degree Newton-cell obstruction, not a
bounded coefficient search.

## 7. Connection and loss ledger

| field | exact content |
|---|---|
| source | the middle-coordinate leaf `F2=0`, then its plane filling `A2_(x,q)` |
| target | first `A2_(a,c)`, then the smooth Danielewski surface `Y` |
| map | rational Keller pair `(11)`, polynomial completion `Phi` in `(27)` |
| preserved | determinant on the punctured chart, generic degree three, and the full triple collision |
| destroyed by direct filling | regularity of `a` along `D=0` |
| polynomial repair | retain exactly the intersection ring `(15)` |
| repair cost | `Pic(Y)=Z`, `chi(Y)=2`, and one omitted target point |
| cheapest decisive next test | solve `(45)` in the first genuinely surviving two-by-three weight cells; fixed-weight, one-homogeneous, two-by-two, and generator-slice cells are closed |

The analogy is exact: retaining the puncture keeps an affine-plane target
and a Keller collision; filling the puncture keeps polynomiality and the
collision but transfers the missing coordinate into a target divisor class.

## 8. Exact companion

The companion independently expands `(1)--(17)`, verifies determinant and
both collision forms, checks the smoothness ideal, the open and all three
boundary etale charts, the inverse cubic, and weights `-1` through `-8` of
the intersection normal form with insufficient-cancellation hostiles.  It
also checks the complete affine-linear Poisson cell, the polynomial
representative `(47)` of `db/c`, finite positive/hostile controls for `(55)`,
and the factorization `(68)` for `1<=k<=8`.  The all-weight intersection,
Picard, de Rham, homogeneous nonentry, and complete two-by-two nonentry
statements are structural proofs above, not finite extrapolations from those
rows.

Reproduce with

```bash
python3 04-computation/jc2_rational_keller_danielewski_completion_thm3561.py
python3 -O 04-computation/jc2_rational_keller_danielewski_completion_thm3561.py
```

The two transcripts agree byte-for-byte, and their LF hashes are pinned in
the frontmatter.  Independent hostile audits have checked the structural
intersection, image, Picard, de Rham, and Newton-cell proofs.
