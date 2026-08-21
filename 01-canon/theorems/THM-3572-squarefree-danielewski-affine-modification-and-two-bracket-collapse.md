---
id: THM-3572
title: "Squarefree Danielewski affine modification and two-bracket collapse"
status: >
  PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT AUDIT.
  For every squarefree polynomial Sigma of degree h at least two, the smooth
  surface c^2 e=Sigma(b) has Picard rank h-1 and an everywhere symplectic
  Poisson structure whose physical form db wedge dc/c^2 is globally exact.
  A Bezout identity for Sigma and Sigma' expresses the constant function one
  as a sum of two polynomial Poisson brackets.  Nevertheless no single
  constant bracket can have one homogeneous output, or at most two weights
  in each output; a Darboux pair must use at least five nonzero weight pieces.
  The first genuinely three-arm example supports an explicit degree-five
  etale noninjective polynomial map from A2.  A single Darboux pair on that
  non-A2 target would therefore be a planar Jacobian counterexample.  No such
  pair, and no counterexample to JC(2), is claimed.
source: root, 2026-08-20
audit: PENDING INDEPENDENT HOSTILE AUDIT
depends_on: []
related:
  - THM-3406-affine-modification-power-jets-and-principal-part-transgression
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3566-chebyshev-pell-odd-keller-collision-tower
external:
  - "Dubouloz, Kunyavskii, and Regeta, Bracket Width of Simple Lie Algebras, Documenta Mathematica 26 (2021), 1601-1627."
  - "Dubouloz and Poloni, On a class of Danielewski surfaces in affine 3-space, arXiv:math/0602549."
script: 04-computation/jc2_squarefree_danielewski_two_bracket_thm3572.py
output: 05-knowledge/results/jc2_squarefree_danielewski_two_bracket_thm3572.out
script_sha256: PENDING
output_sha256: PENDING
hash_basis: LF-normalized bytes
---

# THM-3572 -- squarefree Danielewski affine modification and two-bracket collapse

**PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; PENDING INDEPENDENT AUDIT.**
This theorem does **not** refute `JC(2)`.  It isolates a concrete family of
non-`A2` targets on which the topological obstruction has already failed, the
constant has Poisson-bracket length at most two, and the remaining debt is
the compression of those two brackets to one.

All rings and varieties are over `C`.  Fix a squarefree polynomial
`Sigma in C[b]` of degree `h>=2`, with distinct roots
`beta_1,...,beta_h`, and put

```text
A_Sigma=C[b,c,e]/(c^2 e-Sigma(b)),
Y_Sigma=Spec A_Sigma.                                  (1)
```

## 1. Geometry of the exponent-two affine modification

The gradient of the defining equation is

```text
(-Sigma'(b), 2ce, c^2).                                (2)
```

At `c=0`, equation `(1)` forces `b=beta_i`, where squarefreeness gives
`Sigma'(beta_i)!=0`.  Hence `Y_Sigma` is smooth.  Its boundary consists of
the `h` disjoint affine lines

```text
L_i=V(c,b-beta_i) ~= A1_e.                             (3)
```

On `c!=0`, eliminate `e` to obtain

```text
Y_Sigma\V(c) ~= A1_b x Gm_c.                           (4)
```

Nagata localization now gives

```text
Cl(Y_Sigma)=<L_1,...,L_h>/<L_1+...+L_h>
           ~= Z^(h-1).                                (5)
```

Indeed the localized ring is factorial with unit group `C^*c^Z`, and the
only boundary relation supplied by a localized unit is
`div(c)=sum_i L_i`.  A global unit restricts to `lambda c^n`; its boundary
divisor vanishes only for `n=0`, so `A_Sigma^*=C^*`.  Smoothness gives

```text
Pic(Y_Sigma)=Cl(Y_Sigma)=Z^(h-1).                      (6)
```

Additivity of compactly supported Euler characteristic in `(3)--(4)` gives

```text
chi_c(Y_Sigma)=h.                                      (7)
```

Thus `Y_Sigma` is not `A2` whenever `h>=2`.

There is a useful comparison with the standard Danielewski surface

```text
D_Sigma=Spec C[b,c,d]/(cd-Sigma(b)).                   (8)
```

The morphism

```text
pi:Y_Sigma -> D_Sigma,       (b,c,e) |-> (b,c,d=ce),   (9)
```

is an isomorphism on `c!=0` and contracts each line `L_i` to the smooth
point `(beta_i,0,0)`.  It is not finite.  This nonfinite affine modification
is load-bearing below.

## 2. The exact symplectic primitive

Give `A_Sigma` the Poisson bracket

```text
{b,c}=c^2,            {c,e}=-Sigma'(b),
{b,e}=-2ce.                                           (10)
```

The relation in `(1)` is Poisson, so `(10)` is well-defined.  It is
nondegenerate: on `c!=0` this follows from `{b,c}=c^2`, while at `L_i` it
follows from `{c,e}=-Sigma'(beta_i)!=0`.  On `c!=0` the inverse symplectic
form is

```text
omega=db wedge dc/c^2.                                (11)
```

Squarefreeness is equivalent to the existence of `U,T in C[b]` with

```text
U Sigma-T Sigma'=1.                                   (12)
```

For any such pair define the polynomial one-form

```text
alpha=Uce db-2Te dc-Tc de.                            (13)
```

Differentiating `c^2e=Sigma(b)` gives

```text
c^2 de+2ce dc=Sigma' db.                              (14)
```

Equations `(12)--(14)` imply the cleared global identity

```text
c alpha=db.                                           (15)
```

Therefore `alpha=db/c` on `c!=0`, and

```text
d alpha=omega                                         (16)
```

globally.  In particular, the physical symplectic form is exact.

The comparison `(9)` explains exactly what changed.  The standard
Danielewski form is

```text
omega_D=db wedge dc/c,
pi^*omega_D=c omega.                                  (17)
```

The class `[c omega]` remains nonzero; it is the additional division by the
boundary equation `c` that becomes regular and exact after the exponent-two
modification.  Thus the standard topological obstruction was not silently
discarded: it was shifted by one boundary order.

## 3. Complete algebraic de Rham calculation

Give the coordinates weights

```text
wt(b)=0,             wt(c)=1,             wt(e)=-2.   (18)
```

The Euler field `E=c partial_c-2e partial_e` contracts every nonzero-weight
de Rham subcomplex.  In weight zero put

```text
eta=ce dc,           theta=c omega=db wedge dc/c.     (19)
```

Using `(14)`, every weight-zero differential has a unique representative in
the complex

```text
Omega^0_0=C[b],
Omega^1_0={f(b)db+g(b)eta},
Omega^2_0=C[b]theta,

d(f db+g eta)=(Sigma g)' theta.                       (20)
```

It follows that

```text
H^1_dR(Y_Sigma)=0,
H^2_dR(Y_Sigma)=C[b]/{(Sigma g)':g in C[b]},
dim H^2_dR(Y_Sigma)=h-1.                              (21)
```

If `H'=f`, then the `h-1` functionals

```text
f |-> H(beta_i)-H(beta_1),           2<=i<=h,          (22)
```

annihilate every `(Sigma g)'` and form a basis of the dual quotient.  In
particular

```text
[omega]=0,             [theta]=[c omega]!=0.           (23)
```

The Picard and de Rham ranks therefore agree, but the symplectic form relevant
to the planar construction lies in a nonzero weight and is already exact.

## 4. The constant has bracket length at most two

Rewrite `(13)` as

```text
alpha=(U+T')ce db-Te dc-d(Tce).                       (24)
```

Taking exterior derivatives and using `(16)` gives the exact polynomial
identity

```text
1={(U+T')ce,b}+{-Te,c}.                               (25)
```

Thus `1` is always a sum of two Poisson brackets in `A_Sigma`.

For the quadratic surface of THM-3561,

```text
Sigma=b(b+4),       U=-1/4,       T=-(b+2)/8,          (26)
```

and `(25)` becomes

```text
1={-3ce/8,b}+{(b+2)e/8,c}.                            (27)
```

This is a sharper warning than symplectic exactness alone: the residual
Jacobian problem is not to manufacture the constant from brackets, but to
compress a known two-bracket decomposition to one bracket.

Because `H^1_dR(Y_Sigma)=0`, a single Darboux pair `{P,Q}=1` exists if and
only if

```text
alpha=P dQ+d phi                    for some phi in A_Sigma. (28)
```

Equation `(28)` is the exact polynomial-clearing problem.

## 5. Every homogeneous and two-by-two weight cell fails

The weight-`r` part of `A_Sigma` is

```text
(A_Sigma)_r={c^r f(b): f in C[b],
   Sigma^ceil(-r/2) divides f when r<0}.               (29)
```

On `c!=0`, homogeneous brackets satisfy

```text
{c^r f,c^s g}=c^(r+s+1)(s f'g-r f g').               (30)
```

### 5.1 No homogeneous Darboux pair

If homogeneous `P,Q` had `{P,Q}=1`, their weights would be `r,-r-1`.
After swapping, take `r>=0` and write

```text
P=c^r f,
Q=c^(-r-1) Sigma^m g,       m=ceil((r+1)/2).           (31)
```

Putting `H=Sigma^m g`, equation `(30)` gives

```text
{P,Q}=-(r+1)f'H-rfH'.                                 (32)
```

For `r=0`, `(32)` is divisible by `Sigma`; for `r>=2`, it is divisible by
`Sigma^(m-1)`.  For `r=1`, its leading term has positive degree and nonzero
coefficient

```text
-(2 deg f+deg g+deg Sigma)lc(f)lc(g)lc(Sigma).         (33)
```

Hence `(32)` is never a nonzero constant.  If only one of two arbitrary
outputs were homogeneous, weight decomposition of the other would isolate
the same forbidden complementary component.  Every Darboux pair must mix at
least two weights in both outputs.

### 5.2 The complete two-by-two cell is empty

Suppose, after subtracting constants, that each of `P,Q` has exactly two
nonzero weight pieces and `{P,Q}=1`.  A single complementary pair would be a
forbidden homogeneous pair, so ordering forces cross matching.  The unique
extreme brackets must vanish.  Root valuations at any `beta_i` show that a
vanishing extreme cannot mix a negative and a positive weight: if the
coefficient orders are `A>0,B>=0`, its leading Wronskian coefficient is
`qA-pB!=0` for `p<0<q`.  A zero-weight boundary would force the retained
zero-weight coefficient to be scalar, already removed.

Consequently there are `R,T>=1` and nonzero `f,g,F,G in C[b]` such that

```text
P=c^(-R)f+c^(T-1)F,
Q=c^(-T)g+c^(R-1)G,                                   (34)
```

with `Sigma^ceil(R/2)|f` and `Sigma^ceil(T/2)|g`.  The two extreme equations
and the weight-zero coefficient are

```text
Rfg'-Tf'g=0,
(R-1)F'G-(T-1)FG'=0,                                  (35)

H=(R-1)f'G+RfG'-TF'g-(T-1)Fg'=1.                     (36)
```

At a root of `Sigma`, only the first or last term of `(36)` can survive.
The divisibility in `(29)` and the boundary coefficients show that this
requires `R=2` or `T=2`.

Assume `R=2`.  The first equation in `(35)` and `H(beta_i)=1` force
`T=2k` even and make `f` simple at every root.  The two zero Wronskians give

```text
f=A h,            g=B h^k,             Sigma|h,
F=L K^(2k-1),     G=M K,                               (37)
```

for nonzero constants `A,B,L,M` and nonzero polynomials `h,K`.  Substitution
factors `(36)` as

```text
H=(K h'+2h K')
  (AM-k(2k-1)LB K^(2k-2)h^(k-1)).                     (38)
```

If `H=1`, the first factor must be a unit.  This is impossible because
`Sigma|h`, `deg Sigma>=2`, and

```text
deg(Kh'+2hK')=deg K+deg h-1>=1,                       (39)
```

with nonzero leading coefficient
`(deg h+2deg K)lc(h)lc(K)`.  The case `T=2` is symmetric.

Therefore no constant bracket exists when both outputs use at most two
weights.  Any solution has support sizes

```text
(2,>=3),           (>=3,2),           or (>=3,>=3),   (40)
```

and at least five nonzero homogeneous pieces after constants are removed.

## 6. A degree-five, three-arm near-counterexample

The general system is not merely formal.  Put

```text
t=x^2q,              P=1+t^2,              S=1+5t^2,
a=q/S^2,             c=xSP,                 b=tP^2,
W=125t^6+450t^4+565t^2+256,
e=qW.                                                   (41)
```

The rational pair has

```text
Jac_(x,q)(a,c)=-1.                                    (42)
```

Moreover

```text
(3125b^2+256)=S^2W,
c^2e=b(3125b^2+256).                                  (43)
```

Thus `(b,c,e)` defines a polynomial morphism

```text
Phi_5:A2_(x,q) -> Y_3,
Y_3: c^2e=Sigma_3(b),       Sigma_3=b(3125b^2+256).    (44)
```

On the dense set `S!=0`, `b=ac^2` and
`e=a(3125b^2+256)`, so `(42)` implies

```text
Jac(P(Phi_5),Q(Phi_5))=-{P,Q}(Phi_5)                  (45)
```

for all `P,Q in C[Y_3]`.  Both sides are polynomial, hence `(45)` holds
everywhere.  Since the brackets on source and target are nondegenerate,
`Phi_5` is etale.

The equation `b=t(1+t^2)^2` has generic degree five; after `t` is chosen,
`x` and `q` are recovered from `c=xSP` and `t=x^2q`.  Hence `Phi_5` has
generic degree five.  It is visibly noninjective.  For every `lambda!=0`,
the target point

```text
(b,c,e)=(0,0,256lambda)                               (46)
```

has five distinct preimages: `(x,q)=(0,lambda)`, and for each
`t in {i,-i}` the two choices

```text
q=16lambda,              x^2=t/(16lambda).             (47)
```

Therefore a pair `P,Q in C[Y_3]` with `{P,Q}=1` would make
`(P(Phi_5),Q(Phi_5))` a polynomial, etale, noninjective endomorphism of
`A2`: a genuine counterexample to `JC(2)`.  Sections 4--5 show simultaneously
that `1` is already a sum of two brackets and that no one-bracket solution
can lie in the complete two-by-two weight cell.

For this symmetric cubic,

```text
U=-28125b/131072,
T=-(9375b^2+512)/131072,                               (48)
```

satisfy `(12)`, so `(25)` is completely explicit.  The involution
`(b,c,e)|->(-b,c,-e)` is anti-Poisson; it is a useful parity sidecar, not a
solution of the compression problem.

## 7. What the standard bracket-width analogy does and does not give

The bracket-width literature studies standard Danielewski surfaces such as
`D_Sigma` in `(8)`.  The comparison is exact only through the nonfinite map
`pi` in `(9)`:

```text
source object:       Y_Sigma with exponent-two boundary,
target object:       standard D_Sigma,
map:                 d=ce,
preserved:           the open chart c!=0 and the root labels beta_i,
destroyed:           each entire boundary line L_i,
needed sidecar:      boundary order / the extra factor c,
cheapest hostile:    pi^*omega_D=c omega is nonexact while omega is exact.
```

Thus a standard-Danielewski obstruction cannot be transferred by syntax.
Conversely, the two-bracket identity `(25)` is not a Darboux pair: the
missing datum is a common polynomial second coordinate that compresses both
summands.

## 8. Scope and next target

The theorem proves:

1. a universal exact primitive and two-bracket decomposition for every
   squarefree `Sigma` of degree at least two;
2. the complete homogeneous and two-by-two Newton cells are empty;
3. an explicit degree-five etale noninjective map to a three-arm non-`A2`
   target.

It does **not** prove a one-bracket pair exists, that the target is `A2`, or
that `JC(2)` fails.  The first surviving construction cells have weight
support `(2,3)` and `(3,2)`.  The sharp next calculation is the symmetric
cubic `Sigma_3`, where the anti-Poisson involution can organize those five
pieces without erasing either of the two independent boundary classes.

## 9. Exact verification contract

The companion uses exact rational symbolic arithmetic.  It checks:

- smooth squarefree controls and Bezout primitives for several degrees;
- the cleared identity `c alpha=db`, `d alpha=omega`, and `(25)`;
- de Rham quotient ranks and the root-antiderivative functionals;
- the universal homogeneous formula and bounded structural controls for the
  two-by-two factorization, including `1<=k<=8`;
- the quadratic identity `(27)`;
- every polynomial identity in the degree-five packet `(41)--(48)`, its
  rational Jacobian, five collision preimages, and anti-Poisson involution.

The finite rows are hostile and positive controls for the all-degree proofs;
they are not extrapolated universes.
