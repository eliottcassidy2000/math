---
id: THM-3572
title: "Squarefree Danielewski affine modification and two-bracket collapse"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every squarefree polynomial Sigma of degree h at least two, the smooth
  surface c^2 e=Sigma(b) has Picard rank h-1 and an everywhere symplectic
  Poisson structure whose physical form db wedge dc/c^2 is globally exact.
  A Bezout identity for Sigma and Sigma' expresses the constant function one
  as a sum of two polynomial Poisson brackets.  Nevertheless no single
  constant bracket can have one homogeneous output, at most two weights in
  each output, or a two-by-three support; a Darboux pair must use at least six
  nonzero weight pieces.
  The first genuinely three-arm example supports an explicit degree-five
  etale noninjective polynomial map from A2.  A single Darboux pair on that
  non-A2 target would therefore be a planar Jacobian counterexample.  No such
  pair, and no counterexample to JC(2), is claimed.
source: root, 2026-08-20
audit: >
  PASS.  The squarefree-arm geometry, de Rham and Poisson identities,
  homogeneous/two-by-two/two-by-three support obstructions, all-width
  {-1,1} owner obstruction, selected three-by-three obstruction, and the
  degree-five S5 carrier were independently rederived.  Normal, optimized,
  and stored companion transcripts are byte-identical.
depends_on:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
related:
  - THM-1350-equivariant-fixed-locus-JC
  - THM-3406-affine-modification-power-jets-and-principal-part-transgression
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3566-chebyshev-pell-odd-keller-collision-tower
external:
  - "Dubouloz, Kunyavskii, and Regeta, Bracket Width of Simple Lie Algebras, Documenta Mathematica 26 (2021), 1601-1627."
  - "Dubouloz and Poloni, On a class of Danielewski surfaces in affine 3-space, arXiv:math/0602549."
script: 04-computation/jc2_squarefree_danielewski_two_bracket_thm3572.py
output: 05-knowledge/results/jc2_squarefree_danielewski_two_bracket_thm3572.out
script_sha256: a1ce52375c7264400e7f2b6b9970aa75b3ca16ec34f6611686889ed84934c502
output_sha256: d7bb058966dc1534c000ea4e57102156493b1c4664d4eb1b1cf32407385aa52d
hash_basis: LF-normalized bytes
---

# THM-3572 -- squarefree Danielewski affine modification and two-bracket collapse

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
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

There is also an exact Poisson-homology reading.  Multiplication by `omega`
identifies the linear span of polynomial brackets with the exact polynomial
two-forms, because `{P,Q}omega=dP wedge dQ`.  Hence

```text
A_Sigma/span{brackets}
  ~= c C[b] / c{(Sigma g)':g in C[b]}.                (23a)
```

The constant `1` vanishes in this quotient, while the `h-1` arm classes are
genuine.  For the cubic of Section 6 they have representatives `c,bc`; the
two-bracket identity below is therefore not hiding a nonzero commutator
class.

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
weights.  The universal squarefree theorem THM-3569 closes the entire
two-by-three cell as well.  Any solution has support sizes

```text
(2,>=4),           (>=4,2),           or (>=3,>=3),   (40)
```

and at least six nonzero homogeneous pieces after constants are removed.

### 5.3 The canonical two-bracket compression owner fails at every width

A natural two-sided compression ansatz has one output supported at weights
`{-1,1}`.  This entire lane is empty, even if the other output has
arbitrarily many weights.

Indeed write

```text
P=c^(-1)f+cF,                    Sigma|f,              (41)
```

and decompose an arbitrary `Q` into parity chains.  Only the even chain
containing weights `-2,0` can contribute to output weight zero; all odd or
disconnected chains have disjoint nonzero output supports and must vanish
separately before being deleted.  A chain containing only one of weights
`-2,0` reduces its scalar row to the forbidden homogeneous cell.  Write the
unique remaining finite chain as

```text
Q_even=sum_(j=L)^U c^(2j)q_j,       L<=-1<=0<=U.       (42)
```

The coefficient of output weight `2j` is

```text
L_j=fq_j'+2j f'q_j+2(j-1)F'q_(j-1)-Fq_(j-1)'.        (43)
```

At the lower endpoint, `L_L=0` integrates to
`q_L=Cf^(-2L)`.  Put `z=fF`.  Ascending through the zero layers only as far
as `j=-1` gives

```text
q_j=f^(-2j)R_j(z),                  L<=j<=-1.          (44)
```

To see the closure, multiply `(43)` by `f^(2j-1)`: its first half becomes
the derivative of `f^(2j)q_j`, while the second becomes

```text
z'[2(j-1)R_(j-1)-zR_(j-1)'].                         (45)
```

Independently, the upper endpoint gives `q_U=CF^(2U)`.  Descending through
the zero layers `j=U,...,1` gives

```text
q_j=F^(2j)S_j(z),                       0<=j<=U.
```

Here the two halves of `(43)`, after the appropriate power of `F` is
removed, are respectively
`z'[2jS_j+zS_j']` and `-z'S_(j-1)'`.  Thus
`q_(-1)=f^2R_(-1)(z)` and `q_0=S_0(z)`.  Only now is the scalar layer taken:

```text
L_0=f z'[S_0'-2R_(-1)-zR_(-1)'].                     (46)
```

It is divisible by the nonconstant polynomial `f`, so it cannot equal one.
Thus no Darboux pair can have either output supported exactly at weights
`{-1,1}`.  This closes the most direct attempt to compress `(25)`, not merely
a bounded truncation of that attempt.

### 5.4 Locally nilpotent Hamiltonian shears cannot work

The cited Makar--Limanov classification for surfaces
`x^n z-P(y)`, `n>=2`, `deg P>=2`, applies to `(1)` and gives

```text
ML(Y_Sigma)=C[c].                                      (47)
```

Equivalently, every algebraic `Ga`-action preserves the unique reducible
`A1`-fibre `c=0`.  This has a direct Darboux consequence.  If a Hamiltonian
derivation `D_P={P,-}` were locally nilpotent, `(47)` would give
`D_P(c)=0`.  The kernel of `D_c={c,-}` is `C[c]`, hence `P in C[c]`.
Then `{P,Q}=P'(c){c,Q}=1` would force `P'` to be a nonzero constant and
would give a polynomial slice for `D_c`.  A slice would identify
`A_Sigma=C[c,Q]`, contradicting

```text
A_Sigma/(c)=C[b,e]/(Sigma(b)),                         (48)
```

the disjoint union of `h>=2` affine lines.  Therefore both Hamiltonian
derivations in any Darboux pair must be genuinely non-locally-nilpotent.
Tame `Ga` flows cannot clear the rational pole `1/c`.

### 5.5 A natural parity-compatible three-by-three cell also fails

For the symmetric cubic in Section 6, the anti-Poisson involution suggests
the balanced supports

```text
P=c^(-3)f+c^(-1)a+cF,
Q=c^(-2)g+q_0+c^2H.                                  (48a)
```

This entire support pattern is empty for every squarefree `Sigma`.  The two
extreme equations give

```text
f=A k^3,              g=B k^2,              H=C F^2, (48b)
```

with nonzero constants `A,B,C`, a nonconstant `k`, and `Sigma|k`.  The next
two weight equations integrate to

```text
(a/k)'=(3A/(2B))q_0',          q_0'=2C(aF)'.          (48c)
```

Consequently, for constants `lambda!=0,E`,

```text
a(1-lambda kF)=E k.                                   (48d)
```

If `E=0`, then `kF` is a nonzero constant.  If `E!=0`, coprimality of `k`
and `1-lambda kF` forces the latter factor to be a unit, again making `kF`
constant.  Both conclusions contradict the nonconstancy and divisibility of
`k`.  This closes the most symmetric six-piece compression, but it is not a
classification of all three-by-three supports.

## 6. A degree-five, three-arm near-counterexample

The general system is not merely formal.  Put

```text
t=x^2q,              P=1+t^2,              S=1+5t^2,
a=q/S^2,             c=xSP,                 b=tP^2,
W=125t^6+450t^4+565t^2+256,
e=qW.                                                   (49)
```

The rational pair has

```text
Jac_(x,q)(a,c)=-1.                                    (50)
```

Moreover

```text
(3125b^2+256)=S^2W,
c^2e=b(3125b^2+256).                                  (51)
```

This packet is governed by the one-variable critical-value map

```text
B(t)=t(1+t^2)^2,             B'(t)=(1+t^2)(1+5t^2),
Disc_t(B(t)-b)=b^2(3125b^2+256).                      (52)
```

Thus `Sigma_3` is exactly the reduced finite critical-value polynomial.  The
two roots of `3125b^2+256` are the two side critical values; the central
critical value `b=0` has the two collision roots `t=+/-i` and the unramified
root `t=0`.  This derivative/dessin sidecar explains both the etale
completion and its three boundary arms.

The same passport distinguishes this example sharply from the dihedral
Chebyshev tower.  Infinity supplies a 5-cycle, while either side critical
value has exactly one simple critical point and therefore transposition
inertia.  Conjugating that transposition by the 5-cycle gives a connected
transposition graph, so the generic Galois closure of `B(t)-b` has group
`S_5`.  Adjoining the independent coordinate `c` does not change that group.
Thus `Phi_5` is a genuinely nonradical, nonsolvable near-counterexample
carrier rather than a hidden Kummer cover.

Thus `(b,c,e)` defines a polynomial morphism

```text
Phi_5:A2_(x,q) -> Y_3,
Y_3: c^2e=Sigma_3(b),       Sigma_3=b(3125b^2+256).    (53)
```

On the dense set `S!=0`, `b=ac^2` and
`e=a(3125b^2+256)`, so `(50)` implies

```text
Jac(P(Phi_5),Q(Phi_5))=-{P,Q}(Phi_5)                  (54)
```

for all `P,Q in C[Y_3]`.  Both sides are polynomial, hence `(54)` holds
everywhere.  Since the brackets on source and target are nondegenerate,
`Phi_5` is etale.

The rational map `B:P1_t->P1_b` has degree five, so
`[C(t):C(B)]=5`; equivalently `t(1+t^2)^2-b` is irreducible over `C(b)`.
It remains irreducible over the purely transcendental extension `C(b,c)`.
After `t` is chosen, `x` and `q` are recovered from `c=xSP` and `t=x^2q`.
Hence `Phi_5` has generic degree five.  It is visibly noninjective.  For
every `lambda!=0`,
the target point

```text
(b,c,e)=(0,0,256lambda)                               (55)
```

has five distinct preimages: `(x,q)=(0,lambda)`, and for each
`t in {i,-i}` the two choices

```text
q=16lambda,              x^2=t/(16lambda).             (56)
```

The geometric image is exact:

```text
Phi_5(A2(C))=Y_3(C) minus
  {(beta,0,0):3125beta^2+256=0}.                       (57)
```

Indeed, for `c!=0` a root of `B(t)=b` with `PS!=0` recovers
`x=c/(PS)` and `q=t/x^2`.  Such a root is automatic away from the three
critical values.  At `b=0`, take `t=0`.  At a side value `beta=B(u)`,
`u^2=-1/5`, one has

```text
B(t)-beta=(t-u)^2 Q_u(t),
Q_u=t^3+2ut^2+(7/5)t+(16/5)u,
Q_u(u)=4u,               Q_u(-u)=8u/5,
Res_t(Q_u,1+t^2)=-16/125.                              (57a)
```

Thus the three roots of `Q_u` avoid `PS`.  On the central arm `b=c=0`,
the unramified root `t=0` supplies every value of `e`.  On either side arm,
the root `u` of `S` supplies every `e!=0` because `W=160` modulo `S`, whereas
`e=0` would force `q=0` and contradict `u!=0`.  The two missing points in
`(57)` are codimension two; the scheme-theoretic closure of the image is
still all of `Y_3`.

The complete nonproperness locus is nevertheless divisorial.  With
`D=PS`, elimination gives

```text
Res_t(W,B-b)=(3125b^2+256)^3,
Res_t(D,B-b)=b^2(3125b^2+256),
gcd(W,W')=gcd(W,tD)=1.                                (57b)
```

Consequently the Jelonek set is the union of five affine lines:

```text
N(Phi_5)=
 union_(Sigma_3(beta)=0){b=beta,e=0}
 union_(3125beta^2+256=0){b=beta,c=0}.                (57c)
```

In the coordinates `t=x^2q`, bounded image along an escaping source sequence
forces either `x->infinity,D(t)->0,e->0` or
`x->0,W(t)->0,c->0`; the simple roots in `(57b)` realize every point on the
listed lines.  Thus etaleness has removed ramification, not the nonproper
sheet debt.

The fibre invoice is correspondingly nonconstant.  On `c!=0` the fibre has
respectively `5,1,3` points over generic, central-critical, and side-critical
values of `b`.  The central arm has five points for `e!=0` and one at its
origin; either side arm has two points for `e!=0` and none at its omitted
origin.  No branch ramifies--the missing points have escaped to infinity.

Therefore a pair `P,Q in C[Y_3]` with `{P,Q}=1` would make
`(P(Phi_5),Q(Phi_5))` a polynomial, etale, noninjective endomorphism of
`A2`: a genuine counterexample to `JC(2)`.  Sections 4--5 show simultaneously
that `1` is already a sum of two brackets and that no one-bracket solution
can lie in the complete two-by-two weight cell.

For this symmetric cubic,

```text
U=-28125b/131072,
T=-(9375b^2+512)/131072,                               (58)
```

satisfy `(12)`, so `(25)` is completely explicit.

None of the generators `b,c,e` can be one coordinate of a Darboux pair.
The image of `{b,-}` lies in `(c)`.  The derivation `{c,-}` is locally
nilpotent with kernel `C[c]` but has no slice because `c=0` is reducible.
On the slice `e=1`, the equation `c^2=Sigma_3(b)` is an affine elliptic
curve, and `{e,-}` is dual to the holomorphic differential `db/(2c)`.
A polynomial slice would make that nonzero holomorphic differential exact on
the projective elliptic completion, which is impossible.

The involution `iota(b,c,e)=(-b,c,-e)` is anti-Poisson.  Its fixed locus is
the line `b=e=0`.  For an `iota`-eigen Darboux pair the parities must be
opposite; the odd coordinate vanishes on that line, while etaleness forces
the even coordinate to restrict affinely with nonzero `c` coefficient.  This
is the two-dimensional fixed-locus gate of the cited equivariant theorem,
and it is a useful hostile test rather than a construction.

## 7. What the standard bracket-width analogy does and does not give

The cited Danielewski bracket-width paper proves width at most two for the
standard quadratic surface and relates width one to a global one-form
decomposition.  Its width concerns the entire Lie algebra generated by
locally nilpotent vector fields, not merely the bracket length of the
constant `1`.  On the standard surface `[omega_D]!=0`, so `1` is not even in
the derived Poisson space.  Here the nonfinite modification pulls
`omega_D` back to the still-nonexact form `c omega`, while the physical form
`omega` becomes exact and `1` acquires the two-bracket decomposition `(25)`.

The comparison is therefore exact only through the nonfinite map `pi` in
`(9)`:

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
2. the complete homogeneous, two-by-two, and two-by-three Newton cells are
   empty;
3. the canonical `{-1,1}` compression owner fails at every width, locally
   nilpotent Hamiltonian shears fail, and a natural symmetric three-by-three
   cell fails;
4. an explicit degree-five etale noninjective map onto the complement of two
   points in a three-arm non-`A2` target.

It does **not** prove a one-bracket pair exists, that the target is `A2`, or
that `JC(2)` fails.  The first formally surviving support sizes are
`(2,4)`, `(4,2)`, and `(3,3)`; Section 5.5 kills one particularly symmetric
three-by-three support but does not classify any complete surviving family.
Any
eventual Darboux map `Y_3->A2` must itself be nonproper and have generic
degree at least two: a finite etale cover of `A2_C` is trivial, while a
degree-one etale map, by Zariski Main, would identify this normal surface with
an affine open of `A2`.  A divisorial complement would create a nonconstant
unit, contrary to `A_Sigma^*=C^*`; a codimension-two complement has the same
normal coordinate ring as `A2`, so affineness would force the whole plane,
contrary to `(6)`.

## 9. Exact verification contract

The companion uses exact rational symbolic arithmetic.  It checks:

- smooth squarefree controls and Bezout primitives for several degrees;
- the cleared identity `c alpha=db`, `d alpha=omega`, and `(25)`;
- de Rham quotient ranks and the root-antiderivative functionals;
- the universal homogeneous formula and bounded structural controls for the
  two-by-two factorization, including `1<=k<=8`;
- bounded exact controls for the all-width `{-1,1}` recursion and the
  symmetric three-by-three obstruction;
- the quadratic identity `(27)`;
- every polynomial identity in the degree-five packet `(49)--(58)`, its
  critical-value discriminant, rational Jacobian, five collision preimages,
  exact two-point image complement, and anti-Poisson involution.

The finite rows are hostile and positive controls for the all-degree proofs;
they are not extrapolated universes.
