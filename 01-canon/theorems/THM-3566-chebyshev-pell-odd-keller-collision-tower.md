---
id: THM-3566
title: "Chebyshev--Pell odd Keller collision tower and Danielewski completion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every odd
  d=2r+1>=3 there is an explicit equivariant rational plane Keller pair of
  determinant -1 whose degree-d Chebyshev collision cover has passport
  (1,2^r) over each of two values.  Its maximal polynomial-target-observable
  intersection is the smooth generalized Danielewski ring
  C[b,c,e]/(c^2e-b(b-d^-2)).  The induced polynomial map A2 to that surface
  is everywhere etale, has generic degree d, and its geometric image over C
  is the complement of one point.  Its scheme-theoretic image closure is the
  whole target.  This is a non-A2 near-counterexample tower, not a planar
  Jacobian counterexample.
source: root/rational_keller_family, 2026-08-20
audit: >
  An independent hostile audit rederived the master Jacobian identity,
  Chebyshev--Pell construction, maximal intersection, smoothness, etaleness,
  degree, geometric image, fibre passports, three-component nonproperness
  locus, dihedral carrier, residual ramification budget, symplectic primitive,
  two-bracket identity, and fixed-power parity classification.  It added
  exact controls at odd degrees 17,19,21,25 and fixed powers 13 through 16.
  Ordinary and optimized replays are byte-identical; active-gate, no-assert,
  documentation, hash, and diff checks pass.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
related:
  - THM-3404-factorized-danielewski-principal-parts-and-finite-cover-obstruction
  - THM-3554-punctured-kummer-collision-surface-normal-form
script: 04-computation/jc2_chebyshev_pell_odd_keller_collision_tower_thm3566.py
output: 05-knowledge/results/jc2_chebyshev_pell_odd_keller_collision_tower_thm3566.out
script_sha256: 1f7b8cb7e5748ba5c8e2c50d0daef7d5c191fcc01b2f28c302823ef844c24b91
output_sha256: c297317869cf183b9a41a006ee1f1cc640abdedcb065a08b0ed3824aaf63952f
hash_basis: LF-normalized bytes
---

# THM-3566 -- Chebyshev--Pell odd Keller collision tower

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem
does **not** disprove `JC(2)`, does not claim that its formulas are new in the
literature, and does not turn a rational Keller pair into a polynomial
endomorphism of `A2`.  It gives an all-odd tower of exact near-counterexamples
and identifies, without a degree cutoff, the polynomial target information
that survives their poles.

All rings and varieties are over `C`.  Write

```text
Jac(f,g)=f_x g_q-f_q g_x.                              (1)
```

## 1. The equivariant master identity

Put `t=x^2q` and consider the weight-equivariant rational pair

```text
a=q H(t),                 c=x G(t).                    (2)
```

The weight-zero observable

```text
B(t)=t H(t)G(t)^2=ac^2                                 (3)
```

contains the whole Jacobian equation.  Direct differentiation, with no
division until the last displayed equality, gives

```text
Jac(a,c)=-H G-2tH G'-tH'G=-B'(t)/G(t).                (4)
```

Now choose a polynomial `P(t)` with `P(0)!=0`, put

```text
B=tP^2,                    S=P+2tP',                  (5)
a=q/S^2,                   c=xSP=xB'.                 (6)
```

Then `B'=PS`, so `(3)` holds with `H=S^-2,G=PS`, and

```text
Jac(a,c)=-1.                                           (7)
```

Thus every choice of `P` produces a rational Keller pair on `S!=0`.  The
hard part is not the constant Jacobian; it is arranging the remaining
critical values so that the poles admit a small polynomial completion.

## 2. The odd Chebyshev--Pell specialization

Fix

```text
d=2r+1 >= 3,             beta=d^-2,                  (8)
z=1-2t.
```

For the Chebyshev polynomials `T_j,U_j` of the first and second kind, set

```text
P_d(t)=(U_r(z)+U_(r-1)(z))/d,
S_d(t)= U_r(z)-U_(r-1)(z).                             (9)
```

The two half-angle identities are

```text
B_d=tP_d^2=(1-T_d(z))/(2d^2),
S_d=P_d+2tP_d',
B_d-beta=(t-1)S_d^2/d^2,                              (10)
```

and consequently

```text
B_d'=P_dS_d.                                           (11)
```

For example, if `z=cos(theta)` then

```text
P_d=sin(d theta/2)/(d sin(theta/2)),
S_d=cos(d theta/2)/cos(theta/2),                       (12)
```

which proves `(10)` on a Zariski-dense set and hence polynomially.  In
particular

```text
P_d(0)=S_d(0)=1,
P_d(1)=(-1)^r/d,       S_d(1)=(-1)^r d.               (13)
```

Both `P_d` and `S_d` have degree `r`, are squarefree, and are coprime.  Their
roots are respectively

```text
tau_k   = sin^2(k*pi/d),                 1<=k<=r,
sigma_k = sin^2((2k-1)*pi/(2d)),         1<=k<=r.      (14)
```

They are nonzero, are different from `1`, and the two sets are disjoint.
Equations `(10)` therefore give the two exact degree-`d` passports

```text
B_d^-1(0):     t=0 simple,  tau_1,...,tau_r double;
B_d^-1(beta):  t=1 simple,  sigma_1,...,sigma_r double. (15)
```

The specialized rational Keller pair is

```text
a=q/S_d(t)^2,             c=xP_d(t)S_d(t),             (16)
Jac(a,c)=-1.
```

For every `s!=0`, the target point `(b,c,e)=(0,0,s)` constructed below has
`d` distinct inverse images, all lying in `S_d!=0`.  On those points the
rational Keller pair has the common value

```text
(a,c)=(-d^2s,0).                                      (16a)
```

This is the degree-`d` rational Keller collision tower.

## 3. Maximal polynomial target observables

Define the three polynomial pullbacks

```text
b=ac^2=B_d(t)=tP_d(t)^2,
e=a(b-beta)=q(t-1)/d^2,
c=xP_d(t)S_d(t).                                      (17)
```

They satisfy

```text
c^2e=b(b-beta).                                       (18)
```

The exact maximal intersection is

```text
C[a,c] intersect C[x,q] = C[b,c,e].                   (19)
```

Here `C[a,c]` means the **polynomial target ring** of the rational pair
`(16)`.  Equation `(19)` is not a statement about all rational functions in
`a,c`.

### Proof of maximality

Give `x,q` weights `1,-2`.  Then `t,b` have weight zero and `a,c,e` have
weights `-2,1,-2`.  Since `Jac(a,c)=-1`, the functions `a,c` are
algebraically independent.  The weight-`w>=0` part of `C[a,c]` is

```text
c^w C[b],                                              (20)
```

which already lies in `C[b,c]`.  A weight `-s<0` element has the unique form

```text
c^(-s)F(b),                                            (21)
```

and polynomiality in `a,c` forces

```text
b^m divides F(b),                 m=ceil(s/2).         (22)
```

Indeed, a monomial of weight `-s` is
`a^i c^(2i-s)=b^i c^(-s)`, and its `c` exponent is nonnegative exactly when
`i>=m`.

Write `F=b^mF_0`.  On substituting `(17)`, the factor `b^m` cancels all
negative powers of `x` and all possible `P_d` denominators.  The only
remaining possible denominator is `S_d^s`.  At every root of `S_d`, all
other displayed factors are generically units, while `(10)` gives

```text
ord_(S_d=0)(B_d-beta)=2.                              (23)
```

Because the roots of `S_d` are simple and all have the same `B_d`-value,
the valuation inequalities at those `r` divisors are equivalent to

```text
(b-beta)^m divides F_0(b).                             (24)
```

This is the missing boundary coordinate.  Combining `(22)` and `(24)` gives

```text
c^(-s)[b(b-beta)]^m F_1(b)
  = e^m F_1(b),                 s=2m,
  = c e^m F_1(b),               s=2m-1.               (25)
```

Thus every polynomial pullback belongs to `C[b,c,e]`; the reverse inclusion
is visible in `(17)`.  Finally, `b,c` are algebraically independent because

```text
Jac(b,c)=-c^2                                          (26)
```

is generically nonzero.  Hence `(18)` is the full defining relation and
proves `(19)`.

## 4. Smooth polynomial completion and boundary charts

Let

```text
Y_beta=Spec C[b,c,e]/(c^2e-b(b-beta)).                 (27)
```

The gradient of `c^2e-b(b-beta)` is

```text
(beta-2b, 2ce, c^2).                                  (28)
```

Its vanishing would force `c=0,b=beta/2`, contrary to the defining
equation.  Thus `Y_beta` is smooth.  Equations `(17)` define a polynomial
morphism

```text
Phi_d:A2_(x,q) -> Y_beta.                              (29)
```

On `c!=0`, `(b,c)` are target coordinates and `(26)` proves etaleness.  On
`c=0`, the relation forces `b=0` or `b=beta`, and `(c,e)` are target
coordinates.  Directly from `(17)`, or from `(7)` and the Leibniz rule,

```text
Jac(c,e)=2b-beta.                                      (30)
```

It equals `-beta` over `b=0` and `beta` over `b=beta`.  These two charts
cover the boundary, so `Phi_d` is everywhere etale.

For generic `b`, the equation `B_d(t)=b` has `d` roots.  Once a root `t` is
chosen and `c!=0`, the source point is uniquely reconstructed by

```text
x=c/B_d'(t),              q=t B_d'(t)^2/c^2.          (31)
```

Therefore `Phi_d` is quasi-finite of generic degree `d`.

## 5. Exact geometric image and fibre passports

Put

```text
L_0     = {b=0,c=0},
L_beta  = {b=beta,c=0},
p_beta  = (beta,0,0).                                 (32)
```

For a point `(0,0,s)` on `L_0`, one source point is always

```text
(x,q)=(0,-d^2s).                                      (33)
```

If `s!=0`, every root `tau` of `P_d` supplies two further points through

```text
q=d^2s/(tau-1),          x^2=tau(tau-1)/(d^2s).       (34)
```

Thus the fibre has `1+2r=d` points for `s!=0`, while the fibre over
`(0,0,0)` is the single point `(0,0)`.

For `(beta,0,s)` with `s!=0`, every root `sigma` of `S_d` supplies two
points through

```text
q=d^2s/(sigma-1),        x^2=sigma(sigma-1)/(d^2s).   (35)
```

The fibre has `2r=d-1` points.  When `s=0`, `(35)` is impossible and the
simple `t=1` sheet would require simultaneously `x=0` and `x^2q=1`.
Consequently the geometric/topological image on complex points is

```text
Phi_d(A2(C))=Y_beta(C) minus {p_beta}.                 (36)
```

The scheme-theoretic image closure is all of `Y_beta`, because `Phi_d` is
dominant.  Equation `(36)` is not a claim that the omitted codimension-two
point changes the target coordinate ring.

For completeness, every geometric fibre has the following cardinality;
all listed finite points have etale multiplicity one:

| target stratum | fibre size |
|---|---:|
| `c!=0`, `b` not in `{0,beta}` | `d` |
| `c!=0`, `b=0` or `b=beta` | `1` |
| `L_0`, `e!=0` | `d` |
| `(0,0,0)` | `1` |
| `L_beta`, `e!=0` | `d-1` |
| `p_beta` | `0` |

The exact nonproperness (Jelonek) locus is the **three-component** set

```text
C_0     ={b=0,e=0},
C_beta  ={b=beta,e=0},
L_beta  ={b=beta,c=0}.                                (37)
```

This corrects the tempting but false guess that only the omitted-point arm
`L_beta` is nonproper.  Near a root of `P_d`, taking `x` inversely
proportional to `B_d'(t)` produces `C_0`; roots of `S_d` similarly produce
`C_beta`.  Finally, letting `t->1` while keeping
`q(t-1)/d^2=s` produces all of `L_beta`.  Conversely, outside `(37)` the
table gives exactly `d` etale points.  More formally, Zariski's Main Theorem
factors the quasi-finite map through an open immersion of `A2` into the
normalization of `Y_beta` in `C(x,q)`, followed by a finite map of degree
`d`.  The target is smooth and the normalization is Cohen--Macaulay, so
miracle flatness makes this finite map locally free of rank `d`.  A target
fibre already containing `d` distinct etale source points has no degree left
for a normalization-boundary point.  This proves equality in `(37)`, not
merely containment.

## 6. Dihedral Kummer carrier and the ramification budget

The degree-`d` equation is

```text
1-2d^2b=T_d(z),                  z=1-2t.              (38)
```

On the Kummer normalization

```text
u=z+sqrt(z^2-1),                 z=(u+u^-1)/2,         (39)
```

it becomes

```text
T_d(z)=(u^d+u^-d)/2.                                  (40)
```

Thus the Galois closure of `C(t)/C(b)` has dihedral monodromy `D_d`: the
Kummer rotations `u->zeta_d u` and the reflection `u->u^-1` generate it.
The two passports `(15)` are the chain dessin of this quotient.

There is also a useful obstruction that does not require the Chebyshev
choice.  Suppose `P` is squarefree of degree `r`, `P(0)!=0`, and

```text
B=tP^2,                    deg B=2r+1=d.              (41)
```

Then

```text
B'=P(P+2tP')=PS.                                       (42)
```

The `r` roots of `P` already spend `r` units of the finite
Riemann--Hurwitz ramification budget on the prescribed double collision
roots over `0`.  Since `deg B'=2r`, exactly `r` residual critical roots,
counted with multiplicity, remain in `S`.  They cannot be roots of `P` or
`t`, and their `B`-values are nonzero.  Their exceptional boundary arms can
coalesce -- `(10)` makes all of them coalesce at `beta` -- but they cannot
disappear.

This is a scoped obstruction to a **residual-arm-free** completion: no
nontrivial `(1,2^r)` collision member of `(5)--(6)` can make all additional
critical roots disappear.  It does not, by itself, exclude a hypothetical
`A2` ring assembled from further observables.  The Chebyshev member is the
most economical case, coalescing every residual root into one critical value;
for this member `(19)` proves that the exact maximal polynomial-target-
observable completion is `Y_beta`, which is not `A2`.  Indeed a linear
rescaling identifies `(27)` with the surface of THM-3561, whose Picard group
is `Z` (whereas `Pic(A2)=0`).  No obstruction to an unrelated completion or
nonequivariant construction is asserted.

## 7. Exact symplectic form and the two-brackets-to-one collapse

Let

```text
Q(b)=b(b-beta),                 Q'=2b-beta.            (43)
```

The global symplectic form is represented on `c!=0` by

```text
omega=db wedge dc/c^2,            Phi_d^*omega=-dx wedge dq. (44)
```

It is regular and nonvanishing across `c=0`: differentiating
`c^2e=Q(b)` gives

```text
Q'(b)db=2ce dc+c^2de,
omega=de wedge dc/Q'(b)                                  (44a)
```

on either boundary chart, where `Q'(b)=+beta` or `-beta`.

Unlike the standard exponent-one Danielewski form, this exponent-two form
is exact.  Choose `U,T in C[b]` with

```text
UQ-TQ'=1;   for example U=-4/beta^2, T=-Q'/beta^2.    (45)
```

Then the polynomial one-form

```text
alpha=(U+T')ce db-T e dc-d(Tce)                       (46)
```

equals `db/c` on `c!=0`, and hence

```text
d alpha=omega                                         (47)
```

globally by regularity.

For comparison with the residual Darboux problem, use the opposite Poisson
sign convention (the negative of the bracket inverse to `omega`)

```text
{F,G}_Y=-Jac(F o Phi_d,G o Phi_d),                    (48)
```

so that

```text
{b,c}=c^2,       {c,e}=beta-2b,       {b,e}=-2ce.     (49)
```

Equations `(45)` and `(49)` give the exact two-brackets identity

```text
1={(U+T')ce,b}+{-Te,c}.                               (50)
```

Thus every odd tower member has the same bracket-width-two upper bound for
the unit.  A **single** pair `F,G in C[Y_beta]` with `{F,G}_Y=1` remains
open.  Such a pair would compose with `Phi_d` to a noninjective polynomial
Keller map of `A2`, and hence would be a genuine `JC(2)` counterexample;
neither `(50)` nor this theorem supplies one.

## 8. Fixed denominator powers: the sharp parity boundary

There is a separate exact classification inside the narrower ansatz

```text
D=1+x^2q,          a=q/D^n,          c=xG(D),          n>=0. (51)
```

For a prescribed nonzero constant `kappa`, the condition
`Jac(a,c)=-kappa` is equivalent to

```text
2D(D-1)G'+[n+(1-n)D]G=kappa D^(n+1).                 (52)
```

Every rational solution of `(52)` is polynomial.  Indeed, a pole of order
`s>=1` away from `D=0,1` leaves an uncancelled order-`s+1` derivative pole;
at `D=0` its leading coefficient is `(2s+n)`, and at `D=1` it is `1-2s`.
None can vanish.

For fixed nonzero `kappa`, a solution exists if and only if

```text
n=2m                                                    (53)
```

is even.  Equivalently, the nonzero solution pairs `(G,kappa)` form a
one-dimensional family; for fixed `kappa`, `G` is unique.  Writing the family
parameter as `lambda`, one has

```text
G=lambda D^m P_m(D),
P_m(D)=sum_(j=0)^m binom(2j,j) D^j/4^j,
kappa=lambda(2m+1)binom(2m,m)/4^m.                    (54)
```

To prove necessity, first note that a polynomial solution has degree exactly
`n`: if its degree `N` exceeded `n`, the leading coefficient on the left of
`(52)` would be `(2N+1-n)g_N`, which cannot vanish for `N>n`; degree below
`n` cannot supply the right side's `D^(n+1)` term.  Now write
`G=sum_(j=0)^n g_jD^j`.  For `0<=j<=n`, `(52)` gives

```text
(n-2j)g_j+(2j-n-1)g_(j-1)=0,          g_(-1)=0.       (55)
```

If `n` is odd there is no resonance, so `(55)` kills every coefficient,
contradicting the `D^(n+1)` term.  If `n=2m`, the coefficient of `g_m`
vanishes and leaves one free scalar.  Above the resonance, `(55)` has ratio
`(2j-1)/(2j)` after shifting by `m`, which is exactly the central-binomial
ratio in `(54)`; the top coefficient gives the displayed `kappa`.

For `m=1`, taking `lambda=2` recovers `G=D(D+2)` and `kappa=3`, the triple
member behind THM-3561.  For `m>1`, the polynomial-observable valuation
semigroup is richer; no three-generator completion analogous to `(27)` is
claimed here.  The parity classification is only a fixed-power ansatz
theorem.

## 9. Literature boundary and nonclaims

- Dubouloz--Poloni, [*On a class of Danielewski surfaces in affine
  3-space*](https://arxiv.org/abs/math/0602549), provides the generalized
  `x^n z-Q(x,y)` Danielewski-surface context.  The tower formulas and the
  maximal intersection `(19)` are repository derivations, with no priority
  claim.
- Kochetkov, [*Chebyshev polynomials, Zolotarev polynomials and plane
  trees*](https://arxiv.org/abs/1212.6317), is cited only for the historical
  generalized-Chebyshev/two-critical-value and plane-tree passport language.
- Abhyankar--Assi, [*Jacobian pairs*](https://arxiv.org/abs/math/0209159),
  studies meromorphic Jacobian pairs and Newton-polygon methods.  No theorem
  from that paper is used to turn `(16)` into a polynomial plane pair.
- Dubouloz--Kunyavskii--Regeta, [*Bracket width of simple Lie
  algebras*](https://arxiv.org/abs/2102.08674), Proposition 6(4), Proposition
  8, and its following equivalent-formulation lemma, is a comparison point
  for standard exponent-one Danielewski symplectic/bracket width.  Their
  standard form has nonzero de Rham class for a polynomial with at least two
  simple roots.  Our exact exponent-two form `(47)` belongs to a different,
  nonfinite affine modification, so their result does not apply directly.

The honest frontier is therefore the single-bracket equation on `Y_beta`,
or a genuinely different completion that avoids the ramification budget.
`JC(2)` remains open.

## 10. Exact companion and reproduction

The deterministic SymPy companion checks the master identity, all signs and
normalizations, the Chebyshev--Pell identities through odd degree `15`, the
smooth/etale charts, exact `d=3,5,7` fibre and nonproperness passports,
finite hostile controls for the valuation proof, the symplectic primitive,
the two-brackets identity, and the fixed-power parity system.  These finite
rows are controls for, not extrapolations to, the structural all-degree
proofs above.

```bash
python3 04-computation/jc2_chebyshev_pell_odd_keller_collision_tower_thm3566.py \
  > /tmp/thm3566.normal.out
python3 -O 04-computation/jc2_chebyshev_pell_odd_keller_collision_tower_thm3566.py \
  > /tmp/thm3566.optimized.out
cmp /tmp/thm3566.normal.out /tmp/thm3566.optimized.out
cmp /tmp/thm3566.normal.out \
  05-knowledge/results/jc2_chebyshev_pell_odd_keller_collision_tower_thm3566.out
sha256sum \
  04-computation/jc2_chebyshev_pell_odd_keller_collision_tower_thm3566.py \
  05-knowledge/results/jc2_chebyshev_pell_odd_keller_collision_tower_thm3566.out
```

Frozen universes:

```text
structural theorem: C[x,q], d=2r+1>=3, beta=d^-2;
exact symbolic controls: QQ[x,q,t,b,c,e]; d=3,5,...,15;
boundary controls: d in {3,5,7}, algebraic roots typed by P_d or S_d;
fixed-power controls: QQ[D], 0<=n<=12 with kappa!=0.
```
