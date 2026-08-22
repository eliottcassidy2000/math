---
id: THM-3581
title: "Critical-value multiarm Keller compiler and A13 carrier"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A squarefree critical polynomial with nonzero critical values, deduplicated
  by their distinct-value polynomial, compiles into an etale polynomial map
  from A2 to a smooth multiarm
  Danielewski modification.  The maximal polynomial and rational observable
  intersections, exact image, fibres, Jelonek curves, Picard rank, and
  sheet-debt vector are determined.  The cyclic n=3,r=2 instance has degree
  thirteen, geometric monodromy A13 and arithmetic monodromy S13.  Its two
  side origins are omitted and thirteen source points collide over every
  nonzero central-arm point.  This is a near-counterexample carrier: no
  polynomial Darboux pair and no counterexample to JC(2) is claimed.
source: root / delegated critical-value compiler lane, 2026-08-21
audit: >
  PASS.  An independent audit rederived the compiler, repeated-value
  divisibility, maximal polynomial and rational intersections, exact image,
  fibres, Jelonek curves, Picard rank, sheet debts, A13/S13 distinction,
  finite-sidecar obstruction, and reflection scope.  Normal, optimized, and
  stored-output replays are byte-identical across 440 exact gates.
depends_on:
  - THM-3578-zariski-main-boundary-rank-and-sheet-debt
related:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3567-separated-rational-keller-maximal-observable-nodal-completion-no-go
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
  - THM-3576-higher-exponent-belyi-keller-collision-tower
  - THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow
external:
  - "Cheng, McKay, and Wang, Younger mates and the Jacobian conjecture, Proc. Amer. Math. Soc. 123 (1995), 2939-2947, Theorem 1."
  - "Moskowicz, Involutions and the Jacobian conjecture, arXiv:1410.7705v1, Theorem 4.7 (exchange involution only)."
  - "Kambayashi, Automorphism group of a polynomial ring and algebraic group action on an affine space, J. Algebra 60 (1979), 439-451 (finite plane actions are linearizable)."
script: 04-computation/jc2_critical_value_multiarm_keller_a13_thm3581.py
output: 05-knowledge/results/jc2_critical_value_multiarm_keller_a13_thm3581.out
script_sha256: 109bdda2e95346b45ce91da2e060f3166aadd025f9f98971efc2de49fbaeac02
output_sha256: 5c3d3be949aac7e1ade3dde2347c04331e6ac35e8730bba69164fbfb1392df91
hash_basis: LF-normalized bytes
---

# THM-3581 -- critical-value multiarm Keller compiler and A13 carrier

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem does **not** refute `JC(2)`.  It gives a general exact compiler
for rational Keller pairs with polynomial multiarm completions, and a
degree-thirteen carrier whose inverse problem has geometric group `A13`.
The remaining task is still to find one polynomial Darboux pair on the target.

All varieties and geometric monodromy groups are over `C`.  Arithmetic
monodromy in Section 7 is over `Q(z)`.

## 1. The critical-value compiler

Fix `n>=2` and a nonconstant squarefree polynomial `S in C[t]` with
`S(0)!=0`.  There is a unique polynomial `P` satisfying

```text
P+n tP'=S^(n-1),
P(t)=integral_(u=0)^1 S(tu^n)^(n-1) du.                 (1)
```

Put

```text
B=tP^n,                   d=deg B=n(n-1)deg(S)+1.       (2)
```

Then

```text
B'=P^(n-1)S^(n-1).                                      (3)
```

For every root `sigma` of `S`, define

```text
beta_sigma=B(sigma).                                    (4)
```

Assume all these values are nonzero.  Since `S(0)!=0`, this says exactly
that `P` does not vanish at a root of `S`; hence `P,S` are coprime.  Every
root of `P` is then simple: if its multiplicity were `m`, the two sides of
`(3)` would have orders `nm-1` and `(n-1)m`, forcing `m=1`.  Let `V` be the
set of distinct critical values, and let

```text
k_beta=#{sigma:S(sigma)=0 and B(sigma)=beta},
E(z)=product_(beta in V)(z-beta).                       (5)
```

The values need not be pairwise distinct; using the distinct-value
polynomial `E`, rather than multiplying once per critical point, is
load-bearing.  Define

```text
W(t)=E(B(t))/S(t)^n.                                    (6)
```

At a simple root `sigma` of `S`, equation `(3)` gives

```text
ord_sigma(B-beta_sigma)=n,
W(sigma)=E'(beta_sigma)P(sigma)^(n-1)/(nS'(sigma))!=0.  (7)
```

Thus `W` is polynomial, squarefree, and coprime to `PS`.  For the
squarefreeness assertion, a repeated root away from `PS` would also be a
critical point of `B`, contrary to `(3)`; `(7)` handles the roots of `S`,
and a root of `P` maps to `0`, not to a root of `E`.  The fibre over
`beta in V` has exactly

```text
k_beta roots of multiplicity n and d-nk_beta simple roots. (8)
```

The latter number is positive because `d` is `1 mod n`.

## 2. Rational Keller pair and smooth polynomial completion

On `A2_(x,q)`, set

```text
t=x^nq,
a=q/S(t)^n,                 c=xP(t)S(t),
b=B(t),                     e=qW(t).                    (9)
```

For arbitrary rational functions `H,G` of `t=x^nq`, direct differentiation
gives the master identity

```text
Jac_(x,q)(qH,xG)=-(tHG^n)'/G^(n-1).                    (10)
```

Taking `H=S^(-n)` and `G=PS`, equations `(1)--(3)` give

```text
Jac_(x,q)(a,c)=-1,
b=ac^n,                         e=aE(b).                (11)
```

Consequently `(b,c,e)` is polynomial and defines

```text
Phi:A2_(x,q) -> Y,
Y=Spec C[b,c,e]/(c^n e-bE(b)).                         (12)
```

The polynomial `Sigma(b)=bE(b)` is squarefree, so `Y` is smooth.  Give `Y`
the everywhere symplectic Jacobian bracket

```text
{b,c}=c^n,
{c,e}=-(bE(b))',
{b,e}=-n c^(n-1)e.                                    (13)
```

Then `Phi` is anti-Poisson and etale.  It has generic degree `d`: for generic
`b`, the equation `B(t)=b` gives `d` values of `t`, after which
`x=c/(PS)` and `q=t/x^n` are forced.

## 3. The maximal observable intersections

Inside `C(x,q)` one has both

```text
C[a,c] intersect C[x,q]=C[b,c,e],                     (14)
C(a,c) intersect C[x,q]=C[b,c,e].                     (15)
```

For `(14)`, use the torus weights

```text
wt(x)=wt(c)=1,             wt(q)=wt(a)=-n,
wt(t)=wt(b)=0,             wt(e)=-n,       wt(c^n e)=0. (16)
```

A homogeneous element of `C[a,c]` of weight `-s<0` has the unique form
`c^(-s)f(b)`, with `b^m|f` and `m=ceil(s/n)`.  Along every divisor
`t=sigma`, regularity of its pullback and `(7)` force
`(b-beta_sigma)^m|f`.  Distinct critical points with the same value impose
the same factor, not a repeated one.  Hence

```text
[bE(b)]^m divides f,
c^(-s)f(b)=c^(nm-s)e^m h(b) in C[b,c,e].               (17)
```

Nonnegative weight pieces are `c^rC[b]`.  Exact weight decomposition proves
`(14)`, including the converse containment, without a cancellation gap.

For `(15)`, first note from `(11)` that

```text
Frac(Y)=C(b,c)=C(a,c).                                 (18)
```

Section 4 proves that `Phi` misses only finitely many points.  Every prime
divisor of the normal surface `Y` therefore has an upstairs prime divisor.
Etaleness identifies the corresponding divisorial valuations.  If
`f in C(a,c)` is polynomial in `x,q`, all its upstairs valuations are
nonnegative, hence every height-one valuation of `Y` is nonnegative.
Normality and the Krull intersection theorem give `f in O(Y)`.  This proves
the stronger rational intersection `(15)`.

## 4. Exact image, fibres, Jelonek set, and sheet debt

The set-theoretic image is

```text
Phi(A2(C))=Y(C) minus {(beta,0,0):beta in V}.           (19)
```

The complete affine-fibre invoice is:

| target stratum | affine fibre size |
|---|---:|
| `c!=0`, `b` away from `0` and `V` | `d` |
| `c!=0`, `b=0` | `1` |
| `c!=0`, `b=beta` | `d-nk_beta` |
| `b=c=0`, `e!=0` | `d` |
| `(b,c,e)=(0,0,0)` | `1` |
| `b=beta,c=0,e!=0` | `nk_beta` |
| `(b,c,e)=(beta,0,0)` | `0` |

For `c!=0`, this follows immediately by separating the roots of `B-b`
that lie on `P=0` or `S=0`.  On the central arm, `t=0` contributes one
point and every simple root of `P` contributes `n` points, for a total of
`d`; at the origin only `(x,q)=(0,0)` remains.  On a side arm, each of the
`k_beta` roots of `S` contributes `n` points when `e!=0`, while equation
`(7)` makes the side origin impossible.  This proves `(19)` as well.

For every `eta!=0`, all `d` points over `(0,0,eta)` have the same rational
Keller value

```text
(a,c)=(eta/E(0),0).                                    (20)
```

The nonproperness locus is exactly

```text
J(Phi)={b=0,e=0}
 union union_(beta in V)({b=beta,e=0} union {b=beta,c=0}). (21)
```

One direct proof works in the `(t,x)` chart, where
`q=t/x^n`, `c=xPS`, and `e=tW/x^n`.  Bounded `b=B(t)` keeps `t` bounded.
Escape with `x` tending to infinity forces `PS=0` and gives the central or
side `e=0` curves; escape with `x` tending to zero forces `W=0` and gives a
side `c=0` curve.  Conversely the displayed limiting paths realize generic
points of every line in `(21)`.

There are `1+2|V|` Jelonek components.  Since `bE(b)` has `1+|V|` simple
roots, the standard boundary-divisor calculation for
`c^n e=bE(b)` gives

```text
Pic(Y)=Cl(Y)=Z^|V|,                    O(Y)^*=C^*.      (22)
```

For completeness, if `alpha` ranges over the `1+|V|` roots of `bE(b)`, let
`D_alpha={c=0,b=alpha}` and `E_alpha={e=0,b=alpha}`.  Localization at `c`
is the UFD `C[b,c,c^(-1)]`, so the `D_alpha` generate the class group, while

```text
div(c)=sum_alpha D_alpha,
div(e)=sum_alpha E_alpha,
div(b-alpha)=nD_alpha+E_alpha.                         (22a)
```

These give exactly the single relation `sum D_alpha=0`, hence the free
group in `(22)`; the same valuation calculation leaves only constant units.

Applying THM-3578 to the generic fibre sizes along `(21)` gives the exact
missing-sheet debts

```text
central e=0 line:             d-1,
each side beta, e=0 line:     nk_beta,
each side beta, c=0 line:     d-nk_beta.               (23)
```

Their total is `(d-1)+|V|d`.  Picard rank and sheet debt are necessary
boundary resources, not obstructions when the eventual target is `A2`.

## 5. Exact symplectic primitive and two-bracket collapse

The compiler does not restore a de Rham obstruction.  More generally, let
`Sigma in C[b]` be squarefree and choose `U,T in C[b]` with

```text
U Sigma-T Sigma'=1.                                   (24)
```

On `c^n e=Sigma(b)`, put

```text
alpha_0=Uce db-nTe dc-Tc de.                          (25)
```

The surface relation gives `c^(n-1)alpha_0=db`, so on `c!=0`

```text
d(alpha_0/(n-1))=db wedge dc/c^n.                     (26)
```

Both sides are global regular two-forms, so equality on this dense chart is
global.

Rewriting `(25)` yields the polynomial identity

```text
1={((U+T')/(n-1))ce,b}+{-Te,c}.                       (27)
```

Thus the constant has Poisson-bracket length at most two.  A counterexample
still requires the much sharper one-bracket equation `{F,G}=1`.

### 5.1 A Kummer-sidecar one-coordinate no-go

There is nevertheless an all-exponent obstruction that is invisible to a
bounded weight-cell search.  Put `z=ce`; then

```text
{z,b}=(n-1)Sigma(b).                                  (27a)
```

Let `phi,psi in C[b]`, with `psi!=0`, and let `F` be a nonconstant
one-variable polynomial.  No rational function `P in Frac(Y)` can satisfy

```text
{P,F(phi(b)+z psi(b))}=1.                             (27b)
```

Indeed, write `w=phi(b)+zpsi(b)` and base-change the generic `F(w)`-fibre
to `C(w)`.  On it,

```text
c^(n-1)=Sigma(b)psi(b)/(w-phi(b)),
{b,F(w)}=-(n-1)F'(w)Sigma(b)psi(b).                   (27c)
```

This is a finite separable Kummer extension of `C(w)(b)` of degree `n-1`
(degree one when `n=2`).  Equation `(27b)` would make

```text
-db/[(n-1)F'(w)Sigma(b)psi(b)]                        (27d)
```

exact upstairs.  Trace of differentials makes it exact downstairs.  After
removing the nonzero constants, this would give a rational primitive of
`db/H(b)` for `H=Sigma psi`.

But a reciprocal polynomial `1/H` has a rational primitive only when `H`
has one distinct root.  To see this, a simple root already gives a nonzero
logarithmic residue.  Otherwise, if the root multiplicities are `m_i>=2`, a
rational primitive has poles of orders `m_i-1`, hence degree
`deg(H)-#{roots of H}` as a rational map, while its expansion at infinity
has local degree `deg(H)-1`.  With at least two roots the local degree would
exceed the global degree.  Since squarefree `Sigma=bE(b)` has at least the
central and one side root, `(27b)` is impossible.

This kills every coordinate obtained by a one-variable outer composition
after a single `ce`-linear mixing, even after adjoining the full Kummer
sidecar.  It does not cover nonlinear dependence on `ce` or two genuinely
mixed Darboux coordinates.

## 6. A cyclic family with many distinct arms

Take

```text
S=1+t^r,
P_(n,r)=sum_(j=0)^(n-1) binom(n-1,j)t^(rj)/(nrj+1).   (28)
```

At every root `sigma^r=-1`,

```text
p=P_(n,r)(sigma)=product_(j=1)^(n-1) (jnr)/(jnr+1),
beta_sigma=sigma p^n.                                 (29)
```

The product formula follows by evaluating `(1)` as a beta integral,
`integral_0^1(1-u^(nr))^(n-1)du`.  These values are nonzero and pairwise
distinct.  Therefore

```text
E(b)=b^r+p^(nr),
B=tP_(n,r)^n,
W=(B^r+p^(nr))/(1+t^r)^n,
d=n(n-1)r+1.                                          (30)
```

The case `r=1` is THM-3576's Shabat tower.  For `r>=2`, the finite critical
values are `0` and the `r` distinct side values, so `B` is not Belyi.

## 7. The degree-thirteen A13 carrier

The cyclic multiarm row `n=3,r=2` supplies an odd simple monodromy carrier:

```text
S=1+t^2,
P=(7t^4+26t^2+91)/91,
kappa=(72/91)^3=373248/753571,
B=tP^3,
W=(B^2+kappa^2)/(1+t^2)^3,                            (31)
Y_13: c^3e=b(b^2+kappa^2).                            (32)
```

The degree-thirteen cover `B:A1_t->A1_b` has passports

```text
over 0:          (3^4,1),
over +i kappa:   (3,1^10),
over -i kappa:   (3,1^10),
over infinity:   (13).                                (33)
```

Its geometric monodromy is `A13`.  Indeed, prime-degree transitivity makes
the group primitive, a side inertia generator is a literal `3`-cycle, and
Jordan's theorem gives `A13` inside the group.  Every inertia permutation
in `(33)` is even, giving the reverse containment.

Exact elimination gives

```text
Disc_t(B(t)-z)=13^(-23) z^8(z^2+kappa^2)^2.           (34)
```

The rational constant `13^(-23)` is nonsquare.  Hence arithmetic monodromy
over `Q(z)` is not contained in `A13`; since it contains the geometric
`A13`, it is exactly `S13`.  Geometric `A13` and arithmetic `S13` are not
interchangeable labels.

The map misses exactly `(i kappa,0,0)` and `(-i kappa,0,0)`.  Generic and
nonzero central-arm fibres have size `13`; a side value with `c!=0` has size
`10`, a nonzero side-arm fibre has size `3`, the central special fibres have
size `1`, and each omitted side origin has size `0`.  In particular, for
every `eta!=0`, the thirteen points above `(0,0,eta)` are

```text
(x,q)=(0,eta/kappa^2),

and, for each root rho of P, the three choices satisfying
x^3=rho kappa^2/[eta(1+rho^2)^3],
q=eta(1+rho^2)^3/kappa^2.                             (35)
```

Here `Pic(Y_13)=Z^2`.  A convenient exact instance of `(27)` is

```text
1={-15bce/(4kappa^4),b}
  +{(1/kappa^2+3b^2/(2kappa^4))e,c}.                  (36)
```

### 7.1 No partial algebraic sheet-label sidecar

There are two distinct sidecar obstructions.  First, universally, let `L`
be any finite extension of `Frac(Y)` inside a source function field and let
`Y_L` be the normalization of `Y` in `L`.  Finite conorm and norm give

```text
rank Cl(Y_L) >= rank Cl(Y)=|V|.                        (37)
```

Indeed the kernel of conorm is torsion, while `Cl(Y)` is free.  Thus no
finite algebraic sidecar normalization is factorial or `A2`, even when the
inverse monodromy is imprimitive.  To reach the polynomial source plane one
must still delete boundary divisors: the escape is a nonfinite affine
modification/localization, not a finite sheet-label extension.

Second, the A13 row has the sharper field-theoretic dichotomy

```text
Frac(Y_13)=C(a,c),
C(x,q)=Frac(Y_13)(t),                 B(t)=b.           (37a)
```

After the purely transcendental base change by `c`, the geometric Galois
closure is still `A13`.  The stabilizer of one sheet is `A12`, which is
maximal in `A13`.  The degree-thirteen extension in `(37a)` therefore has no
proper intermediate fields.  An algebraic observable field between target
and source cannot retain only some sheet labels: it is either `Frac(Y_13)`
or all of `C(x,q)`.  This does not constrain subfields such as `C(F,G)`
inside `Frac(Y_13)` and does not solve the Darboux equation.

## 8. A reflection-specific symmetry no-go

The A13 map intertwines the explicit involutions

```text
tau(x,q)=(x,-q),
iota(b,c,e)=(-b,c,-e),
Phi tau=iota Phi.                                      (38)
```

The source map `tau` is a reflection.  It is linearly conjugate to the
exchange involution `alpha(u,v)=(v,u)` by

```text
ell(x,q)=(x+q,x-q),                    ell tau=alpha ell. (39)
```

Suppose `{F,G}=1` on `Y_13` and either `F` or `G` is an `iota`-eigenvector.
Then `(F Phi,G Phi)` is a planar Keller pair and one coordinate is a
`tau`-eigenvector.  Conjugating by `(39)` reduces exactly to Moskowicz,
arXiv:1410.7705v1, Theorem 4.7: a planar Keller pair with one coordinate
symmetric or skew under exchange is an automorphism.  Equivalently, the
exchange case follows by the Cheng--McKay--Wang centralizer theorem, a target
shear making the mate opposite parity, and the exchange-equivariant Keller
theorem.  But `(35)` gives thirteen points with the same pair value, a
contradiction.  Therefore

```text
each coordinate of every Darboux pair must mix iota parities.            (40)
```

The reflection hypothesis is load-bearing.  This proof does **not** invoke
Moskowicz's arbitrary-involution Lemma 2.1 or Theorem 4.9, whose conjugacy
scope would be too broad here.

There is a coordinate-free strengthening.  If a hypothetical Darboux
subalgebra `D=C[F,G]` were `iota`-stable, its induced involution on
`Spec D=A2` would have Jacobian determinant `-1`, because `iota` is
anti-Poisson.  The classical linearization theorem for finite polynomial
actions on `A2` conjugates this involution to a linear reflection.  Reflection
eigen-coordinates would then give an `iota`-eigen Darboux coordinate,
contradicting the preceding paragraph.  Hence

```text
iota(D)!=D                                                   (41)
```

for every hypothetical Darboux subalgebra.  This is stronger than saying
that one chosen presentation mixes parity, but it is still not a proof that
no Darboux pair exists.

## 9. The tempting Belyi quotient loses the arm label

For the A13 row,

```text
C(t)=-B(t)^2/kappa^2                                  (42)
```

is a degree-twenty-six Shabat polynomial with passports

```text
over 0:        (6^4,2),
over 1:        (3,3,1^20),
over infinity: (26).                                  (43)
```

But `C=h(B)` factors through a degree-two map.  Its monodromy is imprimitive
with two blocks of thirteen and it merges the `+i kappa` and `-i kappa`
arms.  On the surface, the corresponding involution

```text
(b,c,e) -> (-b,-c,e)                                  (44)
```

is anti-Poisson.  The quotient therefore loses both the side-arm label and
the symplectic orientation.  The Shabat reduction is a useful inverse
carrier, not a faithful replacement for the three-arm target.

## 10. Hypothesis hostiles

Each compiler hypothesis has a cheap exact failure mode.

1. **Repeated critical values are harmless when deduplicated.**  For `n=2`,
   take

   ```text
   P=t^2+5t+5,                 S=5t^2+15t+5.
   ```

   Both roots of `S` have critical value `-4`, and

   ```text
   B+4=(t+4)S^2/25.                                    (45)
   ```

   The correct polynomial is `E=b+4`, not `(b+4)^2`.  The resulting target
   `c^2e=b(b+4)` remains smooth.

2. **Squarefreeness is essential.**  For `n=2,S=(1+t)^2`, one gets

   ```text
   B-beta=(t+1)^3(9t^2+33t+64)/225,                   (46)
   ```

   which is not divisible by `S^2=(1+t)^4`.  Adding enough multiplicity to
   clear the denominator repeats a target factor and destroys smoothness.

3. **A zero side value collides with the central arm.**  For

   ```text
   n=2,       P=(t-1)^2,       S=(t-1)(5t-1),          (47)
   ```

   the root `t=1` has `beta=0`.  Clearing the completion repeats `b`, giving
   a singular model such as `c^n e=b^2(...)`.  Even normalization does not
   recover an etale symplectic plane: for the local model
   `c^(2k)e=b^2`, normalization `b=c^k r,e=r^2` changes the physical form to
   `dr wedge dc/c^k`, equivalently `{r,c}=c^k`, which degenerates on `c=0`.
   Nor can the tempting sidecar `r=b/c` globally cancel the bad roots.
   At a squarefree nonzero root with `beta=0`, `(1)` forces
   `ord_sigma(P)=n`; if every root of `S` behaved this way then `S^n|P`,
   impossible because `deg P=(n-1)deg S`.

4. **The origin condition is not cosmetic.**  If `S(0)=0`, then the central
   point itself is one of the critical arms and its value is zero, reducing
   to the preceding singular collision.

5. **A constant `S` is only the plane boundary.**  Then `B` has degree one,
   there are no collisions or side arms, and the construction is a rational
   reparametrization of `A2`.

6. **Characteristic zero is used.**  The integral/denominators in `(1)`,
   separability, Jordan's theorem, and the reflection theorem are all
   load-bearing.

## 11. Typed connection and exact scope

The construction can be summarized as follows.

```text
source:       A2_(x,q)
target:       Y={c^n e=bE(b)}
map:          (x,q) |-> (B(x^nq), xPS, qW)
preserved:    etaleness, the symplectic bracket up to sign, and Frac(Y)=C(a,c)
lost:         d-valued inverse sheet labels and finitely many side origins
sidecars:     t with B(t)=b; arm labels beta; finite-envelope boundary debts
cheap tests:  ODE/B' identity, exact divisibility and gcds, fibre degrees,
              discriminant/passports, bracket identity, reflection parity
```

The source is `A2`, but the target is the nonfactorial surface `Y`, not
`A2`.  A polynomial Darboux pair `{F,G}=1` on `Y` would compose with `Phi`
to give a genuine planar Keller counterexample, but no such pair is
constructed.  Exactness, bracket width at most two, large monodromy, missing
points, and many collisions are all near-counterexample data; none is a
counterexample to `JC(2)`.

## 12. Exact verification contract

The companion performs exact symbolic checks, with no Python `assert`, for:

- cyclic rows `2<=n<=5`, `1<=r<=3`, plus noncyclic `S=1+t+t^2`;
- `(1)--(13)`, divisibility, gcd, degree, squarefreeness, target, Jacobian,
  and special-fibre degree invoices;
- the complete A13 row, discriminant `(34)`, passports, collision formulas,
  two-bracket identity `(36)`, and symmetry gates `(38)--(44)`;
- the repeated-value positive control `(45)` and the nonsquarefree and
  zero-value hostiles `(46)--(47)`.

The universal rational-intersection, image, Jelonek, Picard, and Zariski-main
claims are proof-driven; the finite exact atlas is a hostile control, not an
extrapolation.

Reproduce with

```bash
python3 04-computation/jc2_critical_value_multiarm_keller_a13_thm3581.py
python3 -O 04-computation/jc2_critical_value_multiarm_keller_a13_thm3581.py
```
