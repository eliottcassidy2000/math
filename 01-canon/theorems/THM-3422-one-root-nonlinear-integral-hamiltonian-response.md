---
id: THM-3422
title: "One-root nonlinear integral Hamiltonian response"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let K have characteristic
  zero, d>=2, and
  P=ax+b+c(x-alpha)^e z^d with a,c nonzero.  If e=1, every Hamiltonian
  target sector is a free rank-one K[P]-module.  If e>1, target sector s-1
  is K[lambda] plus one lambda-Pruefer arm when d divides s(e-1), and is the
  Laurent line K[lambda,lambda^-1] otherwise, where
  lambda=(P-(a alpha+b))/a.  Thus exactly gcd(d,e-1) characters carry
  Pruefer torsion, all with filtration slope e-1.  The unit class [1] has
  annihilator (lambda^((e-1)/d)) when d divides e-1 and annihilator zero
  otherwise.  A separately typed split-root first-window formula is included;
  no global multiple-root decomposition, new Keller case, or JC(2) conclusion
  is claimed.
source: root-2608-jc-integral-response-2026-08-15
audit: independent line rederivation of the orbit, wrap, Pruefer, thickness, observer, and split-window arguments clean; exact companion normal/optimized/stored outputs byte-identical; hashes, routing, and documentation audit clean
depends_on:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3419-generic-kummer-response-regular-sector-rank
related:
  - THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms
  - THM-3348-linear-z-generic-puncture-response-and-one-root-valuation
script: 04-computation/jc_one_root_nonlinear_integral_response_thm3422.py
output: 05-knowledge/results/jc_one_root_nonlinear_integral_response_thm3422.out
script_sha256: 1c862dd5dfaee00a3ee5827a9d004c160a574deb103cd3b878a6615aeac7a766
output_sha256: 1e0abbff4bbd2ab601cbf2aee4c840789f4ebf36a1444bf2aa8c8506579e485d
hash_basis: LF-normalized bytes
---

# THM-3422 -- one-root nonlinear integral Hamiltonian response

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Inheritance and connection contract

[THM-3418](THM-3418-one-monomial-nonlinear-fiber-keller-classification.md)
gives an exact colimit presentation for every fiber-exponent sector of

```text
C_P=K[x,z]/D_P(K[x,z]),             D_P=Jac(P,-),          (1)
```

when `P=ax+b+g(x)z^d`.  [THM-3419](THM-3419-generic-kummer-response-regular-sector-rank.md)
proves that, for `N=deg(rad(g))`, all `d` sectors have generic `K[P]`-rank
`N`.  Generic rank does not determine the integral lattice.  The closest
corrected hostile is MISTAKE-374: localization may kill the distinguished
unit class while leaving it nonzero integrally.

For one geometric root the missing integral coordinate is elementary:

```text
ell = (base monomial degree) - (fiber depth)(e-1).         (2)
```

It turns the sector colimit into a weighted bilateral chain.  Multiplication
by `P-beta` either has no break or has exactly one break; those alternatives
are respectively a Laurent line and a free line plus a Prüfer arm.

| item | exact content |
|---|---|
| source | a THM-3418 integral sector for one-root `g` |
| target | its full `K[P]`-module, generic line, and torsion filtration |
| map | the orbit decomposition `(2)` and multiplication by `P-beta` |
| preserved by generic localization | sector label and rank one |
| destroyed by generic localization | primary torsion, integral divisibility, and depth |
| required sidecar | `ell` and the wrap-sector quotient bound |
| cheapest hostiles | `(d,e)=(2,2)` and `(2,3)` have equal generic ranks but different integral sectors |

The applicable `META-PATTERNS` card is **Separate descent, ambient scale, and
regularity debt**: the regular generic packet must be followed by a separate
integral intersection calculation.

## 2. Main one-root statement

Let `K` be a field of characteristic zero.  Fix `d>=2` and

```text
P=ax+b+c(x-alpha)^e z^d,
a,c in K*,                    alpha,b in K,               e>=1.     (3)
```

Put

```text
y=x-alpha,                   beta=a alpha+b,
lambda=(P-beta)/a=y+gamma y^e z^d,
gamma=c/a,                   h=e-1.                       (4)
```

Thus `K[P]=K[lambda]`.  Let `C_(s-1)` be the Hamiltonian target sector of
fiber weight `s-1 mod d`, with `1<=s<=d`; `s=d` names the wrap sector.

### Simple root

If `e=1`, every sector is free:

```text
C_(s-1) ~= K[lambda].                                      (5)
```

### Repeated root

If `e>1`, then

```text
C_(s-1) ~=
  K[lambda] direct-sum K[lambda,lambda^(-1)]/K[lambda]
      if d divides s(e-1),

  K[lambda,lambda^(-1)]
      if d does not divide s(e-1).                         (6)
```

The direct-sum coordinates in `(5)--(6)` require nonzero scalar choices along
the monomial chains.  The sector decomposition, its torsion submodule, and
the filtration below are canonical.  Writing

```text
g_0=gcd(d,e-1),
Pr_lambda=K[lambda,lambda^(-1)]/K[lambda],
```

one obtains the abstract module decomposition

```text
C_P ~= (K[lambda] direct-sum Pr_lambda)^g_0
       direct-sum K[lambda,lambda^(-1)]^(d-g_0),            (7)

tors_(K[P])(C_P) ~= Pr_lambda^g_0.                          (8)
```

Every displayed sector has `K(lambda)`-rank one.  Thus `(7)` recovers
THM-3419's one-copy regular generic packet while showing that the integral
characters are not uniform.

### Exact torsion thickness

In a selected sector set

```text
q_s=s(e-1)/d.                                               (9)
```

Let `F_(s,k)` be the image of the depth-`k` THM-3418 stage, whose target
fiber exponent is `s-1+kd`.  Then

```text
dim_K(F_(s,k) intersect tors(C_(s-1)))
 =k(e-1)+q_s.                                              (10)
```

Hence every selected character has slope `e-1`, while the intercept depends
on `s`.  The selected target characters are the shifted subgroup

```text
{s-1 mod d : s(e-1)=0 mod d},                              (11)
```

of size `gcd(d,e-1)`.

### The unit observer

For `theta=[1] in C_0`,

```text
Ann_(K[lambda])(theta)=
  (lambda^((e-1)/d))       if e>1 and d divides e-1,
  (0)                      otherwise.                     (12)
```

In particular `theta` is never zero for nonconstant one-root `g`, but

```text
theta tensor K(lambda)=0
 iff e>1 and d divides e-1.                                (13)
```

Thus generic vanishing in `(13)` is not a polynomial mate.  It is a precise
nonlinear instance of the localization loss corrected by MISTAKE-374.

If `e-1=qd`, the annihilating power has the explicit primitive

```text
Q_q=sum_(j=0)^(q-1)
      [binom(q-1,j)/(1+jd)] gamma^j
      y^(q+j(e-1)) z^(1+jd),                               (14)

(D_P/a)(Q_q)=lambda^q.                                    (15)
```

The proof below shows that `(12)` is minimal, not just an upper bound.

## 3. Monomial-orbit proof

Divide `D_P` by the unit `a`.  Direct differentiation yields

```text
delta=D_P/a
 =(1+gamma e y^(e-1)z^d) partial_z
  -gamma d y^e z^(d-1) partial_y.                          (16)
```

For `n=s+kd`, write

```text
E_(k,m)=[y^m z^(s-1+kd)].                                  (17)
```

Applying `(16)` to `y^m z^n` gives exactly

```text
n E_(k,m)+gamma(ne-dm)E_(k+1,m+e-1)=0.                    (18)
```

For `s<d`, the stage has `m>=0`.  In the wrap sector, the exponent-zero
boundary relation gives the exact quotient

```text
0<=m<e(k+1).                                               (19)
```

This is the load-bearing wrap sidecar; omitting `(19)` creates false
resonances.

Assume first that `h=e-1>0`.  Relation `(18)` preserves

```text
ell=m-kh.                                                  (20)
```

Every `ell in Z` occurs at all sufficiently large stages, including under
`(19)`.  On its one-dimensional chain, the transition numerator is

```text
dm-en=d ell-es-dk.                                         (21)
```

It vanishes at most once.  If it vanishes, the initial segment dies and the
tail remains; if it does not, the whole chain survives.  Either way, the
colimit contributes exactly one line for each `ell in Z`.

Multiplication by `lambda` sends the `ell` line to the `ell+1` line.  Indeed,
put

```text
A=E_(k,m+1),                   B=E_(k+1,m+e).
```

Then `lambda E_(k,m)=A+gamma B`, while `(18)` applied to the input
`y^(m+1)z^n` gives, at a sufficiently late stage,

```text
lambda E_(k,m)
 =[d(m+1)-n(e-1)]/[d(m+1)-ne] A.                          (22)
```

The numerator is independent of stage:

```text
d(m+1)-n(e-1)=d(ell+1)-s(e-1).                            (23)
```

Therefore a `lambda` arrow vanishes exactly at

```text
ell_0=s(e-1)/d-1                                           (24)
```

when that number is integral.  With no zero arrow, the weighted bilateral
chain rescales to `K[lambda,lambda^(-1)]`.  With one zero arrow, the half
`ell<=ell_0` is `Pr_lambda` and the half `ell>=ell_0+1` is `K[lambda]`.
This proves `(6)`.

If `e=1`, the label is instead `ell=m>=0`.  No ordinary-sector transition
can vanish because `1<=s<d`; in the wrap sector the `m` line enters at depth
`m` and persists.  Formula `(23)` becomes `d(m+1)`, so every arrow is
nonzero.  The resulting unilateral chain is `K[lambda]`, proving `(5)`.

The congruence in `(24)` has `gcd(d,e-1)` solutions modulo `d`, proving
`(7)--(8)` and `(11)`.  In a selected sector, its torsion endpoint has the
direct identity

```text
tau_s=[y^(q_s-1)z^(s-1)],
lambda tau_s=delta(y^q_s z^s/s)=0.                         (25)
```

For a selected nonwrap sector, `(21)` cannot also vanish: divisibility of
both `s(e-1)` and `se` by `d` would imply `d|s`.  In the wrap sector its
formal zero is the excluded monomial `m=e(k+1)`.  Thus the depth-`k` image
meets the torsion half-chain in

```text
-k(e-1)<=ell<=ell_0,
```

which has the length `(10)`.

Finally `[1]` lies at `s=1,ell=0`.  If `d` does not divide `e-1`, it belongs
to a Laurent line and has zero annihilator.  If `e-1=qd`, the sole break is
at `ell=q-1`, so exactly `q` applications of `lambda` kill it.  This proves
the minimality in `(12)`.  In `(14)`, both the low coefficient of the `j`th
term and the high coefficient of the preceding term are the adjacent
`(q-1)`st binomial coefficients.  Pascal's identity gives `(15)` directly.
This completes the one-root proof candidate.

## 4. Separately typed split-root first-window sidecar

This section is an exact quotient calculation, not a global multiple-root
module decomposition.  Let `K'/K` be a faithful splitting extension which
also contains the `d`th roots of unity, put `C'=C_P tensor_K K'`, and write

```text
g(x)=c product_(j=1)^N (x-alpha_j)^e_j                  (26)
```

with distinct roots.  Put `beta_i=a alpha_i+b` and fix a repeated root
`e_i>1`.  For `1<=s<=d`, set

```text
v_i=e_i-1,                 v_j=e_j for j!=i,
c_i=gcd(d,v_1,...,v_N).                                  (27)
```

Then

```text
dim_(K') C'_(s-1)/(P-beta_i)C'_(s-1)
 =N       if d divides s v_j for every j,
 =N-1     otherwise.                                     (28)
```

Exactly `c_i` characters take the value `N`.

To prove this, write `y=x-alpha_i`.  The special fiber factors as

```text
P-beta_i
 =y [a+c y^(e_i-1) product_(j!=i)(x-alpha_j)^e_j z^d].    (29)
```

The factors are comaximal because `e_i>1`.  On the vertical component,
`D_P=a partial_z`, so its cokernel is zero.  On the other component, `y`,
`z`, and all `x-alpha_j` are units, and one has the finite etale Kummer cover
of `A1\{alpha_1,...,alpha_N}`

```text
z^d=-a/[c y^(e_i-1) product_(j!=i)(x-alpha_j)^e_j].        (30)
```

Here `P_z` is a unit, so the THM-3419 one-form/Hamiltonian map is an
isomorphism.  Target weight `s-1` corresponds to de Rham character `s`.
That rank-one local system has monodromy exponents `s v_j mod d`.  It is
trivial exactly in the first case of `(28)`.  The base is a wedge of `N`
circles: the trivial local system has `H^1` dimension `N`, while every
nontrivial one has dimension `N-1`.  This proves `(28)`.  The equivalent
covering graph has `c_i` components and Betti number

```text
dN-d+c_i=d(N-1)+c_i,                                      (31)
```

which is the sum of the `d` character dimensions.

The smallest gluing hostile is

```text
d=4,                    (e_i,e_j)=(3,1).                  (32)
```

Character `s=2` is resonant for the chosen root because
`4|s(e_i-1)`, but the other exponent `s e_j=2` is nontrivial.  Formula `(28)`
therefore gives dimension `N-1=1`, not the dimension two predicted by a
naive direct sum of independent root channels.  Only `s=4` has dimension
two.  Distinct boundary values do not collide because
`beta_i=a alpha_i+b` is injective; the obstruction is global Kummer-period
coupling.

No higher `(P-beta_i)`-power, torsion-kernel, Prüfer-persistence, or global
multiple-root decomposition is inferred from `(28)`.

## 5. Boundaries and JC scope

- If `g` is a nonzero constant, or if `g=0`, then `K[x,z]=K[P,z]` and
  `D_P=a partial_z`; the response is zero.  These cases are not `e=0`
  instances of the theorem.
- The same proof at `d=1` gives one free line plus one Prüfer arm for `e>1`,
  recovering the one-root boundary of THM-3412.  The theorem is stated for
  `d>=2` to refine THM-3418/3419.
- The split first-window proof requires `e_i>1`.  At a simple chosen root the
  components in `(29)` meet at critical points and the proof does not apply.
- Generic vanishing of `[1]` in `(13)` is not integral vanishing and gives no
  mate.  THM-3418 already proves the complete Keller classification for this
  sparse family.
- Nothing here treats intermediate fiber coefficients, constructs a general
  affine-modification chart, or proves `JC(2)` or `DC(2)`.

## 6. Exact companion

The standard-library companion uses only integer and `Fraction` arithmetic.
It checks:

```text
34,504  direct bivariate monomial relations,
   280  direct torsion-endpoint identities,
    50  explicit unit-observer primitives,
30,576  orbit-transition identities,
30,576  lambda-arrow and break identities,
 2,268  finite torsion-thickness windows,
 1,260  simple-root controls,
21,060  repeated split-root packets, and
126,360 split-root special-fiber character dimensions.    (33)
```

The last family is independently matched to the Kummer covering-graph
component and Betti counts.  Ordinary and optimized executions are
byte-identical.

Reproduce with

```text
python3 04-computation/jc_one_root_nonlinear_integral_response_thm3422.py
python3 -O 04-computation/jc_one_root_nonlinear_integral_response_thm3422.py
```

Artifact hashes are pinned in the frontmatter.
