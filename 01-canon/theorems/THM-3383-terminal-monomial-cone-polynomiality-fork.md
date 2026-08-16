---
id: THM-3383
title: "Terminal monomial-cone polynomiality fork and initial-ring module"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let e,g be coprime
  positive integers, a>=0, and n=g-ae nonzero.  The normalized terminal
  monomial chart R=u^a r^g, M=1+[e/n]u r^e, x=R^(-1), y=MR satisfies
  dx wedge dy=du wedge d(r^e).  Its residue coordinate generates C(u), and
  both targets rationally decode in C(x,y), exactly when |n|=1.  On that
  locus the sign of n gives a sharp polynomiality fork and an exact
  intersection module: for n=1 the source-polynomial target subring is
  C[u,ut,y^e], while for n=-1 it is C[r^e,ut,y^e], with one explicit
  hypersurface relation in either case.  Thus no polynomial target
  automorphism makes both target coordinates polynomial in this family.
  On the additive quotient of the target polynomial ring by this
  intersection, multiplication by the surviving polynomial coordinate is
  locally nilpotent but has unbounded exact response lengths: t^q has length
  q in the positive orientation, and u^q has length q in the negative one.
  The normalized C3 locus has an explicit two-congruence-class atlas.  The
  4+3 and 4+4+3 controls are opposite sides of the fork after polynomial
  source shears, and depth composition itself is not source-automorphism
  invariant.  No arbitrary terminal module, full C3, or JC(2) exclusion is
  asserted.
source: >
  root/jc2-other-angles-2026-08-14; operation-response strengthening
  root/factorial-jacobian-alternation/2026-08-15
depends_on:
  - THM-3081-terminal-toric-residue-parameter-mobius-rigidity-and-autonomous-decoder
related:
  - THM-3070-polynomial-c3-one-face-escape-leading-cancellation-gate
  - THM-3074-c3-two-pole-binomial-cancellation-and-first-key-form-depth-lattice
  - THM-3080-c3-finite-toric-key-tower-depth-partition-and-gcd-descent
  - THM-2690-normal-crossing-cyclic-cubic-resolvent-exclusion-and-reflection-completion-boundary
  - THM-3397-torsor-killing-versus-effective-boundary-valuations
  - THM-3466-factorial-face-stokes-and-keller-boundary-current
script: 04-computation/jc_terminal_monomial_cone_polynomiality_fork_thm3383.py
output: 05-knowledge/results/jc_terminal_monomial_cone_polynomiality_fork_thm3383.out
script_sha256: b72a196e7fda4a2de33df26e1f26f06f654908a9eb42f453e210728559ad4b10
output_sha256: 514af913271c9fa283bf7c969f9166730a3a187e93880b65d9a480c8966f2abc
semantic_sha256: 883e0b685f2cac9c329b3181e5af86009cd5425d3ebe5d3252f065ac92254dd0
hash_basis: LF-normalized bytes
---

# THM-3383 -- terminal monomial cones retain an exact polynomial ring sidecar

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Inheritance and result

THM-3070 proves that simultaneous polynomiality eliminates the first
nonresonant one-face `C3` escape.  THM-3074/3080/3081 reaches a primitive
terminal Laurent chart and an
autonomous Mobius residue decoder, but deliberately forgets which rational
initial coefficients actually lie in `C[x,y]`.

This theorem restores that ring sidecar on the simplest complete terminal
family.  The depth, gcd, residue-field and Keller-form tests all pass, yet the
two decoded targets lie on opposite sides of the polynomial exponent cone.
The complete intersection with `C[x,y]` is an explicit hypersurface ring.
The missing datum is therefore not another scalar depth: it is the oriented
initial subalgebra and its two boundary divisibilities.

## 2. The normalized terminal monomial chart

Let `e,g>=1` be coprime, let `a>=0`, and put

```text
n=g-ae !=0,                  c=e/n.                    (1)
```

In the rational function field `C(u,r)`, define

```text
R=u^a r^g,
M=1+c u r^e,
x=R^(-1),                    y=MR,
t=r^e,                       T=xy-1.                  (2)
```

The parameters are normalized: the `r`-adic value group of `C(x,y)` is
`Z` when the rational-decoding gate below holds, and the target parameter
then has tame index `e`.  An unnormalized presentation with ambient
exponents `E=de` and `G=dg` is the same packet after replacing `r` by
`s^d`; it must not be counted as a new ramification-index-`E` branch.

The functions `x,y` are algebraically independent because their Jacobian is
nonzero.  Directly,

```text
T=M-1=c u r^e,                                      (3)

dx wedge dy
 =dM wedge dlog R
 =c(g-ae)r^(e-1) du wedge dr
 =e r^(e-1) du wedge dr
 =du wedge d(r^e).                                  (4)
```

Thus `(2)` is an exact ambient symplectic parametrization.  Before the
gate in Section 3, `u,t` need not belong to `C(x,y)` at all.  In the toric
presentation

```text
dx wedge dy=M dlog M wedge dlog R,                    (5)
```

the prefactor has value zero, `w(R)=g`, and `w(M-1)=e`.  Hence the
deviation consumes the whole index-`e` differential budget at one terminal
stage, and `gcd(g,e)=1` is exactly THM-3080's terminal primitivity.

The polynomial shear `(x,y)->(x,x+y)` turns `(2)` into a two-pole packet
without changing the source ring or the two-form.  Accordingly the result is
not tied to the visual distinction between one-pole and equal-leading-pole
coordinates.

## 3. Residue and rational-decoding gates coincide

Put

```text
v=T/c=ut,                    L=1+cv=xy.                (6)
```

The primitive terminal value-zero ratio is

```text
Theta=T^g/R^e,
theta=res(Theta)=c^g u^n.                              (7)
```

Consequently

```text
C(theta)=C(u)  iff  |n|=1.                             (8)
```

Indeed `C(u^n)=C(u)` for a nonzero integer `n` exactly when `|n|=1`.
The same equation is the rational-decoding gate.  The exponent vectors of
`x,v` in the torus `(u,r)` are `(-a,-g)` and `(1,e)`, with determinant
`n`.  Their inverse gives

```text
u=x^(e/n)v^(g/n),
t=x^(-e/n)v^(-ae/n).                                  (9)
```

Because `gcd(e,g)=1`, all four exponents in `(9)` are integers exactly when
`|n|=1`.  Equivalently, both target coordinates belong to `C(x,y)` exactly
on the residue-coordinate locus `(8)`.  For `|n|>1`, `(2)` remains a valid
symplectic parametrization in `C(u,r)`, but it is a finite torus cover rather
than a rational inverse chart.

## 4. The sign gives the first polynomiality fork

Assume `|n|=1`.  If `n=1`, both exponents in the formula for `u` are
nonnegative:

```text
u=x^e v^g in C[x,y],
t=x^(-e)v^(-ae).                                      (10)
```

If `n=-1`, the signs reverse:

```text
t=x^e v^(ae) in C[x,y],
u=x^(-e)v^(-g).                                       (11)
```

The negative `x` powers in `(10)--(11)` alone do not prove nonmembership,
because `y=L/x` is polynomial.  The exact source-ring intersection below
supplies the missing boundary divisibility and proves that the other target
is genuinely nonpolynomial.

## 5. Exact initial-ring intersection

Let `A=C[x,y]` and `B=C[u,t]`, regarded inside `C(x,v)` on the gate
`|n|=1`.  Since `y=L/x`, there is an exact Laurent description

```text
A = direct_sum_(m>=0) x^m C[v]
    direct_sum_(q>=1) x^(-q)L^q C[v].                  (12)
```

To prove `(12)`, write `x^i y^j=x^(i-j)L^j`.  At a fixed exponent
`m=i-j`, the allowed `j` begin at `max(0,-m)`, and the powers
`L^j` span `L^max(0,-m) C[v]`.  The converse is termwise.  Thus a negative
`x` coefficient must both be regular at `v=0` and vanish to the prescribed
order at `L=0`.

### 5.1 The full rational initial-coefficient ring

There is a sign-independent description before intersecting with target
polynomials.  From `(9)`,

```text
C(u,t)=C(v,x^e).                                        (13)
```

Thus `C(x,y)/C(u,t)` is the cyclic degree-`e` extension on which

```text
zeta in mu_e:      x |-> zeta x,      y |-> zeta^(-1)y,
                   v |-> v.                            (14)
```

The source-polynomial elements which are rational target functions are
exactly the invariants:

```text
A intersect C(u,t)
 =A^(mu_e)
 =C[v,X=x^e,Y=y^e]
 =C[V,X,Y]/(XY-(1+cV)^e).                              (15)
```

Indeed invariant monomials have exponent difference divisible by `e`; after
removing powers of `xy=L`, they are generated by `x^e,y^e,L`, equivalently
by `X,Y,v`.  This is the complete polynomial initial-coefficient ring for
rational target functions.  Its orientation is invisible until one asks
which target *polynomials* land in it.

For `e=3`, putting `z=L` identifies `(15)` with
`C[X,Y,z]/(XY-z^3)`, exactly the `m=2` normal-crossing cubic algebra of
THM-2690.  Its class group is `Z/3` and it has no binary Kummer carrier on
the regular locus.  This is a consistency bridge only: no result here
identifies the terminal initial ring with a physical cyclic-cubic resolvent
normalization or supplies the common three-view atom needed for that use.

### 5.2 Positive orientation

If `n=1`, then `g=ae+1`.  Grade a polynomial in `B` by the difference
`deg_u-deg_t`.  Formula `(12)` gives

```text
B intersect A
 = direct_sum_(q>=0) u^q C[v]
   direct_sum_(q>=1) G_+^q C[v]
 = C[u,v,G_+],                                          (16)

G_+=u^(-1)v^g L^e
   =u^(g-1)t^g(1+cut)^e
   =y^e.                                                (17)
```

Indeed a negative piece `u^(-q)f(v)` pulls back to
`x^(-eq)v^(-gq)f(v)`.  By `(12)` it lies in `A` exactly when
`v^(gq)L^(eq)` divides `f`.  Those are precisely the multiples of
`G_+^q`.  The presentation has the single relation

```text
u G_+ = v^g(1+cv)^e.                                   (18)
```

After localizing at `u`, `(18)` embeds in `C[u,u^(-1),v]`, so there is no
second relation.  The missing coordinate `t=u^(-1)v` lacks the required
`v^g L^e` factor and is not in `A`.

### 5.3 Negative orientation

If `n=-1`, then `g=ae-1` and `ae>=2`.  Grading instead by
`deg_t-deg_u` gives

```text
B intersect A
 = direct_sum_(q>=0) t^q C[v]
   direct_sum_(q>=1) G_-^q C[v]
 = C[t,v,G_-],                                          (19)

G_-=t^(-1)v^(ae)L^e
   =u^(ae)t^(ae-1)(1+cut)^e
   =y^e,                                                (20)

t G_- = v^(ae)(1+cv)^e.                                (21)
```

Now a negative piece `t^(-q)f(v)` belongs to `A` exactly when
`v^(aeq)L^(eq)` divides `f`.  Hence `(19)--(21)` are exact and the missing
coordinate `u=t^(-1)v` is not in `A`.

Together the two cases prove

```text
exactly one of u,t is in C[x,y],
and B intersect A is an explicit hypersurface ring.     (22)
```

No polynomial target automorphism can repair this failure.  If
`Phi in Aut_C C[u,t]` had both components in `A`, applying its polynomial
inverse would put both `u` and `t` in `A`, contradicting `(22)`.  This
corollary concerns target automorphisms; it says nothing about arbitrary
noninvertible Keller retargetings.

The two divisibilities at `v=0` and `L=0` are the promised initial-ring
sidecar.  They retain strictly more than the depth, gcd, and residue-field
ledgers.

## 6. Complete normalized arithmetic atlas for `C3`

Set the normalized tame index `e=3`.  Terminal primitivity requires
`gcd(g,3)=1`.  Solving `g-3a=+/-1` gives exactly

```text
g=3k+1:   a=k,       n= 1,       u polynomial;
g=3k+2:   a=k+1,     n=-1,       t polynomial.          (23)
```

There is no `g=3k` row: such a row would violate terminal primitivity.
In an ambient presentation with exponents multiplied by three it has
effective normalized index one, not three.  This is the sharp normalization
boundary.

For example,

```text
(g,a,n)=(2,1,-1):
  u=9/(x^3T^2),                 t=-x^3T^3/27;

(g,a,n)=(4,1,1):
  u=x^3T^4/81,                 t=27/(x^3T^3).           (24)
```

The two cells have the same tame index and primitive terminal lattice but
opposite polynomial targets.

## 7. Strict-depth compositions are source-coordinate gauges

Write the standard terminal packet as

```text
x=R^(-1),                       W=MR.                   (25)
```

Both triangular changes

```text
y_(2)=x+W,
y_(1,1)=x+1+W                                             (26)
```

are polynomial source automorphisms: each preserves `C[x,W]` and
`dx wedge dW`.  Yet their first toric ratios are

```text
y_(2)/x   =1+R^2 M,              depths (2g,e),
y_(1,1)/x =1+R(1+RM),            depths (g,g,e).        (27)
```

The successive normalized deviations in the second row are `1+RM` and
`M`; hence the displayed depth compositions are exact.  The base chart
`(x,W)` is the one-stage terminal packet of depth `e`.  Therefore the
differential depth composition is not invariant under polynomial triangular
changes of source coordinates, while the embedded rings in Section 5 are.

THM-3080's two multistage controls are precisely the two `C3` instances:

```text
g=2, a=1, n=-1:       (2g,e)=(4,3),
g=4, a=1, n= 1:       (g,g,e)=(4,4,3).                 (28)
```

In the first packet `W=y_0-x`; in the second `W=y_0-x-1`, agreeing with
the explicit shears in THM-3080.  Several strict differential keys can
therefore be a coordinate presentation of the same terminal initial module.
A genuinely surviving obstruction must use that module (or an equivalent
boundary sidecar), not the depth word alone.

## 7b. Polynomiality debt is an unbounded response staircase

Let `E=B intersect A`, viewed only as an additive subspace of `B`.  In the
positive orientation, `(16)` gives `E=C[u,v,G_+]`.  Because `uE` is contained
in `E`, multiplication by `u` induces a linear operator `M_u` on the additive
quotient `B/E`.  For every `q>=1` it has the exact response string

```text
[t^q], [u t^q], ..., [u^(q-1)t^q] != 0,
M_u^q[t^q]=[u^q t^q]=[v^q]=0.                         (29)
```

Indeed, for `0<=k<q` put `r=q-k`.  Then

```text
u^k t^q=u^(-r)v^q.                                    (30)
```

The exact membership test in Section 5.2 would require
`v^(gr)L^(er)` to divide `v^q`.  This is impossible: `r,e>0` and
`L=1+cv` is coprime to `v`.  At `k=q`, however, the monomial is `v^q` and
belongs to `E`.  Thus the least killing time is exactly `q`.

In the negative orientation, `(19)` gives `E=C[t,v,G_-]`.  Multiplication by
`t` acts on `B/E`, and the symmetric calculation gives

```text
[u^q], [t u^q], ..., [t^(q-1)u^q] != 0,
M_t^q[u^q]=[t^q u^q]=[v^q]=0.                         (31)
```

Every individual polynomial class is eventually killed, since this holds
monomial by monomial, but the exact response lengths are unbounded.  Hence
the Laurent effectivity debt is an unbounded torsion staircase even though
the residue field, cyclic torsor, and class-group ledgers are finite.  This
is precisely THM-3397's all-power denominator filtration, recast as an
operation-response module, and is the effectivity analogue of THM-3466: a
bounded response bank cannot certify polynomiality for every Laurent power.
It is not a new
obstruction to an arbitrary Keller map and supplies no FC moment nullity.

## 8. Information contract and boundary

```text
source:     terminal Laurent key chart plus the polynomial source ring
target:     residue coordinate and the two decoded target functions
map:        (R,M) -> (x=R^-1,y=MR), then T=xy-1
preserved:  exact inverse Keller form, normalized tame index e, terminal depth,
            primitive gcd lattice and residue field
destroyed:  under field-only decoding, the oriented polynomial exponent cone
sidecar:    A intersect C(u,t), B intersect A, and the v=0/L=0 divisibilities
test:       the opposite C3 cells (g,a)=(2,1) and (4,1).  (32)
```

The theorem excludes the complete family `(1)--(2)` and the two displayed
multistage controls after their source shears.  It does **not** prove that an
arbitrary terminal key module has this form, that every Laurent key is a
polynomial approximate root, or that a general Jelonek component is a target
coordinate line.  It excludes no remaining nonmonomial `C3`, `A4/S4`, `G1`,
`JC(2)`, or `DC(2)` branch by itself.

## 9. Exact verification

The standard-library companion implements Laurent polynomials over `Q`,
checks `(3)--(9)` on `34,987` normalized coprime parameter cells, verifies
`363` rational-decoding cells, checks the cyclic invariant hypersurface and
both oriented target-polynomial generators in every such cell, proves the
normalized `C3` atlas through `g=24`, and reconstructs both multistage source
shears.  An independent replay found ordinary and optimized runs
byte-identical to the stored transcript.

Reproduce with

```text
python3 04-computation/jc_terminal_monomial_cone_polynomiality_fork_thm3383.py
python3 -O 04-computation/jc_terminal_monomial_cone_polynomiality_fork_thm3383.py
```

Artifact and semantic hashes are pinned in the frontmatter.

**QED.**
