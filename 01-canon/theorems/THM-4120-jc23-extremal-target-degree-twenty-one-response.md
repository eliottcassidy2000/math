---
id: THM-4120
title: "JC(2) reduced 2:3 extremal target, unique degree-21 response, and DE principal kernel"
status: >
  PROVED RELATIVE TO THM-3992/4053/4103 + CITED SHIODA--TATE +
  VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT on the smooth nonresonant
  theta-only maximum-weight-eight survivor. Its target pencil is the
  extremal rational elliptic surface II*+I1+I1 with Mordell--Weil group
  {O}. All seven source-infinity punctures therefore map to O, the generic
  fibre degree is exactly 21, and the five labelled responses of THM-4103
  collapse to one. The DE edge is an evaluation quotient with principal
  kernel (x^6t^7+gamma/theta); Keller cancellation locks the two output
  multiplicities and raises the necessary total-degree floors to
  deg(A)>=28 and deg(C)>=31. THM-4130 later excludes this smooth seam;
  THM-4134 excludes the Delta_V full-boundary degrees 20/19 but leaves
  horizontal-BC degrees 16/15. JC(2), the other walls, other cells, and
  maximum residual weight at least nine remain OPEN.
source: planar-jacobian-squeeze / 2026-08-25
audit: >
  PASS. The primary SymPy certificate and independent standard-library
  sparse-polynomial audit agree on the target discriminant and fibre ledger,
  Shioda--Tate rank/discriminant bookkeeping, boundary residue fields,
  degree response, pullback pole degrees, DE layer parametrizations, general
  fixed-layer bracket, re-entry coefficients, square/cube response, and
  degree floors. Both normal and optimized executions byte-match the frozen
  outputs. An independent agent rederived the DE kernel, the phi=0 forced
  second-order sidecar, the simultaneous target coefficient h(T)=-7, and
  the both-baseline/both-higher dichotomy.
depends_on:
  - THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual
  - THM-4053-jc2-live-max-eight-trichotomy-and-eisenstein-survivor
  - THM-4103-jc23-theta-boundary-ramification-and-degree-response
related:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-3996-two-three-companion-node-address-balance-and-jelonek-alternative
  - THM-4130-theta-only-extremal-seam-critical-monodromy-obstruction
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
external:
  - "T. Shioda, On the Mordell--Weil lattices, Comment. Math. Univ. St. Pauli 39 (1990), 211--240, DOI:10.14992/00009994."
  - "R. Miranda and U. Persson, On extremal rational elliptic surfaces, Math. Z. 193 (1986), 537--558, Theorem 4.1 and Table 5.2, DOI:10.1007/BF01160474."
script: 04-computation/jc23_extremal_target_de_response_thm4120.py
output: 05-knowledge/results/jc23_extremal_target_de_response_thm4120.out
independent_audit_script: 04-computation/jc23_extremal_target_de_response_thm4120_independent_audit.py
independent_audit_output: 05-knowledge/results/jc23_extremal_target_de_response_thm4120_independent_audit.out
script_sha256: a03716cc84df1391f7ba8e15ea362133405edc93b127e650fa68cbb9c51327c5
output_sha256: a81cb731233ac7367c624fce3429ba5db0d5a292d5aa4e0034a5b5e1ea5eab32
independent_audit_script_sha256: 1f60a83be3a93778de19bc755cff266b786912686b1b7e591e19e48911bfaa01
independent_audit_output_sha256: 8ad5b294af06770b53348e300dd8090f832896ee1c776710a8540eba6a239b3a
hash_basis: raw LF bytes
---

# THM-4120 -- the extremal target collapses the theta response

**PROVED RELATIVE TO THM-3992/4053/4103 + CITED SHIODA--TATE +
VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT on one reduced seam; JC(2)
OPEN.** Work over an algebraically closed field `k` of characteristic zero.
Retain every hypothesis and notation of THM-4103's smooth nonresonant
theta-only survivor, put `K=k(q)`, and let

```text
varphi_q:C_q -> E_q
```

be its generic-fibre map over `K`. Then

```text
deg(varphi_q)=21,                                         (1)
```

and every one of the seven source-infinity punctures maps to the target
origin `O`. More precisely,

```text
varphi_q^*(O)=P_AB+2(P_BC,+ + P_BC,-)+3P_CD+7P_DE
              +3(P_EF,+ + P_EF,-).                      (2)
```

Consequently the pole divisors of the polynomial coordinates satisfy

```text
deg (A)_infinity=42,             deg (C)_infinity=63.    (3)
```

At the index-seven `DE` puncture, define

```text
wt_DE(x^m t^n)=7m-6n,       w=x^6t^7=s^6t,
T=-gamma/theta.                                           (4)
```

Every fixed `DE`-weight layer is a monomial times a polynomial in `w`, and
restriction to the edge is evaluation `w -> T`, with kernel `(w-T)`. The
Keller equation couples the cancellation orders of the top `A` and `C`
layers. Together with `(1)--(3)` and THM-3992's exact Laurent depths this
forces

```text
deg A>=28,                         deg C>=31.             (5)
```

These are total polynomial degrees. They sharpen THM-4103's prior floors
`15,16`; they do not construct or exclude the survivor.

## 1. The target pencil has only its zero section

THM-4103's target curve is

```text
E_q: V^2=U^3-(3a^2/4)U+q-a^3/4,              a!=0.      (6)
```

Its short-Weierstrass invariants are

```text
c4=36a^2,
Delta=-432q(q-a^3/2).                                  (7)
```

Thus the finite fibres at `q=0,a^3/2` have simple discriminant zero and
nonzero `c4`; both are type `I1`. At infinity put `u=1/q` and make the
integral change `U=u^-2 X,V=u^-3 Y`. The local coefficients become

```text
a4'=-(3a^2/4)u^4,             a6'=u^5-(a^3/4)u^6,
(v_u(c4),v_u(c6),v_u(Delta))=(4,5,10).                  (8)
```

The model is minimal and `(8)` is type `II*`. Its fundamental line bundle
is `O_P1(1)`, so the relatively minimal surface is rational and has
Neron--Severi rank ten. The section and fibre span the odd unimodular block

```text
[[-1,1],[1,0]],                                         (9)
```

while the components of `II*` not meeting the section span `E8`. These
blocks already have rank ten and absolute discriminant one. Shioda--Tate
therefore gives Mordell--Weil rank zero, and the discriminant/index formula
gives trivial torsion:

```text
E_q(K)={O}.                                             (10)
```

This is also the `X_211` surface in Miranda--Persson's classification: their
`II*+I1+I1` row has section group of order one. Notice that `(9)` is not the
even hyperbolic plane; only its rank and unimodularity are used.

## 2. Five rational punctures force seventeen sheets

The complete Newton boundary in THM-4103 has edges

```text
AB, BC, CD, DE, EF
```

of lengths `1,2,1,1,2`. The linear edge points `AB,CD,DE` are `K`-rational.
The two `EF` points are the distinct roots of the constant-field polynomial

```text
theta+phi Z+epsilon Z^2,                               (11)
```

whose discriminant is `Delta_V!=0`; since `k` is algebraically closed, both
are also `K`-rational. Their ramification indices sum to

```text
1+3+7+3+3=17.                                          (12)
```

The map `varphi_q` is defined over `K`. Equations `(10)--(12)` force all
five rational punctures to `O`.

The two `BC` punctures instead form the single quadratic closed point

```text
kappa Y^2=q+gamma.                                     (13)
```

It is irreducible over `K`, because `q+gamma` has odd valuation at
`q=-gamma`. Galois conjugacy makes its two index-two points respond
together. Hence the only degrees left by target rationality are

```text
n=17 or 21.                                            (14)
```

THM-4053 proves that `n` is an Eisenstein norm `r^2-rs+s^2`. The prime
`17=2 mod 3` occurs to odd exponent, so `17` is not such a norm. This proves
`(1)`, and the equality of the total index with `21` proves `(2)`. Since
`U,V` have pole orders `2,3` at `O`, `(3)` follows.

The norm gate in this last step is genuinely load-bearing. After the
quadratic base change

```text
q=a^3/2+r^2,                                           (15)
```

the target acquires the point `(U,V)=(a/2,r)`. Thus `(10)` cannot be
silently extended from `K` to the `BC` residue field.

## 3. The DE edge is a principal evaluation kernel

Along `DE`, THM-4103 gives

```text
z=1/s,        t~Tz^6,        x~T^-1z^-7.               (16)
```

Solutions of `7m-6n=D` differ by `(m,n)->(m+6,n+7)`.
Therefore a fixed weight layer has the form

```text
A_D=x^a t^b f(w),          C_E=x^c t^d g(w),            (17)
```

for polynomials `f,g`. Equation `(16)` restricts `(17)` by evaluating
`w` at `T`; its kernel is exactly `(w-T)`.

This kernel has a known re-entry cost. Expanding the exact source equation
`t(1-q^-1 H)=gamma q^-1 s^2` gives

```text
w=T+(gamma phi/theta^2)z
   +(gamma(epsilon+kappa)/theta^2-gamma phi^2/theta^3)z^2+... . (18)
```

Hence

```text
ord_z(w-T)=1                         if phi!=0.          (19)
```

If `phi=0`, the inherited theta-only row is essential:

```text
epsilon+kappa=-14336 gamma/(135 A5^3)!=0,
A5=a^5,                                                    (20)
```

so `(18)` instead gives `ord_z(w-T)=2`. The assertion would be false for
arbitrary coefficients satisfying only `phi=0`; `(20)` is load-bearing.

For general layers `(17)`, direct differentiation gives

```text
J(A_D,C_E)=x^(a+c-1)t^(b+d-1)
 [(ad-bc)fg+D w f g'-E w f'g].                          (21)
```

Let `D_A,D_C` be the maximum `DE` weights occurring in `A,C`. The exact
extreme terms `gamma^2x^2,gamma^3x^3` give

```text
D_A>=14,                         D_C>=21.                (22)
```

The top weight of the constant Keller bracket is positive, so its
coefficient `(21)` vanishes. If the top `f,g` vanish at `T` to orders
`r,s`, the term of order `r+s-1` in `(21)` is a nonzero scalar times

```text
D_A s-D_C r.                                           (23)
```

Thus

```text
D_A s=D_C r.                                           (24)
```

The two outputs cannot cancel independently.

## 4. The degree-floor dichotomy

By `(1)`, the `DE` point maps to `O`. THM-4103's target parameter therefore
fixes the actual leading terms

```text
A~49theta^2 z^-14,             C~343theta^3 z^-21.     (25)
```

There are two cases.

### 4.1 Both maximum weights are minimal

If `(D_A,D_C)=(14,21)`, THM-3992's depth bounds and the two exceptional
extremes give

```text
A_14=x^2 f(w),       f(0)=gamma^2,
C_21=x^3 g(w),       g(0)=gamma^3.                     (26)
```

Equation `(21)` becomes

```text
7x^10t^6(2fg'-3f'g)=0.                                 (27)
```

Unique factorization in `k[w]` and the constants in `(26)` give

```text
f=gamma^2 h^2,        g=gamma^3 h^3,        h(0)=1.    (28)
```

At `w=T`, equations `(16),(25),(28)` read

```text
theta^2 h(T)^2=49theta^2,
-theta^3 h(T)^3=343theta^3.                            (29)
```

The two equations jointly force `h(T)=-7`. Since `h(0)=1` and `T!=0`,
`h` is nonconstant. Its square and cube in `(28)` force monomials of total
degrees at least

```text
2+2(6+7)=28,                  3+3(6+7)=42.             (30)
```

The first equation in `(29)` alone would leave both signs; the simultaneous
`A,C` response is essential.

### 4.2 Both maximum weights are higher

Suppose `D_A>14`. Its top layer must vanish at `T`, or it would give a pole
deeper than `(25)`. Equations `(24)--(25)` then imply `D_C>21`: if
`D_C=21`, its top layer would also vanish at `T`, after which neither that
layer nor any lower one could supply the required pole order `-21`. The
same argument with `A,C` reversed proves

```text
D_A>14 iff D_C>21.                                     (31)
```

Apart from `x^2`, every monomial of `A` obeys `m<=n+1`. On a layer
`D_A>14`, the allowed ladder begins

```text
(m,n)=(D_A-6,D_A-7), (D_A,D_A), ... .                  (32)
```

Vanishing at the nonzero point `T` requires at least two ladder monomials,
so

```text
deg A>=2D_A>=30.                                       (33)
```

Similarly, the `C` depth `m<=n+2` makes a layer `D_C>21` begin

```text
(D_C-12,D_C-14), (D_C-6,D_C-7), ...,
deg C>=2D_C-13>=31.                                    (34)
```

Taking the coordinatewise minimum of `(30)` and `(33)--(34)` proves `(5)`.

## 5. Exact audit and failure boundaries

The primary certificate verifies the target invariants, fibre valuations,
rank/discriminant ledger, boundary weights, Eisenstein filter, quadratic
hostile, edge expansion, general bracket `(21)`, Diophantine layers, target
coefficients, and both degree cases symbolically. A separate
standard-library audit reconstructs the same ledgers with sparse polynomial
differentiation and exact fractions. Normal and optimized executions match
the frozen outputs byte for byte.

Factors `(w-T)` can re-enter across weight layers: at re-entry cost one a
weight-fifteen layer may contribute at effective weight fourteen, and at
cost two a weight-sixteen layer may do the same. Therefore this theorem does
not claim that the raw weight-fourteen coefficient alone equals the target
response. Equations `(18)--(24)` are the required sidecar.

Finally, the theorem has a narrow entry contract. It concerns only the
smooth nonresonant theta-only `M=8`, `b=d=0`, oriented reduced `(2,3)` cell
inside the completion used by THM-3992. It does not apply automatically to
the three collision walls, `(b,d)!=(0,0)`, residual weight at least nine,
other reduced depth pairs, arbitrary completions, or all of `JC(2)`.

THM-4130 later bypasses this coefficient program and excludes the smooth
seam: its exact 20-point critical ledger forces a three-cycle commutator,
contradicting the inherited boundary packet. THM-4134 later normalizes
`Delta_V=0`, excludes degrees `20,19`, and leaves horizontal-BC degrees
`16,15`; the other two collision walls remain open.

**QED.**
