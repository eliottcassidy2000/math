---
id: THM-2671
title: "Degree-twenty-two complete support-three weighted Hensel closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the inherited
  genuine nonsplit degree-twenty-two branch on the open first-flux chart,
  all ten strata with exactly three active coefficients among B,C,D,E,W are
  empty.  THM-2617 and THM-2636 supply BDW and BCD.  For each of the other
  eight triples, a coefficient root supplies a signed scale t and a monic
  quintic eliminant with the same irreducible squarefree fixed section L_5.
  Its five roots and ten unordered root pairs exhaust all possible line and
  quadratic factors.  At the first order beyond the specialized t-degree,
  exact Hensel lifting gives five equations whose coefficient matrix has
  full column rank in both fields, uniformly contradicting the active
  parameter monomials.  The retained square T^2=Z has five odd fixed places,
  so its connected double cover has genus at least two.  Two uniform
  quintic boundary eliminants close y=0 for B nonzero and B=0.  Hence every
  remaining degree-twenty-two trajectory in this inherited genuine nonsplit
  branch has coefficient support at least four.  Split/even descent, support four or five, integral raising, JC(2),
  and DC(2) remain open.
source: root-long-frontiers-2026-07-28-complete-support-three
depends_on:
  - THM-2411-degree-twenty-two-first-flux-pole-divisor-square-class-reduction
  - THM-2480-degree-twenty-two-BC-plane-hensel-ramification-closure
  - THM-2617-degree-twenty-two-BDW-triple-fixed-section-and-last-quadratic-type
  - THM-2636-degree-twenty-two-BCD-triple-spectral-square-Kummer-closure
related:
  - THM-2446-twojet-zgraded-jacobian-decomposition-and-cone-system
script: 04-computation/jc2_degree22_complete_support_three_hensel_thm2671.py
output: 05-knowledge/results/jc2_degree22_complete_support_three_hensel_thm2671.out
script_sha256: 6245dd4cc85d0a70bdbc8e56a0511ffad7889b6274130aed759e5729f92472e6
output_sha256: 189cc0b80d282363e899b5f4112910a09ace904579eda130313afec43b3d8657
hash_basis: working-tree bytes (LF)
---

# THM-2671 -- every degree-twenty-two support-three stratum is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2411 leaves five weighted constant coefficients

```text
(B,C,D,E,W),                    weights=(2,3,4,5,6),       (1)
```

on the genuine nonsplit degree-twenty-two branch.  The ten support-two
planes are empty, while THM-2617 and THM-2636 close the support-three tori
`BDW` and `BCD`.  This theorem closes the other eight triples by one common
fixed-section mechanism.  The result is a genuine frontier reduction:

```text
every surviving trajectory in this inherited branch has support at least four. (2)
```

It is not a closure of degree twenty two or of the planar Jacobian problem.

## 1. One normalized flux system, three scale choices

Retain THM-2411's coordinates

```text
y=11s,                    u=d_0T,                    Z=T^2,
mathcal A=616B-1089u+63y^2 != 0.                       (3)
```

First take `y!=0` and put

```text
v=u/y^2,                  zeta=Z/y^3,
(b_2,b_3,b_4,b_5,b_6)=(B/y^2,C/y^3,D/y^4,E/y^5,W/y^6). (4)
```

The two normalized fluxes are exactly the pre-scale equations `F_1,F_2` in
the companion.  Eliminating `zeta` before choosing a scale gives a primitive
sixty-term polynomial

```text
R_gen(b_2,b_3,b_4,b_5,b_6;v),                           (5)
```

with raw integer content `28,344,976`.  It specializes to THM-2636's
universal eliminant under

```text
b_i=a_i t^i.                                             (6)
```

Over the algebraically closed constant field, every required coefficient
root exists.  The exhaustive scale atlas is:

| active triples | choice of `rho` | normalized nonzero parameters |
|---|---|---|
| `BCD,BCE,BCW,BDE,BDW,BEW` | `rho^2=B` | the other two coefficients divided by `rho^i` |
| `CDE,CDW,CEW` | `rho^3=C` | `D/rho^4,E/rho^5,W/rho^6` as present |
| `DEW` | `rho^4=D` | `E/rho^5,W/rho^6` |

In every chart put `t=rho/y`.  Thus the scale supplier becomes `t^2`, `t^3`,
or `t^4`, and each remaining active coefficient is a nonzero constant times
the corresponding `t^i`.  These ten charts are disjoint and exhaustive; no
choice of a square, cube, or fourth root is asserted to be canonical.

## 2. The fixed root and root-pair atlas is uniform

For each specialization, primitive-part normalization gives a polynomial

```text
P(t,v) in C[t,v],                 deg_v P=5,              (7)
```

whose `v^5` coefficient is a nonzero constant.  It is therefore monic after
constant rescaling.  Every chart has the same `t=0` section, namely the monic
form of

```text
L_5(v)=
 155624547606v^5+3215383215v^4-1700698560v^3
 +58124770v^2-855470v+2583.                              (8)
```

This quintic is irreducible and squarefree over `Q`.

If `P=QS` is reducible over `C(t)[v]`, monicity and the integral closure of
`C[t]` (equivalently Gauss integrality here) let both factors be chosen monic
in `v` inside `C[t,v]`.  Their `v`-degrees remain positive at `t=0`
and add to five.  After swapping them, the smaller fixed factor has degree
one or two.  A line selects one of the five roots of `(8)`; a quadratic
selects one of its ten unordered root pairs.  As in THM-2636, these form one
irreducible degree-five root field `A_1` and one irreducible degree-ten pair
field `A_2`.  They cover every complex root and every unordered pair under
their embeddings, so no factor component is omitted.

## 3. The eight terminal-rank obstructions

Write

```text
P=sum_n r_n(v)t^n,
Q=sum_n q_n(v)t^n,              S=sum_n s_n(v)t^n.        (9)
```

In `A_1` or `A_2`, the fixed factors are coprime.  Hensel uniqueness gives
the exact recursion

```text
f_n=r_n-sum_(i=1)^(n-1)q_i s_(n-i),
q_n=rem_(q_0)(f_n (s_0 mod q_0)^(-1)),
s_n=(f_n-q_n s_0)/q_0.                                  (10)
```

Let `D_t=deg_t P` and `N=D_t+1`.  In an actual polynomial factorization,
degree additivity gives `deg_t Q+deg_t S=D_t`; hence

```text
q_N=s_N=0.                                               (11)
```

The companion evaluates `(10)` exactly through `N`, checks every division
and every reconstructed product coefficient, and writes the five scalar
equations in `(11)` as a matrix times their complete parameter-monomial
vector.  The result is:

| stratum | parameters | terms | `(D_t,N)` | terminal monomials | rank in `(A_1,A_2)` |
|---|---|---:|---:|---|---:|
| `BCE` | `(C,E)` | 41 | `(10,11)` | `E,C,C^2E,C^3` | `(4,4)` |
| `BCW` | `(C,W)` | 38 | `(10,11)` | `C,CW,C^3` | `(3,3)` |
| `BDE` | `(D,E)` | 36 | `(10,11)` | `E,DE` | `(2,2)` |
| `BEW` | `(E,W)` | 31 | `(10,11)` | `E,EW` | `(2,2)` |
| `CDE` | `(D,E)` after `C=1` | 25 | `(10,11)` | `E,D^2` | `(2,2)` |
| `CDW` | `(D,W)` after `C=1` | 22 | `(8,9)` | `1,W` | `(2,2)` |
| `CEW` | `(E,W)` after `C=1` | 21 | `(10,11)` | `E,EW` | `(2,2)` |
| `DEW` | `(E,W)` after `D=1` | 19 | `(10,11)` | `EW` | `(1,1)` |

Every listed matrix has full column rank at every embedding.  In the root
field, the first nonzero maximal minor always has numerator degree four and
five terms; in the pair field it has degree nine and ten terms.  Each
numerator is coprime to its field modulus.  The numbers of nonzero maximal
minors in both fields are respectively

```text
BCE:2, BCW:7, BDE:9, BEW:9, CDE:9, CDW:9, CEW:9, DEW:5. (12)
```

Consequently `(11)` forces every coordinate in the corresponding monomial
vector to vanish.  This is impossible on the physical torus: for example
`CDE` would force `(E,D^2)=0`, `CDW` would force `(1,W)=0`, and `DEW` would
force `EW=0`.  The other five rows are immediate from their displayed
nonzero active factors.  Thus none of the eight eliminants has a line or a
quadratic factor.  Section 2 exhausts the smaller side of every possible
factorization, so all eight are absolutely irreducible.

Together with THM-2636, this proves the required signed-`t` absolute
irreducibility on the eight new tori and on `BCD`.  The remaining `BDW`
torus is closed separately by THM-2617 in its even coordinate
`p=B/y^2`; no irreducibility transfer from `R(p,v)` to `R(t^2,v)` is used.

## 4. The same five odd places close every new nonconstant signed-scale chart

At every root of `(8)`, squarefreeness makes `t` a uniformizer.  The fixed
first flux is independent of the chosen scale:

```text
(83853-1449459v)zeta
 +(3689532v^2-101640v+252)=0.                            (13)
```

Both coefficient polynomials are coprime to `L_5`; hence `zeta` is a unit at
all five fixed points.  From `t=rho/y` and `zeta=Z/y^3`, in every chart,

```text
H=rho^3 zeta/t^3=Z=T^2.                                 (14)
```

Thus, on each of the eight new charts (and on `BCD` by THM-2636),
`ord(H)=-3` at each fixed point.  The Kummer cover `X^2=H` is connected
and has at least five visible odd branch places; parity raises the total to
at least six.  Riemann--Hurwitz gives genus at least two, independently of
the base-curve genus.

If `t` is nonconstant, a trajectory would embed this connected cover into
`C(x)`, producing a nonconstant map from `P^1` to a genus-at-least-two curve,
which is impossible.  If `t` is constant, then `y` is constant; monicity of
`P(t,v)` makes `v`, then the open first flux makes `zeta`, and `(14)` makes
`Z,T,q` constant, contradicting the genuine nonsplit deck.  Hence the
`y!=0` branch is empty in the eight new charts; THM-2636 closes `BCD` and
THM-2617 separately closes the even-coordinate `BDW` chart.

## 5. Two uniform quintics close `y=0`

Treat `y=0` before defining `t`.  If `B!=0`, normalize `B=1` and retain all
four remaining coefficients, setting the inactive ones to zero only after
elimination.  The primitive resultant of the first two fluxes is

```text
-2264031u^5+5305608u^4-(3763584D+3829056)u^3
 +(-1267200C^2+4257792D-2509056W+878080)u^2
 +(1433600C^2-1204224D+2838528W)u
 -2293760CE+3244032E^2-802816W.                          (15)
```

Its raw content is `1,104,726,788,605,792` and its constant leading
coefficient is `-2,264,031`.

If `B=0`, the open chart first gives

```text
-1089u!=0.                                               (16)
```

Only then solve the first flux for `Z`.  The primitive resultant is

```text
-22869u^5-38016Du^3-(12800C^2+25344W)u^2+32768E^2,      (17)
```

with raw content `109,367,952,071,973,408` and constant leading coefficient
`-22,869`.  Equations `(15)` and `(17)` make `u` algebraic over the
algebraically closed constant field, hence constant.  The open first flux
then makes `Z`, and therefore `T,q`, constant, again contradicting the deck.
This closes `y=0` uniformly, including the cases where a lower coefficient
or constant term of the displayed quintic vanishes.

## 6. Consequence, reproduction, and boundary

The ten support-three strata are

```text
BCD, BCE, BCW, BDE, BDW, BEW, CDE, CDW, CEW, DEW.        (18)
```

Sections 1--5 and THM-2617/2636 empty all ten.  Combined with the existing
axis and plane closures, any degree-twenty-two trajectory in the inherited
genuine nonsplit branch has at least four active coefficients.

Run

```bash
python3 04-computation/jc2_degree22_complete_support_three_hensel_thm2671.py
python3 -O 04-computation/jc2_degree22_complete_support_three_hensel_thm2671.py
```

Both executions byte-match the declared output.  The companion pins the
audited THM-2636 engine by SHA-256, independently reconstructs the pre-scale
eliminant and all eight specializations, performs both root- and pair-field
Hensel recursions, checks every terminal support, rank, maximal-minor
coprimality, and product coefficient, and derives both boundary quintics.

An independent hostile audit reconstructed the scale atlas, factor
exhaustion, terminal-degree implication, eight torus contradictions, fixed
unit and Kummer argument, and the `B=0` order of division.  It independently
derived both boundary resultants and replayed normal and optimized modes
against the frozen transcript and declared hashes.

The theorem does not address support four or five, split/even short-edge
descent, integral `2`-adic raising, another degree, `JC(2)`, or `DC(2)`.

QED.
