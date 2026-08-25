---
id: THM-4089
title: "Hybrid p-adic zeta margin optimization and next-case obstruction"
status: >
  PROVED ELEMENTARY OPTIMIZATION OF THE DISPLAYED EXTERNAL MARGIN +
  VERIFIED-EXACT + INDEPENDENTLY CALCULUS/NUMERIC-AUDITED. The cutoff cost is
  globally convex and has an explicit
  unique minimizer; the analytic-radius contribution is globally strictly
  concave and has a unique critical maximizer. The 22 displayed witnesses
  have positive formula margins, while exact derivative brackets and tangent
  bounds prove negative global maxima for the immediate next cases (2,31),
  (3,13), (5,7), and (7,5). This is a sharp obstruction to extending the
  manuscript by parameter retuning alone. The first-chamber stationary point
  is interior exactly for p<11; at p=13,s=3 even an idealized one-power Hasse
  cost has negative margin for every continuation radius. The manuscript's
  22 p-adic-zeta irrationality theorem remains AUTHOR-CLAIMED / UNREFEREED; its
  new geometric and adelic bridges are not certified here.
source: codex-padic-zeta-tournament-20260825
audit: >
  PASS. The frontier certificate uses exact derivative brackets and rigorous
  tangent bounds for the four next-weight cells, with a separate 100-digit
  implementation checking the optima. The boundary certificate checks 1,064
  general stationary-point identities, all five p=13 arithmetic chambers,
  the exact exponential bound, and both all-Y obstructions. The p=13 algebra,
  endpoints, quantifiers, p=11 boundary, ideal/actual costs, and collision
  optimizer were also independently checked line by line. Both scripts have
  byte-identical normal/optimized output, no assert statements, and no
  floating literals.
depends_on: []
related:
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
  - THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density
  - THM-4091-integral-coordinate-change-lcm-depth-boundary
script: 04-computation/hybrid_padic_zeta_margin_frontier_thm4089.py
output: 05-knowledge/results/hybrid_padic_zeta_margin_frontier_thm4089.out
script_sha256: 90232af6c30b71ae7a4ae2a3f777effb32d6706c5a1413b948478b910c3a139e
output_sha256: 393893b53c76580a8b891f61f22dc94e50ac2c3496cfd32e562b7b75ad6d4932
independent_script: 04-computation/padic_zeta_next_cell_margin_obstruction_thm4089.py
independent_output: 05-knowledge/results/padic_zeta_next_cell_margin_obstruction_thm4089.out
primary_source_audit: .scratch/padic_repo_audit_20260825/REPORT.md
independent_audit: .scratch/padic_margin_referee_20260825/REPORT.md
independent_script_sha256: 43d0313595ed4f8704572ce6589b0aaa9764ba2c39cb47f75131bdfd012df0b6
independent_output_sha256: ca51232dec53c195fe26157a6cc597668031f900c6ba04acdfc3a7582529066e
boundary_script: 04-computation/hybrid_padic_margin_boundary_thm4089.py
boundary_output: 05-knowledge/results/hybrid_padic_margin_boundary_thm4089.out
boundary_script_sha256: 28baa8f3e529c728a1ac4580c51a16b5bb20a07139ef215b22659be0f37e6a35
boundary_output_sha256: cd7cce306e7b551f680e0fc107212d5b1b440f2408fcc1de0a660b899ae2f60d
hash_basis: raw LF bytes
external_source: https://github.com/octonion/p-adic-zeta-irrationality/commit/b46a1770901551961710e155d775aae7c5ea39e7
---

# THM-4089 -- exact frontier of the displayed hybrid margin

**PROVED elementary optimization of the displayed formula +
VERIFIED-EXACT. External irrationality theorem: AUTHOR-CLAIMED / UNREFEREED.
The optimization was independently calculus- and
100-digit-numerically audited.**

This theorem does not import the external manuscript's claimed arithmetic
holonomy theorem. It takes only its displayed elementary margin function as a
defined real function, optimizes both free parameters globally, identifies
the first four cases which parameter retuning cannot reach, and proves the
first-chamber and formal `p=13,s=3` boundaries.

The external source is Christopher D. Long's
[`p-adic-zeta-irrationality` repository at commit `b46a177`](https://github.com/octonion/p-adic-zeta-irrationality/commit/b46a1770901551961710e155d775aae7c5ea39e7).
The frozen provenance and proof-risk ledger are in the
[local source audit](../../05-knowledge/reference/p-adic-zeta-irrationality-source-audit-20260825.md).
The repository and manuscript themselves request independent specialist
review. Their exact interval program checks the final numerical margins, not
the new fixed-bundle, neat-descent, Hasse-kernel, or torsor-Cartier proof
bridges.

## 1. The displayed optimization problem

Fix

```text
p in {2,3,5,7},       s>=3 odd,       m=s+1,
c_p=(p+1)/12.                                             (1)
```

For `1<xi<m`, put

```text
K=floor(s/xi),
I_s(xi)=(s+1/2)H_K-K xi
        +(m-(K+1)xi)_+^2/(2(K+1)),                        (2)
J_p(xi)=integral_0^xi (1-c_p x)_+ dx,
S_(p,s)(xi)=s^2 xi-(s-1)J_p(xi),                          (3)
tau_(p,s)(xi)=2(S_(p,s)(xi)+s I_s(xi))/(s+1)^2.           (4)
```

For `0<Y<1/p`, define

```text
Lambda_p(Y)=12 log(p)/(p-1)-2 pi Y,
C_p(Y)=4 pi sum_(p|c, cY<1) phi(c)/c^2
       *[acos(cY)-cY sqrt(1-c^2Y^2)],                     (5)
M_(p,s)(xi,Y)=s Lambda_p(Y)-(s+1)tau_(p,s)(xi)-C_p(Y).     (6)
```

The finite sum in `(5)` is the manuscript's displayed collision energy. A
positive value of `(6)` contradicts the manuscript's claimed rationality
inequality **if all of its prior geometric and adelic inputs are valid**.

## 2. Exact global minimizer in the cutoff variable

### Theorem 2.1 (convex cutoff cost)

The function

```text
F_(p,s)(xi)=S_(p,s)(xi)+s I_s(xi)                         (7)
```

is convex on `(1,m)`. Its unique global minimizer is

```text
xi* = 12(s(s+1)-1)/(12s^2+(s-1)(p+1)).                   (8)
```

It lies in the first smooth cell

```text
1<xi*<(s+1)/s,       K=s-1,       c_p xi*<1.             (9)
```

### Proof

Away from the finitely many breakpoints, differentiation of `(2)` gives

```text
I_s'(xi)=-K-(m-(K+1)xi)_+.                               (10)
```

At a floor change `xi=s/k`, the two one-sided derivatives agree. Within the
new cell `(10)` increases linearly until the positive part vanishes and is
then constant. Thus `I_s'` is nondecreasing. Likewise

```text
S_(p,s)'(xi)=s^2-(s-1)(1-c_p xi)_+                        (11)
```

is nondecreasing. This proves convexity.

For `p<11`, the right derivative at `xi=1` is negative. On the cell `(9)`,
equations `(10)--(11)` give

```text
F_(p,s)'(xi)
 =s^2-(s-1)(1-c_p xi)+s^2(xi-2).                         (12)
```

Its unique zero is `(8)`. Direct cross-multiplication using
`p in {2,3,5,7}` proves all inequalities in `(9)`. Convexity makes this zero
the unique global minimizer. This also explains why the rational `xi`
witnesses in the external certificate have a uniform closed form rather than
being search artifacts.

## 3. Exact global maximizer in the analytic-radius variable

Write

```text
G_(p,s)(Y)=s Lambda_p(Y)-C_p(Y).                          (13)
```

For `0<x<1`,

```text
d/dx [acos(x)-x sqrt(1-x^2)]=-2 sqrt(1-x^2).              (14)
```

Therefore

```text
G_(p,s)'(Y)
 =-2 pi s+8 pi sum_(p|c, cY<1) phi(c)/c sqrt(1-c^2Y^2).  (15)
```

Each summand decreases strictly until it reaches zero at its layer boundary.
The derivative is continuous across those boundaries, tends to a positive
infinite limit as `Y` tends to zero, and tends to `-2 pi s` as `Y` tends to
`1/p`. Hence `G_(p,s)` is strictly concave and has a unique global maximizer,
the unique zero of `(15)`.

Equations `(8)` and `(15)` reduce the two-dimensional margin search to one
closed rational point and one rigorously bracketable scalar root.

## 4. The 22 positive formula checks

The derivative companion reconstructs the formulas with exact `Fraction`
arithmetic and outward-rounded `10^80` integer intervals. It verifies all 22
fixed rational witnesses printed by the source. The smallest lower endpoint is

```text
(p,s)=(5,5):
M_(5,5)>0.131799356827016832557664457131.                 (16)
```

This agrees with the source certificate. A separate 100-digit direct
implementation independently checked the four optima and all derivative
brackets. Neither path is an independent proof
that any `p`-adic zeta value is irrational; it is a conditional numerical
check of `(6)`.

## 5. Four global next-case obstructions

For each immediate next case, the exact certificate brackets the unique root
of `(15)` between rational endpoints. Concavity makes every tangent an upper
bound on `G`. A tangent at either bracketing endpoint, evaluated across the
bracket, therefore gives a rigorous upper bound at the unknown maximizer.
After inserting the exact cutoff minimizer `(8)`, the results are

| next case `(p,s)` | certified root bracket for `Y` | global upper bound for `max M_(p,s)` |
|---|---|---:|
| `(2,31)` | `[102749/5000000,205499/10000000]` | `-1.943953720720382239456347491909` |
| `(3,13)` | `[182689/5000000,365379/10000000]` | `-1.655957196700841065227181181457` |
| `(5,7)` | `[107943/2500000,431773/10000000]` | `-3.841753001087929215971397147627` |
| `(7,5)` | `[466407/10000000,58301/1250000]` | `-3.609764970272283548875193762173` |

A second, no-import exact `Fraction` implementation uses different
transcendental enclosures and narrower derivative brackets. It independently
tightens the same four global upper bounds to

```text
(2,31) < -1.943953720741976711,
(3,13) < -1.655957196706988949,
(5, 7) < -3.841753001089363045,
(7, 5) < -3.609764970278347646.                          (16a)
```

The two displays are compatible upper enclosures, not competing numerical
claims; `(16a)` is the sharper independent enclosure.

Thus no choice of the two displayed parameters can make the current margin
positive in any of these four cases. If every other term is held fixed, the
hybrid cost `tau` must decrease by at least, respectively,

```text
0.0607485537725,  0.1182826569072,
0.4802191251360,  0.6016274950454.                        (17)
```

The obstruction is one-step only: it does not classify all larger odd `s`.
It says exactly that extending each row of the claimed list by its next odd
weight requires a genuinely stronger proof ingredient, not finer rational
parameter search.

## 6. The general first-chamber wall is p=11

On the active unsaturated subchamber

```text
1<xi<(s+1)/s,                    c_p xi<1,               (18)
```

formula `(12)` is affine with positive slope, and its stationary point is
still `(8)`. Its distance from the domain boundary factors exactly:

```text
boxed: xi*-1=(s-1)(11-p)/(12s^2+(s-1)(p+1)).             (19)
```

Moreover `xi*<(s+1)/s`, and `c_p xi*<1` exactly when `p<11`. Hence this
stationary point lies inside `(18)` exactly for `p<11`; it is `xi=1` at
`p=11` and lies below the admissible interval for `p>11`.

The source's prime levels `2,3,5,7` are exactly its four selected levels, but
this is not a genus-zero classification: `X_0(13)` is also genus zero.
Equation `(19)` is an optimizer boundary and nothing more.

## 7. A p=13,s=3 all-Y obstruction

The `p=13` formal extrapolation lies beyond the first-chamber regime, so it
requires a separate global chamber calculation.

### 7.1 Optimistic one-power arithmetic cost

First replace `J_13(xi)` by its absolute upper bound `xi`, making the Hasse
term unrealistically favorable. Then

```text
S_id(xi)=7xi,
tau_id(xi)=(7xi+3I_3(xi))/8.                             (20)
```

The floor and positive-part breakpoints split `1<xi<4` at

```text
4/3, 3/2, 2, 3.                                         (21)
```

Across the five successive chambers, the derivative of the numerator in
`(20)` is

```text
9xi-11,       1,       6xi-8,       4,       3xi-5.     (22)
```

Only the first expression vanishes in its chamber. Endpoint comparison proves
the unique global minimum

```text
boxed:
xi=11/9,       I_3=305/108,       S_id=77/9,
min tau_id=613/288.                                      (23)
```

Since `(1-c_p x)_+<=1`, every actual one-power Hasse integral satisfies
`J_p(xi)<=xi`; therefore `(23)` is a lower bound for every cost of this form,
not merely for the manuscript's particular `c_13`.

### 7.2 Actual formal arithmetic cost

For `p=13`, `c_13=7/6`. For every admissible `xi>1`, the positive part in
`J_13` has ended at `6/7`, so

```text
J_13(xi)=integral_0^(6/7)(1-7x/6)dx=3/7.                (24)
```

The five chamber derivatives of the numerator `S+3I_3` are

```text
9xi-9,        3,       6xi-6,       6,       3xi-3,     (25)
```

all strictly positive on `1<xi<4`. Hence the infimum is approached but not
attained as `xi` decreases to one:

```text
boxed: inf_(1<xi<4) tau_(13,3)(xi)=515/224.             (26)
```

### 7.3 Global negativity for every continuation radius

For `s=3`, equations `(4)--(5)` give, for every `0<Y<1/13`,

```text
M=3log(13)-6pi Y-4tau-C_13(Y).                           (27)
```

The positive exponential series supplies the exact elementary bound

```text
exp(8/3)>sum_(j=0)^5 (8/3)^j/j!
        =49621/3645>13,
log(13)<8/3.                                            (28)
```

Dropping the two nonpositive terms in `(27)` and using `(23)` proves the
robust idealized obstruction

```text
boxed: M < 8-613/72=-37/72<0.                           (29)
```

Using the actual cost `(26)` sharpens the rational upper bound to

```text
boxed: M < 8-515/56=-67/56<0.                           (30)
```

Thus no `Y` rescues this formal `p=13,s=3` margin, even after idealizing away
the one-power Hasse codimension.

### 7.4 Exact analytic optimizer

In the one-collision-layer sector `1/26<=Y<1/13`, put `x=13Y`. Then

```text
C_13(Y)=(48pi/169)(acos(x)-x sqrt(1-x^2)).               (31)
```

The derivative of `-6pi Y-C_13(Y)` vanishes exactly at

```text
sqrt(1-x^2)=13/16,
x=sqrt(87)/16,
Y*=sqrt(87)/208.                                        (32)
```

Below `1/26`, the `c=13` layer alone makes the derivative positive; after
`Y*`, the one-layer derivative is negative. Hence `(32)` is the global
analytic optimizer. The optimistic exact maximum is

```text
3log(13)-613/72-(48pi/169)acos(sqrt(87)/16)
                    -(3pi sqrt(87))/208,                (33)
```

numerically `-2.0879472412217...`. Replacing `613/72` by the actual infimum
cost `515/56` gives the actual-cost supremum, numerically
`-2.770486923761...`. The rational bounds `(29)--(30)`, not these decimals,
certify negativity.

## 8. What must improve, and what does not transfer

The source audit identifies the load-bearing internal gates:

1. fixed-bundle comparison from scalar CDT sources to the vector-bundle
   source, including integral pivot control;
2. uniform neat-level descent, spreading, support bounds, and the claimed
   single unipotent action;
3. genuine-source Hasse-kernel divisibility at small auxiliary primes;
4. torsor-valued pole-grade/Cartier strictness at large primes;
5. the modular Jensen and adelic-template assembly.

The cheapest decisive positive control is an explicit `(p,s)=(5,5)` worked
case naming the frame, torsor chart, divisor, exceptional primes, normalizer,
source bases, one small-prime Hasse matrix, and one large-prime pole-grade
pivot. Re-running `(16)` does not test these gates.

The lawful repo connections are narrower:

- the exponents `v_l(lcm(1,...,N))` are the same valuation clock as in
  THM-4056, but they retain no LRC phase, owner, carrier, or safe time;
- binomial residues modulo prime powers can support a height-aware Sun
  branch-and-bound, but THM-4027 forbids any fixed congruence obstruction;
- Apéry-style denominator clearing shares the `LCM` valuation budget, while
  the Bost pivots do not supply an explicit approximant recurrence,
  nonvanishing linear form, or characteristic-root decay;
- primes versus pivots form a weighted incidence structure with ties, not an
  intrinsic tournament.

THM-4091 gives the exact elementary coordinate-change boundary behind the
LCM warning: cumulative clearing is functorial, while literal coefficient
depth above one is not.

## 9. Replay, provenance, and scope

```bash
python3 -B 04-computation/hybrid_padic_zeta_margin_frontier_thm4089.py
python3 -B -O 04-computation/hybrid_padic_zeta_margin_frontier_thm4089.py
python3 -B 04-computation/padic_zeta_next_cell_margin_obstruction_thm4089.py
python3 -B -O 04-computation/padic_zeta_next_cell_margin_obstruction_thm4089.py
python3 -B 04-computation/hybrid_padic_margin_boundary_thm4089.py
python3 -B -O 04-computation/hybrid_padic_margin_boundary_thm4089.py
```

All three normal/optimized pairs produce their frozen outputs. The primary interval primitives are adapted,
with the complete MIT copyright and permission notice in the script, from the source repository's certificate;
the derivative enclosure, concavity/tangent proof, four next cases, and output
are new to that companion. The standalone next-cell companion imports no
source code and serves as an independent exact path; the boundary companion
proves the first-chamber and `p=13` statements by exact rational chamber
algebra.

This theorem proves only the optimization and sign statements for the
explicit real function `(6)`. The external claim that 22 Kubota--Leopoldt
values are irrational remains **AUTHOR-CLAIMED / UNREFEREED**.
No LRC(14), planar Jacobian, Sun-classification, transcendence, or tournament
invariant result follows.
