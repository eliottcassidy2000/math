---
id: THM-2210
title: "Nested adaptive moment LP and sharp binomial-prefix minorants"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For every finite
  multiplicity packet, the least zero-mass compatible with factorial moments
  through degree d is a finite primal/dual LP value L_d. For one fixed
  realizable full packet these values are attained, nested, and exact at full
  Boolean depth by binomial inversion. With the alternating lower-degree
  coefficients fixed, the sharp top coefficient is d/j for even d and -1
  for odd d. Those explicit prefix minorants are not themselves nested.
  Degree two is exactly THM-2209's quadratic bound clipped at zero; degree
  three is the first adaptive level that can improve it. This supplies an
  exact escalation framework, not a uniform positive estimate or a proof of
  LRC(14).
source: codex-2026-07-24-peel-moment-hierarchy
depends_on: []
related:
  - THM-2209-sharp-quadratic-reversed-peel-and-joint-fourier-ledger
  - THM-2216-residual-capacity-hinge-gram-law
script: 04-computation/lrc14_adaptive_moment_lp_hierarchy_thm2210.py
output: 05-knowledge/results/lrc14_adaptive_moment_lp_hierarchy_thm2210.out
script_sha256: 0fa36f21c161b36997b7a419c5c44ecf4d3f48b98cdcfa26eabee6c5078fffe4
output_sha256: cf618191143679b65806a4ed88dfd54aaaa291918225e4e5ec6099b6572533a4
hash_basis: working-tree bytes (LF)
---

# THM-2210 -- the adaptive factorial-moment LP

The word **nested** applies to the optimal LP values below. It does not
apply to the displayed alternating binomial-prefix polynomials. Confusing
those two statements produces a false monotonicity claim, whose minimal
positive-degree witness is recorded in Section 5.

## 1. Multiplicity packets and factorial moments

Fix `j>=0` and a nonnegative packet

```text
p=(p_0,...,p_j).
```

Its factorial moments are

```text
M_r=sum_(k=0)^j binom(k,r)p_k,          0<=r<=j.     (1)
```

No probability normalization is required. If `M_0>0`, dividing every
`p_k`, `M_r`, and bound below by `M_0` gives the normalized version. If
`M_0=0`, nonnegativity forces the zero packet and every assertion is
trivial.

The LRC interpretation is exact. Let `G` be a measurable body set, let
`D_1,...,D_j` be peel danger sets, and put

```text
K=sum_(i=1)^j 1_(D_i),
p_k=measure(G intersection {K=k}).                   (2)
```

Then

```text
p_0=measure(G intersection {K=0}),                  (3)

M_r=sum_(|A|=r)
 measure(G intersection intersection_(i in A)D_i).  (4)
```

Equation (4) follows by double-counting: a point of multiplicity `k`
belongs to exactly `binom(k,r)` of the `r`-fold intersections. Thus the
moments retain the symmetric overlap histogram, but not the labels of the
subsets producing it.

## 2. The primal and dual hierarchy

For `0<=d<=j`, define

```text
L_d(M)
 =min x_0

subject to
 x_k>=0                                      (0<=k<=j),
 sum_(k=0)^j binom(k,r)x_k=M_r              (0<=r<=d).
                                                               (5)
```

The original packet `p` is feasible. The row `r=0` fixes

```text
sum_k x_k=M_0,
```

so the feasible set is closed and bounded. The minimum is therefore
attained and

```text
0<=L_d(M)<=p_0.                                      (6)
```

The exact dual is

```text
L_d(M)
 =max sum_(r=0)^d a_r M_r

subject to
 q_a(k):=sum_(r=0)^d a_r binom(k,r)
          <=1_(k=0)                         (0<=k<=j),
 a_r unrestricted.                                  (7)
```

To see the sign, write the Lagrangian as

```text
x_0+sum_r a_r(M_r-sum_k binom(k,r)x_k).
```

Its infimum over `x_k>=0` is finite exactly when the coefficient of every
`x_k` is nonnegative, which is the minorant inequality in (7). The zero
polynomial is dual-feasible. Finite-dimensional LP strong duality gives
equality and dual attainment.

Equivalently, every feasible `q_a` satisfies

```text
sum_r a_rM_r
 =sum_k q_a(k)p_k
 <=p_0.                                              (8)
```

The LP selects the best degree-`d` binomial-basis minorant for the actual
moment vector rather than fixing its coefficients in advance.

## 3. Nesting and full-depth exactness

For one fixed realizable full moment vector,

```text
0<=L_0<=L_1<=...<=L_j=p_0.                           (9)
```

The primal feasible sets shrink when one more moment is imposed.
Equivalently, every degree-`d` dual polynomial embeds at degree `d+1` by
setting its new coefficient to zero. This is the only nesting theorem
claimed here; unrelated truncated moment vectors cannot be compared.

At full depth the binomial moment matrix is unit triangular. Binomial
inversion gives

```text
p_k
 =sum_(r=k)^j (-1)^(r-k)binom(r,k)M_r,               (10)

p_0=sum_(r=0)^j(-1)^rM_r.                            (11)
```

Thus the full packet has one feasible histogram and `L_j=p_0`. This is
finite inclusion--exclusion, not a uniform positivity theorem.

## 4. Sharp fixed-prefix minorants

Fix `1<=d<=j` and consider only the one-parameter family

```text
q_(d,c)(k)
 =sum_(r=0)^(d-1)(-1)^r binom(k,r)+c binom(k,d).
                                                               (12)
```

At `k=0` it equals one, and at `1<=k<d` it equals zero. For `k>=d`,

```text
sum_(r=0)^(d-1)(-1)^r binom(k,r)
 =(-1)^(d-1)binom(k-1,d-1)
 =(-1)^(d-1)(d/k)binom(k,d).                         (13)
```

Consequently the largest coefficient for which

```text
q_(d,c)(k)<=1_(k=0)             for 0<=k<=j          (14)
```

is

```text
c_(j,d)=
 d/j,                  d even,
 -1,                   d odd.                        (15)
```

For even `d`, the binding node is `k=j`; for odd `d`, it is `k=d`.
Explicitly, for `k>=1`,

```text
q_(d,c_(j,d))(k)
 =-(j-k)/j binom(k-1,d-1),          d even,
 =-binom(k-1,d),                    d odd.            (16)
```

The resulting valid fixed-prefix bound is

```text
B_d
 =sum_(r=0)^(d-1)(-1)^rM_r+c_(j,d)M_d
 <=L_d<=p_0.                                          (17)
```

At `d=j`, coefficient (15) is `(-1)^j`, so (17) becomes the exact
inversion (11).

There are two degree-zero edge cases. If `d=0<j`, the sharp constant
minorant is zero and `L_0=0`. If `d=j=0`, the sharp constant is one and
`L_0=M_0=p_0`. The expression `d/j` must not be used at `(0,0)`.

Sharpness in (15) concerns only the fixed-prefix family (12). It does not
say that `q_(d,c_(j,d))` is the instance-optimal dual polynomial.

## 5. The fixed prefixes are not nested

Take `j=4` and the normalized packet

```text
p_0=p_4=1/2,                    p_1=p_2=p_3=0.
```

Its moments begin

```text
(M_0,M_1,M_2,M_3)=(1,2,3,2).                        (18)
```

The fixed-prefix values decrease:

```text
B_2=1/2,                       B_3=0.                (19)
```

Pointwise, this is already visible at the terminal node:

```text
q_2(4)=0,                      q_3(4)=-1.             (20)
```

Nevertheless `B_2=p_0` forces `L_2=p_0`, and nesting (9) forces

```text
L_2=L_3=1/2.                                         (21)
```

Thus adaptive optimal values are nested while the canonical alternating
certificates are not. If degree zero is included, the smaller example
`j=2,p_2=1` has `B_0=0>B_1=-1`.

## 6. What happens at low degree

For `j>=1`,

```text
L_1=max(0,M_0-M_1).                                 (22)
```

The two displayed terms are dual values from `q=0` and `q(k)=1-k`.
For attainability, if `M_1<=M_0`, put mass `M_0-M_1` at zero and `M_1`
at one. If `M_1>=M_0`, realizability puts `M_1/M_0` in `[1,j]`; interpolate
between its two neighboring positive integers to realize the same mass
and mean with no mass at zero.

For `j>=2`, THM-2209's quadratic expression is

```text
B_2=M_0-M_1+(2/j)M_2,

L_2=max(0,B_2).                                      (23)
```

For `j=2` this is full-depth inversion. For `j>=3`, inspect the vertices
of the three-dimensional dual polyhedron. Three active positive nodes
force the zero quadratic. If the node zero and two positive roots `a<b`
are active, `q(0)=1` makes the leading coefficient positive. Feasibility
on every integer `1,...,j` then forces `a=1,b=j`, which is exactly the
quadratic in (23). These are the only vertices, proving the formula.

Thus degree-two adaptation adds only the trivial nonnegativity clip to
THM-2209. The first possible genuine gain is cubic. For example, with

```text
j=4,                  p_0=p_3=1/2,                  (24)
```

one has

```text
(M_0,M_1,M_2,M_3)=(1,3/2,3/2,1/2),

L_2=1/4,                       L_3=1/2=p_0.          (25)
```

This is an abstract Boolean packet demonstrating strict cubic gain. It is
not asserted to arise from circle danger combs.

## 7. LRC scope and the adaptive stopping rule

For a concrete reversed-peel row, compute the body-relative moments (4)
in increasing degree and solve (5), stopping as soon as `L_d>0`. The
canonical THM-2209 hostile row already stops at `d=2`; THM-2210 does not
close a new LRC row by itself.

The escalation retains the multiplicity histogram but forgets:

```text
which labelled peel subsets meet,
endpoint owner and current,
guard-hole/root-sheet alignment,
and continuation data between depths.               (26)
```

Those sidecars remain mandatory. Computing every `M_r` through `r=j` can
cost essentially the direct Boolean union calculation. Full-depth
exactness therefore proves termination of the information hierarchy, not
its efficiency or positivity.

THM-2216's hinge/Gram law addresses a different order statistic: top
terminal-label capacities over residual owners. The present LP addresses
the peel multiplicity histogram. Their shared binomial/hinge syntax does
not identify their underlying packets.

## 8. Exact referee

Run

```bash
python3 04-computation/lrc14_adaptive_moment_lp_hierarchy_thm2210.py
python3 -O 04-computation/lrc14_adaptive_moment_lp_hierarchy_thm2210.py
```

Using exact rational arithmetic, the companion independently enumerates
primal bases and dual active-node bases for all `1,086` nonzero packets

```text
p in {0,1,2}^(j+1),                  0<=j<=5.
```

At every degree it checks primal-dual equality, (9), (11), (17), and the
closed forms (22)--(23). A separate sweep checks the sharp coefficient
and its failure after adding `1/100` for all `231` pairs
`0<=d<=j<=20`. It also freezes witnesses (18)--(25). Ordinary and
optimized runs are byte-identical to the stored output.

The finite referee audits indexing and edge cases. Finite LP duality and
the binomial identities above supply the all-`j` proof. QED.
