---
id: THM-2236
title: "Pointwise-nested binomial minorants and the cubic vertex fan"
status: >
  PROVED + FINITE-EXACT + HOSTILE-AUDITED. For integers 1<=r<=j, the
  canonical degree-r minorant of 1_(K=0) on K=0,...,j is
  P_(r,j)(K)=(1-K)binom(j-K,r-1)/binom(j,r-1). These minorants increase
  pointwise from 1-K to the exact indicator. For 2<=r<=j, P_(r,j) is the
  unique pointwise greatest degree-r minorant which dominates P_(r-1,j).
  Without that nesting condition there is no greatest cubic for j>=4. The complete
  coefficient domain with coefficients of 1 and K fixed at 1 and -1 is an
  explicit rational polyhedron, so its normalized LP optimal values are
  monotone in degree and reach equality at r=j.
  THM-2210's unrestricted adaptive LP may be stronger. On THM-2179's
  six-peel hostile row, the nested family rises from the negative linear
  estimate through THM-2209's positive quadratic estimate to the exact
  safe mass.
source: codex-2026-07-24-pointwise-nested-binomial-minorants
depends_on: []
related:
  - THM-2210-nested-binomial-minorant-and-adaptive-moment-lp-hierarchy
  - THM-2209-sharp-quadratic-reversed-peel-and-joint-fourier-ledger
  - THM-2179-reversed-peel-relative-jackson-relation-packet
  - THM-735-bonferroni-simultaneous-multi-peel-defeats-the-clustered-non-isolated-wall
  - THM-2216-residual-capacity-hinge-gram-law
---

# THM-2236 -- pointwise-nested binomial minorants

Let `j>=1` and let `K` take values in `{0,...,j}`. For `1<=r<=j`,
consider factorial-basis polynomials

```text
P_c(K)=1-K+sum_(s=2)^r c_s binom(K,s).                (1)
```

They retain intersection moments through order `r` after integration.
The normalization `a_0=1,a_1=-1` is fixed throughout. THM-2210 instead
optimizes all dual coefficients and clips the low-degree answer by the
zero polynomial; the present pointwise-nested family is a structured
subfamily, not a replacement for that adaptive LP.

## 1. The complete coefficient polyhedron

The inequality

```text
P_c(K)<=1_(K=0)             for K=0,...,j             (2)
```

is automatic at `K=0,1` and is equivalent to the finite rational system

```text
sum_(s=2)^min(r,k) c_s binom(k,s)<=k-1
                                      for k=2,...,j.  (3)
```

Write `C_(r,j)` for this coefficient polyhedron. Equality at an integer
`k>=2` occurs exactly when its corresponding inequality in (3) is tight.
This is the whole admissible domain; no sign assumption on the `c_s` is
needed.

If `A` is a measurable set and `K(t)` counts membership in `j` measurable
peel sets, put

```text
I_s=integral_A binom(K(t),s)dt.                       (4)
```

Integration of (2) gives, for every `c in C_(r,j)`,

```text
measure({t in A:K(t)=0})
 >=I_0-I_1+sum_(s=2)^r c_s I_s.                      (5)
```

Maximizing the right side over `C_(r,j)` is therefore an LP over an exact
rational polyhedron; it is a rational LP whenever the moment data are
rational, as in the interval applications here. The feasible set for
degree `r` embeds in degree `r+1` by setting `c_(r+1)=0`, so its optimum is
monotone in `r`. At `r=j`, the exact indicator polynomial is feasible;
hence the hierarchy reaches equality.

## 2. The canonical nested hierarchy

Define

```text
P_(1,j)(K)=1-K,

P_(r,j)(K)
 =(1-K) product_(s=j-r+2)^j (1-K/s)
 =(1-K) binom(j-K,r-1)/binom(j,r-1).                 (6)
```

The binomial in (6) is its polynomial extension. On the integer domain,

```text
P_(r,j)(0)=1,
P_(r,j)(1)=0,
P_(r,j)(K)<0       for 2<=K<=j-r+1,
P_(r,j)(K)=0       for j-r+2<=K<=j.                  (7)
```

Thus every `P_(r,j)` satisfies (2), and `P_(j,j)=1_(K=0)` exactly.
Moreover,

```text
P_(r,j)(K)>=P_(r-1,j)(K)       for K in {0,...,j}.   (8)
```

Indeed, with `s=j-r+2`,

```text
P_(r,j)=P_(r-1,j)(1-K/s).
```

On `2<=K<s`, the old polynomial is negative and the new factor lies in
`(0,1)`; at `K=s` the new value is zero; for `s<K<=j` the inherited high
roots make both sides zero.

For `2<=r<=j`, the nesting condition makes (6) canonical:

> **Nested optimality.** Among degree-at-most-`r` polynomials of the form
> (1) which satisfy (2) and dominate `P_(r-1,j)` on every integer
> `0,...,j`, `P_(r,j)` is the unique pointwise greatest one.

At the `r-1` roots

```text
K=1,j-r+3,...,j
```

the old polynomial is zero. Domination and (2) force any competitor `Q`
to vanish there too. Since `Q(0)=P_(r-1,j)(0)=1`, polynomial division gives

```text
Q(K)=P_(r-1,j)(K)(1-cK).                              (9)
```

On the remaining integers `2<=K<=s=j-r+2`, domination of the negative old
polynomial forces `c>=0`, while the minorant inequality forces
`c<=1/K`. Thus `0<=c<=1/s`. The expression in (9) increases pointwise with
`c` on precisely this range, so the unique maximum is `c=1/s`, which is
`P_(r,j)`.

The factorial-basis coefficients are explicit. Put `m=r-1`,

```text
a_s=(-1)^s binom(j-s,m-s)/binom(j,m)       (0<=s<=m),
a_s=0                                      otherwise.
```

Then the coefficient of `binom(K,s)` in (6) is

```text
c_s=(1-s)a_s-s a_(s-1),                              (10)
```

with `a_(-1)=0`. In particular `c_0=1`, `c_1=-1`, and

```text
c_r=(-1)^r r/binom(j,r-1).                           (11)
```

## 3. Cubics: canonical versus data-adaptive

Write a cubic as

```text
1-K+a binom(K,2)+b binom(K,3).                       (12)
```

Its complete feasible domain is

```text
a<=1,
a+(k-2)b/3<=2/k              for k=3,...,j.          (13)
```

For `j>=3`, the finite Pareto vertices relevant to nonnegative factorial
moments are indexed by `k=2,...,j-1`:

```text
a_k=2(2k-1)/(k(k+1)),
b_k=-6/(k(k+1)).                                    (14)
```

The corresponding polynomial is

```text
-(K-1)(K-k)(K-k-1)/(k(k+1)),                         (15)
```

and is tight at `K=0,1,k,k+1`. The vertices arise by intersecting the
successive active constraints for `k` and `k+1`; direct subtraction of the
lines in (13) shows that these are all breakpoints of the upper boundary.

Neither unbounded edge can improve an actual moment objective. On the
initial edge `a=1,b<=-1`, decreasing `b` cannot improve the objective
because `I_3>=0`; if `I_3=0`, its value is already attained at the vertex
`k=2`. Pointwise,

```text
binom(K,3)<=((j-2)/3)binom(K,2),
```

so

```text
I_3<=((j-2)/3)I_2.                                  (16)
```

Along the final edge the recession coefficient is
`I_3-(j-2)I_2/3`, which is nonpositive by (16); in the equality case its
value is already attained at the last vertex. Hence the data-optimal cubic
bound is

```text
max_(2<=k<=j-1) [I_0-I_1+a_k I_2+b_k I_3].           (17)
```

The canonical nested cubic is the final vertex `k=j-1`:

```text
P_(3,j)
 =1-K+[2(2j-3)/(j(j-1))]binom(K,2)
       -[6/(j(j-1))]binom(K,3).                      (18)
```

It is the unique greatest cubic which dominates THM-2209's quadratic.
There is no pointwise greatest unrestricted cubic when `j>=4`. The vertex
`k=2` is zero at `K=2` and negative at `K=j`, while (18) is negative at
`K=2` and zero at `K=j`. A polynomial dominating both would be forced to
vanish at `K=1,2,j`; its normalization at zero makes it

```text
(1-K)(1-K/2)(1-K/j),
```

which is positive for `3<=K<=j-1`, contradicting (2). Thus "optimal cubic"
is incomplete language unless it names either the nested order or the
actual moment objective.

For `j=2`, the quadratic is already exact. For `j=3`, (18) is ordinary
exact inclusion--exclusion.

## 4. Exact hostile-row calibration

On THM-2179's six-peel row, exact interval integration gives the nested
bounds for `r=1,...,6`:

```text
-727156708364069/33883663627876920,
 733141701884261/8470915906969230,
 4875380795456321/42354579534846150,
 533835654816241/4132154100960600,
 418801891650499/3080333057079720,
 470973614624713/3388366362787692.                   (19)
```

They increase from the failed signed linear estimate to the true safe mass
in the last entry. The adaptive cubic chooses `k=4`, rather than the nested
vertex `k=5`, and improves the cubic value to

```text
40182011991889519/338836636278769200.                 (20)
```

An independent hostile audit checked (3), every endpoint and equality set
in (7), the uniqueness argument (9), the cubic recession boundary (16),
the exceptional cases `j=2,3`, all fractions in (19)--(20), and the
distinction between nested and data-adaptive optimality. QED.
