---
id: THM-4088
title: "Order-tournament arithmetic-type blindness and LRC margin density"
status: >
  PROVED ELEMENTARY OBSERVER-QUOTIENT, LIOUVILLE-CONSTANT TRANSCENDENCE, AND
  P-ADIC SHELL LEMMA + VERIFIED-EXACT + INDEPENDENTLY OBSERVER/TYPING-AUDITED.
  Ordering pairwise-distinct rational approximants,
  or ordering pairwise-distinct p-adic approximation qualities, always gives
  a transitive tournament and cannot distinguish a rational,
  algebraic-irrational, or transcendental target. Explicit monotone rational
  sequences with limits 1/2, sqrt(2), and
  the Liouville constant have the identical labelled tournament on N given
  by i->j iff i<j;
  every finite increasing rational prefix admits continuations of all three
  types. Arithmetic conclusions must retain quantitative sidecars such as
  determinants, heights, clearing, residuals, decay, and degree fences. A
  strict LRC witness interval contains times of all three types, so arithmetic
  type cannot separate the interior branch of LRC(14). No p-adic-zeta
  irrationality, transcendence beyond the displayed Liouville constant,
  LRC(14), or tournament-invariant inequality follows.
source: codex-padic-zeta-tournament-20260825
depends_on:
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
  - THM-613-margin-measure-slope-bridge
related:
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
  - HYP-3114-lrc14-irrational-transcendental-approximation-sidecar
script: 04-computation/order_tournament_arithmetic_type_thm4088.py
output: 05-knowledge/results/order_tournament_arithmetic_type_thm4088.out
script_sha256: f4e85a03b3e9e608106883fccc7d86e5cbc6ed174fa110329b0e0ca41195b764
output_sha256: 9544d89a7d290012ec71a71ba6535dd454465aed60b97209ea1f5cbc43123731
hash_basis: raw bytes
---

# THM-4088 -- a bare order tournament cannot see arithmetic type

**PROVED elementary observer quotient, Liouville-constant transcendence, and
p-adic shell lemma + VERIFIED-EXACT + INDEPENDENTLY OBSERVER/TYPING-AUDITED.**

This theorem answers the proposed rational/irrational/transcendental
tournament correspondence in the first natural model and identifies the
missing coordinates. Ordinary order is already a potential, so its
tournament has no curl. Arithmetic type lives in asymptotic magnitude and
height data that the orientation deletes.

## 1. Inheritance and the intrinsic observable

[THM-4057](THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge.md)
proves that reduced rationals are ordered coprime arcs, not tournament edges,
and that completing by ordinary order gives the transitive tournament. Its
closest hostile is the Apéry firewall: a signed Farey area is not a
denominator-cleared small nonzero linear form. The least-used relevant
sidecar is [HYP-3114](../../05-knowledge/hypotheses/HYP-3114-lrc14-irrational-transcendental-approximation-sidecar.md),
whose interval-margin observation is upgraded below from rational grid points
to all three arithmetic types.

Let

```text
r_i=a_i/b_i in Q,        b_i>0,        gcd(a_i,b_i)=1.       (1)
```

The intrinsic pairwise order observable is

```text
i -> j  iff  r_i<r_j.                                      (2)
```

### Theorem 1.1 (order collapse)

For every finite family of distinct rationals, `(2)` is a transitive
tournament. If `r_1<r_2<...`, its labelled orientation is exactly

```text
i -> j  iff  i<j.                                          (3)
```

This follows immediately from transitivity of `<`. In homogeneous
coordinates the arc sign is

```text
Delta_ij=a_j b_i-a_i b_j,
r_j-r_i=Delta_ij/(b_i b_j).                               (4)
```

The tournament keeps only `sgn(Delta_ij)`. It destroys the determinant
magnitude and both heights, which together reconstruct the rational gap in
`(4)`.

## 2. One labelled tournament, three different limit types

There are explicit strictly increasing rational sequences with each of the
three possible arithmetic limit types.

### Rational limit

For `n>=1`, put

```text
u_n=1/2-1/(n+2).                                           (5)
```

Then `u_n` increases strictly to the rational number `1/2`.

### Quadratic-irrational limit

Put `a_0=b_0=1` and

```text
a_(k+1)=3a_k+4b_k,       b_(k+1)=2a_k+3b_k.               (6)
```

Multiplication by `3+2sqrt(2)` gives

```text
a_k^2-2b_k^2=-1,
a_(k+1)b_k-a_kb_(k+1)=2.                                 (7)
```

Hence `a_k/b_k` is strictly increasing. Also

```text
2-(a_k/b_k)^2=1/b_k^2,                                   (8)
```

so it tends to the algebraic irrational `sqrt(2)`.

### Transcendental limit

Let

```text
ell_n=sum_(k=1)^n 10^(-k!),
L=sum_(k=1)^infinity 10^(-k!).                            (9)
```

The rational truncations increase strictly to `L`. With
`q_n=10^(n!)`,

```text
0<L-ell_n<2*10^(-(n+1)!)=2/q_n^(n+1).                   (10)
```

If `L=A/B` were rational, every distinct `p/q` would satisfy
`|L-p/q|>=1/(Bq)`, contradicting `(10)` for all sufficiently large `n`. If
`L` were algebraic irrational of degree `d`, let `P in Z[X]` be its primitive
minimal polynomial. For every rational `a/q!=L`, the nonzero integer
`q^d P(a/q)` has magnitude at least one. On a fixed compact neighborhood of
`L`, the mean-value theorem gives `|P(a/q)|<=C|a/q-L|`; hence
`|L-a/q|>=C^(-1)q^(-d)`. For all sufficiently large `n` with `n+1>d`, the
upper bound in `(10)` is smaller than this lower bound, again a contradiction.
Thus `L` is transcendental.

All three sequences satisfy `(3)`. Therefore, after the rational decorations
are forgotten, even their complete countably infinite labelled index
tournaments are identical. No invariant of that bare tournament can
distinguish their limit types. The full decorated approximation sequence is a
different object and can of course determine its limit.

## 3. Every finite prefix has all three continuations

Let

```text
r_1<...<r_N<B                                              (11)
```

be rationals. The interval `(r_N,B)` contains a rational, an algebraic
irrational (a small rational translate and scale of `sqrt(2)`), and a
transcendental (the same construction with `L`). For any selected target
`alpha` in that interval, set `s_0=r_N`, and recursively choose a rational
`s_k` with

```text
max(s_(k-1),alpha-2^(-k)) < s_k < alpha.                  (12)
```

Density of `Q` supplies the choice. Then `(s_k)` is strictly increasing and
converges to `alpha`. It appends `(11)` without changing the labelled
transitive order pattern. Consequently no finite ordinary-order or depth-free
tournament prefix can certify an arithmetic type.

This does not say that every continued-fraction statistic is useless. An
infinite tail plus a degree/height fence or irrationality measure is precisely
the missing sidecar.

## 4. What an irrationality proof must retain

For real `xi`, integer pairs `(p_n,q_n)` with `q_n!=0` give an elementary
irrationality criterion:

```text
0<|q_n xi-p_n| -> 0  implies  xi notin Q.                 (13)
```

Indeed, if `xi=A/B` in lowest terms, every nonzero residual in `(13)` has
absolute value at least `1/B`. An Apéry construction may begin with rational
coefficients rather than integer pairs, in which case its LCM or valuation
clock is needed to reach `(13)`. Transcendence needs still more: a degree- and
height-dependent lower bound compared against an upper bound such as `(10)`.
Signs of `Delta_ij` provide none of these magnitudes.

Thus the faithful hierarchy is

```text
orientation sign
  < sign + |Delta_ij| + (b_i,b_j)
  < target residual + denominator clearing
  < residual rate + algebraic degree/height fence.        (14)
```

## 5. The p-adic quality tournament collapses too

There is no field-compatible order on `Q_p`, so `(2)` is not a native
observable for p-adic values. The nearest lawful substitute ranks rational
approximants by error valuation.

Fix a prime `p`, normalize `v_p(p)=1`, fix `xi in Q_p`, and take any strictly
increasing integer sequence `(h_n)`. Density of `Q` in `Q_p` gives `q_n in Q`
with `v_p(xi-q_n)>h_n`. Set

```text
r_n=q_n+p^(h_n).                                         (15)
```

The ultrametric inequality gives the exact equality

```text
v_p(xi-r_n)=h_n.                                         (16)
```

Therefore orienting `i->j` when the `j`th approximation has larger valuation
again gives `(3)`, for every p-adic target. Rational, algebraic-irrational,
and transcendental p-adic targets can all carry the same quality tournament.
Here arithmetic type is measured over `Q` inside `Q_p`.
All three types occur. Rational targets come from `Q`. For odd `p`, Hensel's
lemma applied at zero to `X^2-X-p` gives a root in `Q_p`, and its discriminant
`1+4p` is not a rational square; for `p=2`, use `X^2-X-4` with discriminant
`17`. Finally, the elements of `Q_p` algebraic over `Q` form a countable set,
whereas `Q_p` is uncountable, so transcendental targets exist.
Small p-adic error alone cannot contradict rationality; an adelic height or
product-formula sidecar is load-bearing.

## 6. LRC witness times of all three types

For a finite nonempty set of nonzero integer speeds `S`, define

```text
F_S(t)=min_(s in S) ||s t||,       M=max_(s in S)|s|.     (17)
```

Distance to the nearest integer is one-Lipschitz. Hence every term in the
minimum is `|s|`-Lipschitz and

```text
|F_S(t)-F_S(u)|<=M|t-u|.                                 (18)
```

If, for some target threshold `c`,

```text
F_S(t_0)>=c+delta,       delta>0,                         (19)
```

then every `u` with `|u-t_0|<delta/M` satisfies `F_S(u)>c`.
Every nonempty real interval contains rational, algebraic-irrational, and
transcendental numbers, so the strict witness component contains times of all
three types.

For LRC(14), this proves that the arithmetic type of an **interior** lonely
time cannot distinguish a speed family or close the conjecture. The survivor
is the boundary branch `F_S(t)=1/14`, where the positive margin vanishes, and
the quantitative task of landing on a prescribed finite denominator shell.
The interval-local Stern sign sum suggested by THM-4057/4071 remains a lawful
open problem because it retains denominator and interval labels; it is not a
classification by arithmetic type.

## 7. Transfer and loss ledger

| source | target/map | preserved | destroyed / required sidecar |
|---|---|---|---|
| rational approximants | order tournament `(2)` | pairwise order | determinant magnitude, height, target, residual rate |
| monotone sequence | index-order tournament on `N` | index order type | limit value and arithmetic type |
| p-adic approximants | distinct valuation-quality ranking | relative error depth | error units, archimedean height, product-formula budget |
| strict LRC witness | open interval via `(18)` | loneliness above threshold | boundary owner/equality and denominator shell |
| linear form | its sign | one orientation bit | nonzero magnitude, clearing, decay, degree/height fence |

If the claimed p-adic-zeta certificate margins are pairwise distinct, sorting
them creates a transitive priority tournament; with ties it creates a
transitive preorder until an explicit tie-break is supplied. Either object
ranks proof slack, not values. The matching-logic reachability graph in
arXiv:2608.13306v1 is not a tournament at all: it may have loops, missing
pairs, and both directions.

## 8. Replay and scope

```bash
python3 -B 04-computation/order_tournament_arithmetic_type_thm4088.py
python3 -B -O 04-computation/order_tournament_arithmetic_type_thm4088.py
```

The companion checks identical order matrices for exact prefixes of `(5)`,
`(6)`, and `(9)`, determinant reconstruction, `28,812` exact Lipschitz pairs,
and `45` rational p-adic controls. It explicitly does not infer transcendence
by computation.

This theorem proves a no-go and identifies the positive sidecars. It proves
no irrationality or transcendence result for a new named constant, no claim
from the external p-adic-zeta manuscript, no LRC(14), no planar Jacobian
case, and no inequality for `H`, `disc`, or rooted Pfaffian energy.
