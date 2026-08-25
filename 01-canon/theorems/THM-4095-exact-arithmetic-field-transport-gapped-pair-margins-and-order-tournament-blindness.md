---
id: THM-4095
title: "Exact arithmetic-type transport, gapped pair margins, and order-tournament blindness"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every fixed finite LRC observer
  F_S(t)=min_(v in S)||vt|| preserves the singly generated field exactly:
  Q(F_S(t)-beta)=Q(t) for rational beta. Consequently each arithmetic field
  fiber has its exact typed image and is dense below the rational optimum.
  In contrast, primitive two-speed optimized 1/3-margins form the discrete
  set {1/6} union {1/6-1/(2q): q odd, q>=3}, with gap (0,1/15).
  Order tournaments see neither distinction: even infinite transitive prefix
  towers realize every reciprocal mass in (0,1]. This is an information-loss
  theorem, not LRC(14) or a rejection of intrinsically defined tournaments.
source: codex-padic-zeta-tournament-20260825
depends_on: []
related:
  - THM-1155-threespeed-exhaustive-and-ceiling
  - THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density
  - THM-4096-twisted-padic-affine-lrc-ray-optimizer-and-next-case-obstruction
script: 04-computation/lrc_margin_arithmetic_type_tournament_blindness_thm4095.py
output: 05-knowledge/results/lrc_margin_arithmetic_type_tournament_blindness_thm4095.out
script_sha256: 15ff1d8e382de6f9183f148619e20a32e163a10f4bed71a306e950f6aa62f591
output_sha256: 35e8b49afb170ee032feb0166c2c87325579f39273a4410068100928010c4ae3
hash_basis: raw LF bytes
audit: >
  ACCEPT after explicit greedy-remainder inequalities and a one-speed
  minimality qualifier. The independent hostile audit checked field-fiber
  surjectivity/density, rational maximizers, the primitive-pair upper bound,
  the complete optimized spectrum, and normal/-O/frozen transcript identity.
---

# THM-4095 -- field transport is exact, optimized pair margins are gapped,
# and order tournaments see neither

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** Three
superficially conflicting observations become compatible once the
quantifiers are separated:

1. varying the time for one fixed finite speed set gives a full interval;
2. varying the speed set and then optimizing can give a discrete spectrum;
3. replacing the metric data by its order tournament destroys both facts.

## 1. A fixed finite observer preserves the generated field

Let `S` be a finite nonempty set of positive integers, choose the canonical
representative `t in [0,1)` of a point of `R/Z`, and put

```text
F_S(t)=min_(v in S)||vt||,
M(S)=max_(t in R/Z) F_S(t),
Delta_(S,beta)(t)=F_S(t)-beta,          beta in Q.           (1)
```

The range is exactly

```text
boxed: Delta_(S,beta)(R/Z)=[-beta,M(S)-beta].               (2)
```

Indeed, `F_S` is continuous, the circle is connected, `F_S(0)=0`, and the
maximum exists.

More rigidly, for every `t`, choose an active speed `v` and an integer `k`
such that

```text
F_S(t)=||vt||=epsilon(vt-k),        epsilon in {+1,-1}.      (3)
```

Then

```text
t=(k+epsilon F_S(t))/v.                                    (4)
```

Since the reverse inclusion `F_S(t) in Q(t)` is immediate from `(3)`, and
`beta` is rational,

```text
boxed: Q(Delta_(S,beta)(t))=Q(F_S(t))=Q(t).                 (5)
```

This is pointwise field equality, not merely preservation of the coarse
labels rational/algebraic/transcendental.

The tent functions have rational walls `k/(2v)`. On every complementary
cell they are affine over `Q`, so their finite lower envelope has rational
vertices. It cannot have a constant open segment because every active affine
piece has nonzero slope. Hence

```text
M(S) is rational, and every maximizing time is rational.                 (6)
```

The second assertion also follows from `(5)` after the first. Similarly,
the lower endpoint in `(2)` is attained only at rational times.

## 2. Exact typed images and density

For any singly generated real field `K`, write

```text
T_K={t in [0,1):Q(t)=K}.                                   (7)
```

Equations `(2)` and `(5)` give the exact image, in both directions:

```text
boxed:
Delta_(S,beta)(T_K)
 ={y in [-beta,M(S)-beta]:Q(y)=K}.                          (8)
```

The reverse direction in `(8)` is worth recording. Given such a `y`, the
intermediate-value theorem supplies a time with `F_S(t)=y+beta`; equation
`(5)` then forces that time into `T_K`.

If `K=Q(x)`, the affine generators `a+bx`, with rational `a,b` and `b!=0`,
are dense in `R` and still generate `K`. Therefore every field fiber in
`(8)` is dense in the open interval

```text
(-beta,M(S)-beta).                                         (9)
```

For `K=Q`, the typed image also contains both rational endpoints. In
particular, rational, algebraic-irrational, and transcendental times map
exactly to values of the same respective type; all three value types are
dense in the interior, while only rational times attain the optimum.

There is also an effective strict-witness version. Let `V=max(S)`. Every
tent `||vt||` is `v`-Lipschitz in circular distance, hence `F_S` is
`V`-Lipschitz. If

```text
eta=F_S(t_0)-beta>0,                                       (10)
```

then

```text
d_circle(t,t_0)<eta/V  implies  F_S(t)>beta.                (11)
```

Every nonempty component of the strict witness set therefore contains
times from every singly generated field fiber. Moreover, for every integer

```text
q>V/(2eta),                                                 (12)
```

the nearest grid point `r/q` to `t_0` lies in the ball `(11)`. Thus any
strict real witness compiles to an explicit rational grid once its margin is
known.

## 3. Optimization changes the quantifier and destroys density

Now specialize to a primitive pair `S={a,b}`, with `1<=a<b` and
`gcd(a,b)=1`. Its exact optimum is

```text
boxed:
M({a,b}) = 1/2                         if a,b are both odd,
M({a,b}) = (a+b-1)/(2(a+b))            if a,b have mixed parity. (13)
```

The odd case is attained at `t=1/2`. For the mixed case put `q=a+b`, which
is odd. To prove the upper bound, suppose both distances were greater than
`1/2-1/(2q)`. There would be integers `r,s` and errors `u,w` with

```text
at=r+1/2+u,       bt=s+1/2+w,
|u|,|w|<1/(2q).                                           (14)
```

But

```text
0=b(at)-a(bt)=br-as+(b-a)/2+bu-aw.                         (15)
```

Because `b-a` is odd, `(15)` forces `bu-aw` to be a half-integer, whereas
its absolute value is strictly less than `(a+b)/(2q)=1/2`, a contradiction.
For equality, choose `k mod q` with

```text
ak=(q-1)/2 mod q                                           (16)
```

and take `t=k/q`. Since `b=-a mod q`, the two residues have circular
distance `(q-1)/(2q)`. This proves `(13)`.

Subtract the two-runner LRC threshold `1/3`. Every odd `q>=3` occurs through
the primitive mixed pair `{1,q-1}`, while an odd pair supplies `1/6`.
Therefore the complete optimized-margin spectrum is

```text
boxed:
{M({a,b})-1/3:gcd(a,b)=1, a<b}
 ={1/6} union {1/6-1/(2q):q>=3 odd}.                        (17)
```

It contains `0` at `q=3`, next contains `1/15` at `q=5`, and has the open
gap

```text
boxed: (0,1/15).                                           (18)
```

Thus the fixed-observer density theorem `(8)--(9)` does **not** imply that
optimized margins are dense when the speed set varies. The optimization
quantifier is the destroyed coordinate.

## 4. The order tournament is completely blind to scale and type

Let `A=(a_1<a_2<...)` be a strictly increasing sequence of positive
integers. Orient its prefix by

```text
i -> j  iff  a_i<a_j.                                      (19)
```

Every prefix tournament is the same transitive tournament, for every such
sequence. Nevertheless, for every real `x in (0,1]` there is an infinite
strictly increasing sequence with

```text
sum_(n>=1) 1/a_n=x.                                        (20)
```

Apply the greedy Egyptian-fraction algorithm to `x`. If the current remainder
is `x_j>0` and `m_j=ceil(1/x_j)`, then

```text
1/m_j<=x_j<1/(m_j-1).                                     (21)
```

When `x_(j+1)=x_j-1/m_j` is positive,

```text
x_(j+1)<1/[m_j(m_j-1)],
m_(j+1)>m_j(m_j-1)>=m_j.                                  (22)
```

The `m_j=1` case is `x=1` and terminates immediately. Otherwise a
nonterminating denominator sequence tends to infinity, so
`x_j<1/(m_j-1)->0` and its unit fractions sum to `x`. Thus the denominators
are strictly increasing unless the expansion terminates. If the last term is
`1/m`, replace it by the infinite telescoping tail

```text
1/m=sum_(r=m)^infinity 1/(r(r+1)).                         (23)
```

The new denominators start at `m(m+1)>m`, so they remain strictly increasing
after the earlier greedy terms. If the algorithm never terminates, its
remainders tend to zero and it already gives `(20)`.

Hence the identical tournament tower realizes rational, algebraic
irrational, and transcendental reciprocal masses. It even coexists with a
divergent mass, by taking `a_n=n`. For a one-speed set `M({a})=1/2`
regardless of `a`, so two vertices are minimal for an optimized-margin
hostile. The smallest such hostile is

```text
{1,2}: M-1/3=0,
{1,3}: M-1/3=1/6.                                          (24)
```

Both have the same two-vertex order tournament. Parity, metric gaps,
reciprocal mass, optimum, and arithmetic type have all been discarded.

## 5. Boundary and exact audit

This does not reject tournament analysis when an intrinsic binary relation
is present. It says that the cosmetic orientation supplied by numerical
order has no recovery theorem unless metric/parity sidecars are retained.
The connection ledger is:

```text
source:       finite metric observer F_S, or an increasing speed sequence
quotient:     retain only pairwise numerical order
preserved:    cardinality and transitivity
destroyed:    parity, gaps, scale, reciprocal mass, fields, optimized margin
sidecar:      active affine branch (v,k,epsilon), or the actual speed labels
recovery:     exact from Delta for fixed S; impossible from order alone
hostile:      {1,2} versus {1,3}.                            (25)
```

Reproduce the finite exact referee with

```bash
python3 -B 04-computation/lrc_margin_arithmetic_type_tournament_blindness_thm4095.py
python3 -B -O 04-computation/lrc_margin_arithmetic_type_tournament_blindness_thm4095.py
```

It independently reconstructs every primitive pair with `a<b<=80`: `1,965`
pairs split into `651` odd and `1,314` mixed-parity cases. It also checks
`12,956` active affine recovery maps, exact rational-grid witnesses, six
greedy/telescoping controls, the first seventy-nine mixed-parity spectrum
points, and the hostile pair `(24)`. All calculations use exact rational
arithmetic and explicit `require` gates that remain active under `-O`.

**QED.**
