---
id: THM-2169
title: "A bounded relation on every deletion of a zero-safe LRC(14) row"
status: >
  PROVED from THM-2144 and THM-2164. For every coordinate i of a distinct
  positive thirteen-speed row with zero-measure safe set, the other twelve
  speeds carry a nonzero integer relation of height at most 11130. More
  sharply, the height is at most 4558 unless v_i is an integer multiple
  q v_j with 2<=q<=105. Thus every deletion core relevant to a hypothetical
  LRC(14) counterexample lies on an explicit bounded relation hyperplane.
  Under radix descent the deletion relation has l1 norm at most 133560 and
  fewer than 267120 carry states; the quotient-owner masks remain necessary.
  This is a structural reduction, not a finite speed bound, Q108, or LRC(14).
source: opus-2026-07-24-puzzle-atlas
depends_on:
  - THM-2144-anisotropic-selberg-kraft-relation-box
  - THM-2164-relative-packet-rank-harvesting
related:
  - THM-2162
  - THM-2163
  - THM-2166
  - THM-2167
---

# THM-2169 -- bounded relations on every deletion

Put

```text
J=[1/14,13/14],
S(v)={t in T:v_i t in J for every i},

Lambda(v)={m in Z^13:m.v=0}.                         (1)
```

Let `v=(v_1,...,v_13)` have pairwise distinct positive integer coordinates,
and suppose

```text
mu(S(v))=0.                                           (2)
```

## 1. Coordinatewise dichotomy

Fix a coordinate `i`. THM-2144, applied with cap `1` at coordinate `i` and
cap `105` elsewhere, supplies a nonzero relation

```text
r in Lambda(v),
|r_i|<=1,              ||r||_infinity<=105.           (3)
```

Divide by the gcd of its coefficients and call the resulting primitive
relation `w`. The bounds in (3) persist.

If `w_i=0`, then `w` itself is already a relation on the deletion
`v\{v_i}`, of height at most `105`.

Suppose `w_i!=0`; then `w_i=+1` or `-1`.

### Support at least three

If `|supp(w)|>=3`, THM-2164's whole-subtorus packet lemma gives

```text
s in Lambda(v)\Qw,              ||s||_infinity<=43.   (4)
```

Eliminate coordinate `i`:

```text
u=s_i w-w_i s.                                        (5)
```

Then

```text
u in Lambda(v),        u_i=0,
```

and `u!=0`, because otherwise (5) would make `s` rationally proportional to
`w`. Coordinatewise,

```text
|u_k|
 <=|s_i||w_k|+|w_i||s_k|
 <=43*105+43
 =4558.                                               (6)
```

Thus the deletion has a nonzero height-`4558` relation.

### Support two

If `|supp(w)|=2`, write its support as `{i,j}`. Positivity of the speeds
forces the two coefficients to have opposite signs. Since `|w_i|=1`,

```text
v_i=q v_j,                 q=|w_j|.                   (7)
```

Distinctness excludes `q=1`, while (3) gives

```text
2<=q<=105.                                            (8)
```

We have proved the sharp coordinatewise fork:

```text
either
  Lambda(v) contains nonzero u with u_i=0,
  ||u||_infinity<=4558,

or
  v_i=q v_j for some j!=i and 2<=q<=105.              (9)
```

The lock in the second branch is directed: it says the selected speed
`v_i`, not merely one of the pair, is the bounded multiple.

## 2. Uniform deletion relation

The divisibility-lock branch still has a bounded deletion relation. By
THM-2164,

```text
dim_Q span{m in Lambda(v):||m||_infinity<=105}>=2.    (10)
```

Choose `s` of height at most `105` independent of the support-two `w`, and
again use (5). It is nonzero, avoids coordinate `i`, and satisfies

```text
||u||_infinity
 <=105*105+105
 =11130.                                              (11)
```

Combining the three cases proves:

> **Deletion theorem.** For every `i=1,...,13`, there is a nonzero
> `u^(i) in Lambda(v)` such that
>
> ```text
> u^(i)_i=0,             ||u^(i)||_infinity<=11130.   (12)
> ```

Equivalently, every twelve-speed deletion of `v` lies on an integer
relation hyperplane of coefficient height at most `11130`.

## 3. Exact interfaces

### Q108 / deletion induction

Any zero-safe row relevant to LRC(14) therefore has **all thirteen** of its
twelve-speed deletion cores in the bounded-relation class (12). A
deletion-based proof of LRC(14) does not need a Q108-type stability theorem
for arbitrary twelve-sets first; it is enough to prove the required
measure/stability statement on twelve-sets carrying a height-`11130`
relation, together with the directed locks (7).

This is a genuine restriction but not a finite search. A bounded relation
still leaves eleven rational degrees of freedom and unbounded speeds.

### Radix carry

Regard `u^(i)` as a twelve-coordinate relation after deleting its zero
coordinate. Then

```text
||u^(i)||_1<=12*11130=133560.                         (13)
```

THM-2163 converts it, in every base `q>=2`, into an exact terminating carry
path with fewer than

```text
2*133560=267120                                       (14)
```

possible integer carry states. In the generic branch (6), the corresponding
figures improve to

```text
||u||_1<=12*4558=54696,
fewer than 109392 carry states.                       (15)
```

These are state-alphabet bounds, not depth bounds. THM-2163's nested
quotient-owner masks remain indispensable for termination; (12) does not
make the speed search finite by itself.

## 4. Boundary audit

Each hypothesis is used:

- zero safe measure invokes the Selberg relations;
- positivity turns a support-two relation into a divisibility lock;
- distinctness changes `1<=q<=105` into `2<=q<=105`;
- primitivity of `w` makes `w_i` exactly `+/-1` when nonzero;
- independence of `s,w` makes the eliminated relation nonzero.

The constants are literal elimination costs:

```text
43*105+43=4558,
105*105+105=11130.
```

No determinant division, saturation denominator, or unproved sparsity
claim is hidden in them. QED.
