---
id: THM-4331
title: "LRC(14) safe-component endpoint-denominator odd-wall escape"
status: >
  PROVED ELEMENTARY COMPONENT-WALL CERTIFICATE + INDEPENDENTLY AUDITED;
  LRC(14) OPEN. A positive-length component of the 1/14-safe set with
  endpoint denominators Q_L,Q_R accepts every distinct positive odd-tail
  pair a<b after doubling whenever
  bW>=2/7-1/min(Q_L,Q_R). The stronger two-endpoint certificate is valid
  with strict > only; an exact equality containment shows why >= is false
  for that proof operation. This is a sufficient fixed-body/cofinal
  certificate, not arbitrary-body entry or LRC(14).
source: root + entry_scout + repo_connection_scout / LRC14 continuation session, 2026-09-01
depends_on: []
related:
  - THM-2047-phase-height-toric-arrangement-for-lrc
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-4151-scale-sensitive-first-window-odd-tail-lrc14-transfer
audit: >
  PASS / ACCEPT. Two independent symbolic audits checked the circle cut,
  endpoint-wall reduction, divisibility by 14, strict danger convention,
  two physical lifts, tooth containment, endpoint orientation, equality,
  min/max direction, parity collisions, and overlap with THM-4052/4151.
  A separate exact hostile attains equality in the rejected nonstrict
  additive two-endpoint bound.
---

# THM-4331 -- safe-component endpoint-denominator odd-wall escape

**PROVED ELEMENTARY COMPONENT-WALL CERTIFICATE + INDEPENDENTLY AUDITED;
LRC(14) REMAINS OPEN.**

## 1. Statement

For a nonempty finite set `H` of positive integers, put

```text
G_H={y in R/Z:min_(h in H)||hy||>=1/14}.               (1)
```

Let `I` be a positive-length connected component of `G_H`. Since `G_H`
misses a neighborhood of zero, it has a unique nonwrapping closed lift

```text
I=[L,R] subset (0,1),              W=R-L>0.             (2)
```

Write its endpoints in lowest terms as

```text
L=A_L/Q_L,             R=A_R/Q_R,                       (3)
```

and put `Q_*=min(Q_L,Q_R)`.

> **Endpoint-denominator escape theorem.** Let `a<b` be distinct positive
> odd integers. If
>
> ```text
> bW >= 2/7-1/Q_*,                                      (4)
> ```
>
> then there is an `x in R/Z` such that
>
> ```text
> min_(v in 2H union {a,b})||vx|| >= 1/14.              (5)
> ```

There is also a stronger, genuinely two-endpoint sufficient condition:

```text
bW > 2/7-1/Q_L-1/Q_R.                                  (6)
```

The inequality in `(4)` is nonstrict. The inequality in `(6)` is strict;
Section 4 gives an exact equality containment showing that this distinction
cannot be erased from the component-in-tooth argument.

## 2. Endpoint denominators retain the wall address

Every endpoint of `I` is a boundary wall of one of the body constraints.
Thus it has the unreduced form

```text
(14k+-1)/(14h),                  h in H.                (7)
```

The numerator `14k+-1` is coprime to `14`. Reduction can therefore cancel
only a divisor of `h`, leaving

```text
14 | Q_L                 and                 14 | Q_R.  (8)
```

This is the coordinate lost by a width-only component summary.

## 3. Two lifts and the open-tooth proof

For every `y in I`, both physical half-lifts

```text
x_0=y/2,                         x_1=(y+1)/2             (9)
```

are safe for `2H`. For an odd integer `r`, the two values `rx_0,rx_1`
differ by `1/2` modulo one, so `r` can be strictly dangerous on at most one
of the two lifts. Moreover, it is dangerous on one lift exactly when

```text
||ry||<1/7.                                             (10)
```

Suppose `(5)` were false. At every `y in I`, the two tails would have to
kill complementary lifts. In particular the larger tail `b` would be
dangerous on one lift throughout `I`, so the connected closed interval `I`
would lie in one open tooth, for some `k in Z` after choosing the compatible
real lift of the possibly circle-wrapping tooth,

```text
(u,v)=((7k-1)/(7b),(7k+1)/(7b)),       v-u=2/(7b).      (11)
```

Because `7|Q_L`, the positive left gap is

```text
L-u=[bA_L-(Q_L/7)(7k-1)]/(bQ_L) >= 1/(bQ_L).           (12)
```

The numerator is an integer; in fact it is odd because `A_L,b` are odd and
`Q_L/7` is even. Applying the same calculation at the right wall gives

```text
v-R >= 1/(bQ_R).                                       (13)
```

Since both endpoint gaps are positive,

```text
bW < 2/7-1/Q_L,
bW < 2/7-1/Q_R.                                        (14)
```

Either reversed nonstrict inequality contradicts hypothetical failure.
Their disjunction is exactly `(4)`, because

```text
min(2/7-1/Q_L,2/7-1/Q_R)=2/7-1/min(Q_L,Q_R).           (15)
```

Keeping both lower bounds in `(12)--(13)` instead gives

```text
bW <= 2/7-1/Q_L-1/Q_R,                                 (16)
```

so the strict reverse inequality `(6)` also proves escape. This proves both
certificates. **QED.**

## 4. Exact equality boundary for the additive certificate

Take

```text
H={3,10},
I=[43/140,13/42] is a positive-length connected component of G_H,
W=1/420,
(Q_L,Q_R)=(140,42),
b=107.                                                  (17)
```

The interval is exactly the intersection of the relevant safe bands of
runners `10` and `3`. It lies strictly inside the `b`-tooth

```text
(230/749,232/749),                                      (18)
```

with endpoint gaps

```text
43/140-230/749=1/(107*140),
232/749-13/42=1/(107*42).                               (19)
```

Consequently

```text
bW=107/420=2/7-1/140-1/42.                              (20)
```

Thus equality in `(16)` is compatible with strict containment in an open
tooth. This is a hostile to replacing `>` by `>=` in `(6)`; it is not
claimed to be an unsafe full row.

For comparison, the first safe component of `H={1,...,11}` is

```text
[1/14,13/154],                    W=1/77.               (21)
```

Condition `(4)` therefore closes every odd `b>=17`, independently of the
smaller odd tail `a`. This recovers the large-tail odd-wall calculation of
THM-4151 from the full safe-component address rather than a prescribed first
window.

## 5. Scope and relation to prior mechanisms

- THM-4052 compares an inherited Lipschitz arc with a complete joint-spoil
  component. Its sharper pair-dependent term and this endpoint certificate
  are incomparable; the Lipschitz arc need not begin at a body wall.
- THM-4151 is recovered on its anchored first window when
  `L=1/(14m)` and `Q_L=14m`, but its separate small-tail endpoint lift is
  not part of this theorem.
- A positive-length component is required. Isolated threshold witnesses are
  not covered.
- The conditions are sufficient, not necessary. They do not give
  arbitrary-body entry, a uniform bound over all bodies, or LRC(14).
