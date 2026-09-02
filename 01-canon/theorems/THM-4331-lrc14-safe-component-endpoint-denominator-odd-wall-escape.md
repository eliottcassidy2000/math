---
id: THM-4331
title: "LRC(14) safe-component endpoint-denominator wall escapes"
status: >
  PROVED ELEMENTARY COMPONENT-WALL CERTIFICATE + VERIFIED-EXACT CONTROL +
  INDEPENDENTLY AUDITED;
  LRC(14) OPEN. A positive-length component of the 1/14-safe set with
  endpoint denominators Q_L,Q_R accepts every distinct positive odd-tail
  pair a<b after doubling whenever
  bW>=2/7-1/min(Q_L,Q_R). The stronger two-endpoint certificate is valid
  with strict > only. Without doubling, the same component accepts a new
  positive speed s when sW>=1/7-1/min(Q_L,Q_R), again with a strict
  two-endpoint improvement. Consequently mu(G_A)>1/189 implies that every
  s>=27(sum(A)-1) can be appended. Exact equality containments delimit both
  proof operations. These are sufficient certificates, not arbitrary-body
  entry or LRC(14).
source: root + entry_scout + repo_connection_scout + cofinal_transfer_scout + thm4333_referee / LRC14 continuation session, 2026-09-01
depends_on: []
related:
  - THM-2047-phase-height-toric-arrangement-for-lrc
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-4151-scale-sensitive-first-window-odd-tail-lrc14-transfer
  - THM-4333-lrc14-rank-three-surplus-and-cofinal-third-tail-completion
endpoint_control_script: 04-computation/lrc14_rank3_endpoint_appender_control_thm4331.py
endpoint_control_output: 05-knowledge/results/lrc14_rank3_endpoint_appender_control_thm4331.out
endpoint_control_independent_script: 04-computation/lrc14_rank3_endpoint_appender_control_thm4331_independent_audit.py
endpoint_control_independent_output: 05-knowledge/results/lrc14_rank3_endpoint_appender_control_thm4331_independent_audit.out
endpoint_control_script_sha256: 857d6e2607002741ac0d72063d18701237ad098b17e6e332c4217ddb546b66fe
endpoint_control_output_sha256: 9b68149de56b22c319df6a510f6fcbd16bb6ad7d3b5c218bf22ddcd31997eea8
endpoint_control_independent_script_sha256: 30a75fdd72f484cb4ae05ae0059854e45bf5d839230d9b5f5be7d8e6325e91fe
endpoint_control_independent_output_sha256: b25ff0942207df4ed22b0271df19f3d9f8c6fa09d4114b5cb73da119f278b310
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. Independent symbolic audits checked the circle cut,
  endpoint-wall reduction, divisibility by 14, strict danger convention,
  two physical lifts, tooth containment, endpoint orientation, equality,
  min/max direction, parity collisions, the single-appender reduction,
  component counting, denominator bounds, and overlap with THM-4052/4151.
  Separate exact hostiles attain equality in each rejected nonstrict
  additive two-endpoint bound. Global-wall and successive-intersection exact
  implementations independently reproduce the displayed rank-three control,
  normally and under optimized Python.
---

# THM-4331 -- safe-component endpoint-denominator wall escapes

**PROVED ELEMENTARY COMPONENT-WALL CERTIFICATE + VERIFIED-EXACT CONTROL +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

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

## 5. Single-speed component escape and global reserve corollary

The same endpoint arithmetic applies without doubling or a parity
restriction. Let `s` be any positive integer. Either condition

```text
sW >= 1/7-1/Q_*,                                      (22)
```

or

```text
sW > 1/7-1/Q_L-1/Q_R                                  (23)
```

implies

```text
I intersect G_s != empty,
G_(H union {s}) != empty.                              (24)
```

Indeed, hypothetical failure puts the connected closed lift `I` strictly
inside one compatible real lift of an open `s`-danger tooth

```text
(u,v)=((14k-1)/(14s),(14k+1)/(14s)),
v-u=1/(7s).                                            (25)
```

The divisibilities `14|Q_L,Q_R` from Section 2 give

```text
L-u >= 1/(sQ_L),             v-R >= 1/(sQ_R).          (26)
```

Both gaps are positive. Consequently failure forces

```text
sW < 1/7-1/Q_L,
sW < 1/7-1/Q_R,
sW <= 1/7-1/Q_L-1/Q_R.                                (27)
```

The contrapositives prove `(22)--(24)`. Notice again that the one-endpoint
certificate is nonstrict because the corresponding failure bounds are
strict, whereas the additive reverse must be strict.

The additive equality boundary is already visible in

```text
H={1,13},
I=[85/182,97/182],          W=6/91,
(Q_L,Q_R)=(182,182),        s=2.                       (28)
```

This is a positive component of `G_H` and lies strictly inside the tooth

```text
(13/28,15/28),
```

with both endpoint gaps equal to `1/364`. Hence

```text
sW=12/91=1/7-1/182-1/182.                              (29)
```

Thus equality in the additive bound can lose this entire component. It is
not a globally unsafe row: for example, `x=1/4` is safe for `{1,2,13}`.

There is a useful global corollary. Let `A` be any nonempty finite positive
integer set, put

```text
C_A=sum_(a in A)a,                                     (30)
```

and suppose

```text
mu(G_A)>1/189.                                         (31)
```

Then every positive integer

```text
s>=27(C_A-1)                                           (32)
```

satisfies `G_(A union {s})!=empty`. To prove this, note that the danger
teeth of `A` give at most `C_A` positive-length safe components. A largest
one therefore has

```text
W>1/(189C_A).                                          (33)
```

Each of its endpoint walls belongs to a speed at most `C_A`, so
`Q_L,Q_R<=14C_A` and

```text
1/Q_L+1/Q_R>=1/(7C_A).                                 (34)
```

For `(32)`, equations `(33)--(34)` give

```text
sW>27(C_A-1)/(189C_A)
   =1/7-1/(7C_A)
   >=1/7-1/Q_L-1/Q_R.
```

The strict additive certificate `(23)` completes the proof. **QED.**

## 6. Exact rank-three endpoint control

The rank-three hostile in THM-4333 also supplies a sharp component-address
control. Delete the auxiliary pool label `193` from its minimizing body and
put

```text
K={10,80,85,95,120,143,145,168},
H=K union {50,70},
S=2H union {1,9}
 ={1,9,20,100,140,160,170,190,240,286,290,336}.        (35)
```

The provenance body `K union {193}` has mask `011cd402`. Two independent
exact constructions of `G_S` give

```text
positive components=362,
mu(G_S)=64917367/577007200,                             (36)
```

and both contain the addressed component

```text
I=[227/476,937/1960],
W=39/33320,
(Q_L,Q_R)=(476,1960).                                  (37)
```

The left endpoint is owned only by speed `170`, and the right only by speed
`140`. The one-endpoint and additive thresholds are exactly

```text
(1/7-1/476)/W=4690/39,                                 (38)
(1/7-1/476-1/1960)/W=4673/39.                          (39)
```

Thus `(22)` certifies every integer `s>=121`, while the strict condition
`(23)` certifies every integer `s>=120`. In particular `s=120` is fresh and
gives an explicit thirteen-speed safe row. This component-specific bound is
far below both the original `27C=52,434` discrepancy bound and the global
reserve bound `27(C-1)=52,407`; it is a sufficient certificate, not the
minimal literal appended speed.

The primary artifact constructs all `3,882` reduced walls and classifies
their exact midpoints. The independent artifact successively intersects the
twelve closed safe-band lists and never constructs the global arrangement.
Normal and optimized runs of both paths reproduce `(35)--(39)` byte-for-byte.

## 7. Scope and relation to prior mechanisms

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
- The global reserve corollary uses only the widest component and the coarse
  denominator bound. Exact component addresses can give much smaller
  appender thresholds; certificate failure is not danger.
