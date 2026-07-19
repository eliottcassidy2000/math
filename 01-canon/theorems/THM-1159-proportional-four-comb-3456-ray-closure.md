---
id: THM-1159
title: Exact shape erosion closes the proportional four-comb ray m(3,4,5,6)
status: PROVED -- exact closed-endpoint shape atlas plus THM-1148's uniform eight-core component floor. Every legal proportional prefix (3m,4m,5m,6m) leaves a closed interval of length 1/(42m), so an arbitrary fifth killer above 6m cannot cover it. This removes the first infinite method residual named in THM-1148; uniform r=5 remains OPEN
source: codex-2026-07-18-S74 frontier continuation
depends_on:
  - THM-1148  # exact eight-core component floor and identification of the residual ray
related:
  - THM-1097  # sharp one-comb discrepancy
  - THM-1133  # span-at-most-30 cone
  - THM-1137  # exact interval transfer Phi
  - THM-1140  # corrected four-comb reconnaissance
  - THM-1147  # exact two-comb endpoint law
script: 04-computation/lrc14_r5_3456_shape_erosion_referee_codex_S74.py
output: 05-knowledge/results/lrc14_r5_3456_shape_erosion_referee_codex_S74.out
---

# THM-1159 -- the first proportional residual has an exact erosion certificate

For a positive integer `k`, put

```text
D_k={t in R/Z : ||kt||<1/14}.
```

Thus `D_k` is open, while its safe complement is closed.  THM-1148 left

```text
(k1,k2,k3,k4)=m(3,4,5,6)                              (1)
```

as the first primitive ray outside its span, multiplier, `Q4`, and exact
`Phi` cones.  The ray is not a residual of the four-comb theorem itself.  A
direct primitive-shape calculation closes every legal scale on it.

## 1. Exact primitive-shape atlas

Let

```text
S={x in R/Z : ||3x||,||4x||,||5x||,||6x|| >= 1/14}
delta=1/42=1/(7*6).                                   (2)
```

Solving the wall equations

```text
kx=n +/- 1/14,                 k in {3,4,5,6},          (3)
```

and testing the resulting rational cells gives the following complete list
of positive-length closed components of `S`:

```text
[1/42,13/84]       [5/28,13/70]      [3/14,13/56]
[15/56,13/42]      [5/14,27/70]      [29/70,27/56]
[29/56,41/70]      [43/70,9/14]      [29/42,41/56]
[43/56,11/14]      [57/70,23/28]     [71/84,41/42].     (4)
```

There are 36 distinct wall events in `(0,1)` and no ties.  The wall labels at
the two ends of the components in (4) are, in order,

```text
3->6, 6->5, 5->4, 4->3, 3->5, 5->4,
4->5, 5->3, 3->4, 4->5, 5->6, 6->3.                   (5)
```

The closed convention matters: every displayed endpoint is safe at equality.
The finite verification also checks every wall separately, so (4) omits no
isolated safe wall.

For a closed safe component `[a,b]`, the starts `s` for which the oriented
interval `[s,s+delta]` stays in that component form `[a,b-delta]` when
`b-a>=delta`, and form the empty set otherwise.  Eroding (4) by `delta`
therefore gives the exact cyclic start set

```text
T = [1/42,11/84]
  U [15/56,2/7]
  U [5/14,38/105]
  U [29/70,11/24]
  U [29/56,59/105]
  U [43/70,13/21]
  U [29/42,17/24]
  U [71/84,20/21].                                    (6)
```

Its cyclic complementary gap word is

```text
(23/168, 1/14, 11/210, 5/84,
 11/210, 1/14, 23/168, 1/14).                         (7)
```

In particular the largest gap is

```text
G=23/168,             G+delta=27/168=9/56.             (8)
```

> **Shape-erosion lemma.**  Every real interval of length strictly greater
> than `9/56` contains a subinterval of length `1/42` whose image modulo one
> lies in `S`.

**Proof.**  If `I=[A,B]`, the allowed starts of a `delta`-interval inside
`I` form `[A,B-delta]`.  This interval has length greater than `G`.  After
reduction modulo one it cannot avoid `T`, because every connected component
of the cyclic complement of `T` has length at most `G` by (7).  Choose a lift
`s in [A,B-delta]` of an intersection point.  Definition (6) gives
`[s,s+delta] mod 1 subset S`.  The same argument is immediate if the allowed
start interval winds all the way around the circle.  This proves the lemma.
`square`

The strict `>G` avoids any ambiguity from an interval placed exactly in one
of the two largest closed-end gaps.  The produced `delta`-interval itself is
closed and includes safe wall endpoints.

## 2. Every legal scale clears the sharp four-comb target

Let `P` be an eight-element subset of `{1,...,12}`, put `M=max(P)`, let `m`
be a positive integer, and let `ell(P)` be the length of a longest component of

```text
C_P={t : ||pt||>=1/14 for every p in P}.               (9)
```

The exact 495-core atlas in THM-1148 proves

```text
ell(P)(13M+1)>=72/35,
ell(P)>=72/[35(13M+1)].                                (10)
```

Assume that the first proportional killer is legal:

```text
3m>13M.                                                (11)
```

The exact comparison needed by the shape-erosion lemma is

```text
      13M/3 > 5(13M+1)/64,
```

because clearing the positive denominator `192` leaves

```text
637M>15.                                               (12)
```

Thus (10)--(12) give

```text
m ell(P)
 >= m*72/[35(13M+1)]
  > (5(13M+1)/64)*72/[35(13M+1)]
  =9/56.                                               (13)
```

> **Proportional-ray theorem.**  Under (11), the set safe for the core `P`
> and the four killers `(3m,4m,5m,6m)` contains a closed interval of length
>
> ```text
> 1/(42m)=1/[7(6m)].                                   (14)
> ```

**Proof.**  Choose a longest closed core-safe component `J` and a real lift
`[A,B]` of it.  Under the phase map `x=mt`, its lift is the interval
`[mA,mB]`, whose length is strictly greater than `9/56` by (13).  Apply the
shape-erosion lemma to obtain

```text
[x0,x0+1/42] mod 1 subset S.                           (15)
```

Pulling (15) back by `x=mt` gives the closed interval

```text
I=[x0/m,(x0+1/42)/m] subset J.                         (16)
```

It is core-safe because it lies in `J`.  For `a in {3,4,5,6}` and `t=x/m`,

```text
||(am)t||=||ax||>=1/14,                                (17)
```

so it is also safe for all four proportional killers.  Its length is exactly
the value in (14).  `square`

For orientation, an eight-subset has `8<=M<=12`.  At the least legal integer
scales, the five exact phase-length lower bounds are

```text
M       8       9          10          11      12
m      35      40          44          48      53
m ell 24/35   288/413    3168/4585    24/35   3816/5495,
```

all far above `9/56`.  These rows are checks on (13), not a bounded-scale
substitute for it.

## 3. Appending the fifth killer

Let `k5>k4=6m`.  Every connected component of the open danger set `D_k5` is
an open tooth of length

```text
1/(7k5)<1/(7k4)=1/(42m)=|I|.                          (18)
```

If the connected interval `I` from (16) were contained in `D_k5`, it would
have to lie in one connected tooth, contradicting (18).  Hence some point of
`I` is safe for `k5` as well, and

```text
P union {3m,4m,5m,6m,k5}                              (19)
```

is `1/14`-lonely.

There is a slightly stronger endpoint fact: even an open tooth of length
exactly `|I|` cannot contain the closed interval `I`, since strict containment
at both endpoints would force a strict length inequality.  Thus the proof is
fully aligned with the actual open-danger/closed-safe convention; it does not
silently replace equality by a strict component estimate.

This proves the entire legal ray (1), including the `m>=53` tail first named
in THM-1148.  It does not prove the uniform four-comb theorem away from this
shape.

## 4. Tournament and alternate-carrier audit

Ordering the four runner vertices by speed gives the transitive tournament

```text
scores=(0,1,2,3), directed cycles=0,
four singleton SCCs, Hamiltonian paths=1.              (20)
```

It cannot distinguish this closed ray from any unresolved four-speed shape.
Ordering the 36 wall-crossing vertices after cutting the circle at zero is
again transitive, with scores `0,...,35`, no cycles, 36 singleton SCCs, and
one no-tie Hamiltonian path.  Ordering the twelve safe components gives the
same uninformative fingerprint on twelve vertices.

The proof-bearing object is instead the weighted cyclic interval complex

```text
(wall owner word,
 closed safe-component lengths,
 delta-eroded start components,
 cyclic complementary gap word,
 core phase-needle length m*ell(P)).                   (21)
```

The switch is erosion by the target length `delta`; the tie Hamiltonian path
is the grouped chronological wall path (there happen to be no ties for this
shape).  The maximum weighted complementary gap, not a pair orientation, is
the observable used in (8).

Candidate vertex sets explicitly challenged here are runners, combs, wall
crossings, wall endpoints, safe components, eroded start components, and
core components.  Runner and wall tournaments destroy all interval lengths.
Endpoint vertices without their paired components destroy erosion
eligibility.  Start-component vertices without the cyclic gap weights destroy
the covering radius.  Core-component vertices without the phase map destroy
the scale comparison (13).  Thus Tournament Analysis is useful telemetry,
but a labelled weighted cyclic carrier is minimal for this proof.

## 5. Exact replay

The dependency-free referee uses `fractions.Fraction` for every wall,
component, erosion, gap, and core inequality.  It reconstructs all 36 wall
events rather than trusting the displayed atlas; checks that they are
distinct; verifies all endpoints and the absence of isolated safe walls;
reconstructs (4), (6), and (7); and checks the least legal scale for each
`M=8,...,12`.  Ordinary and optimized runs are required to be byte-identical.

```text
04-computation/lrc14_r5_3456_shape_erosion_referee_codex_S74.py
SHA-256 cf21433cbb618d26805e672355a2c37e0c2684c6d4bcf05b55390f97efca2f3f

05-knowledge/results/lrc14_r5_3456_shape_erosion_referee_codex_S74.out
SHA-256 73d1722998bdbce7b2986d4e3cf7ba4fdcbd135d87a24c09df0e80dac4aa7d33
```
