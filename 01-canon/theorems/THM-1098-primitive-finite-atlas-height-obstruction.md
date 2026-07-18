---
id: THM-1098
title: Primitive finite-atlas obstruction with explicit lonely witnesses and logarithmic height cost
status: PROVED (elementary core; logarithmic asymptotic uses the prime number theorem); deterministic exact audit stored
date: 2026-07-18
session: codex-2026-07-18-S67 (finite-sieve frontier audit)
depends_on:
  - THM-566
related:
  - THM-1100
  - THM-1105
  - HYP-2052
  - HYP-2865
  - HYP-2876
results:
  - 05-knowledge/results/lrc14_primitive_finite_atlas_obstruction_codex_s67.out
---

# THM-1098: primitive finite atlases fail, quantitatively

Write

```text
C = {1,2,...,11,13}
```

and call a rational point `t` lonely for a speed set `S` when

```text
||v t|| >= 1/14  for every v in S.
```

The point of this theorem is that primitivity removes a **common** dilation,
but does not control the divisibility carried by one coordinate.

## Theorem A: every finite rational atlas is killed by an actually lonely row

Let `A` be any finite subset of `Q/Z`, and let `D_A` be the least common
multiple of the reduced denominators of its elements.  Choose an integer

```text
M >= 3731,       lcm(84,D_A) | M,
S = C union {M}.
```

Then:

1. `S` is a primitive `13`-speed set;
2. `S` covers every modulus `2,...,14` in the THM-523 sense;
3. no point of `A` is lonely for `S`;
4. nevertheless, `S` has an explicit lonely rational point whose reduced
   denominator is at most `2M`.

### Proof

Primitivity follows from `1 in C`.  The speeds `2,...,11,13` cover their own
moduli, while `84 | M` supplies multiples of `12` and `14`.

If `p/q in A` is reduced, then `q | D_A | M`, so

```text
M p/q is an integer.
```

The `M`-runner is at the observer, and `p/q` is not lonely.  This proves the
finite-atlas obstruction.

It remains to show that the constructed row really is lonely.  Put

```text
t0 = 17/41.
```

For `v in C`, the distances `41 ||v t0||` are

```text
17, 7, 10, 14, 3, 20, 4, 13, 11, 6, 18, 16,
```

in the order `v=1,...,11,13`.  Hence

```text
min_{v in C} ||v t0|| = 3/41 = 1/14 + 1/574.
```

Choose an odd integer `a` with

```text
|a - 2M t0| <= 1
```

(consecutive odd integers are two apart), and set `t=a/(2M)`.  Then

```text
||M t|| = ||a/2|| = 1/2.
```

Since distance to the nearest integer is `1`-Lipschitz, for every `v in C`,

```text
||v t||
  >= ||v t0|| - v |t-t0|
  >= 3/41 - 13/(2M)
  >= 3/41 - 1/574
   = 1/14.
```

The penultimate inequality is exactly `M >= 3731`.  Thus `t` is lonely, and
its reduced denominator is at most `2M`.  This proves all four claims.  `square`

The construction strengthens THM-566 in one useful respect: the obstructing
rows are not merely sieve-blind.  They come with a direct rational loneliness
certificate beyond the killed atlas.

## Corollary B: the least denominator has an unavoidable logarithmic height cost

For `B >= 14`, let

```text
L_B = lcm(1,2,...,B),
S_B = C union {L_B},
q_min(S_B) = least reduced denominator of a lonely rational point.
```

Theorem A applies because `84 | L_B` and `L_B >= L_14 = 360360`.  Every
`q <= B` divides `L_B`, so the last runner kills every rational with reduced
denominator at most `B`.  The explicit witness in Theorem A gives

```text
B < q_min(S_B) <= 2 L_B.                         (1)
```

This is an exact statement, requiring no asymptotics.  To express its height
cost, recall the Chebyshev function

```text
psi(B) = sum_{p^k <= B} log p = log lcm(1,...,B).
```

Thus `log(max S_B)=psi(B)`, and the prime number theorem
`psi(B) ~ B` turns (1) into

```text
liminf_{B -> infinity}
  q_min(S_B) / log(max S_B) >= 1.                (2)
```

Consequently, any height-dependent bound valid for all primitive covering
rows must exceed `B` at height `L_B`.  In particular it cannot be
`o(log H)`; a proposed asymptotic bound `C log H` needs `C >= 1` when natural
logarithms are used.  This does not supply an upper bound of that order: it is
the sharp necessary growth forced by divisor loading.

The same exact audit separates this obstruction from THM-1105's sampled `q0`
position law.  At `B=23`, the prefix row has first unblocked modulus `q0=25`,
but its least lonely denominator is `53` (witness `22/53`).  Thus the exact
excess `q_min-q0` is `28`, already larger than the incoming adversarial-search
record `19`.  This is not an unboundedness theorem for the excess; it is a
guardrail that `q0` remains a useful first address rather than a closure law.

## Corollary C: rational address and component thickness are independent

If a family contains speed `M`, every connected component of its lonely set
is contained in a safe tooth of that runner.  Those teeth have length

```text
(1 - 2/14)/M = 6/(7M).                           (3)
```

Now set

```text
M_k = 84(41k+1),
T_k = C union {M_k},       k >= 0.
```

Each `T_k` is primitive and covering.  Moreover `M_k == 84 (mod 41)`, and the
residue table above together with

```text
min(17*84 mod 41, 41-(17*84 mod 41)) = 7
```

shows that the *same point* `17/41` is lonely for every `T_k`.  Yet (3) gives

```text
largest lonely-component length <= 6/(7 M_k) -> 0.
```

So the slogan from THM-1100, "the gaps are situated, not large", is rigorous:
a fixed low-denominator address can persist while the component containing it
becomes arbitrarily thin.  Denominator address and geometric thickness are
different coordinates of the object.

## Consequence for the LRC(14) sieve frontier

Define the denominator-loading depth

```text
kappa(S) = max {B : for every q <= B, some v in S is divisible by q}.
```

This is `q0(S)-1` in THM-1105's notation, where `q0` is the first modulus
dividing no speed.

Whenever `q_min(S)` exists,

```text
q_min(S) > kappa(S).
```

The quantity `kappa` is unbounded even on primitive covering sets, by the
family `S_B`.  Therefore a rational-point sieve can only be an adaptive or
height-dependent atlas.  It must first delete the exact denominator packets
killed by divisibility and then prove that one of the remaining residue bands
is live.  For `q <= 14`, absence of a divisible speed is sufficient; for
`q > 14`, it is only necessary, and the genuine residual is the covering of
the numerator group by the thirteen danger bands.

## Assumption challenge and tournament view

The failed implicit assumption was that tournament vertices should be runners
modulo their common gcd.  That quotient preserves common dilation but destroys
the decisive fact that one runner can carry an arbitrarily large finite ideal
of denominator divisors.

A predicate-preserving alternative uses vertices

```text
reduced denominator packets, exact rational addresses, safe teeth,
divisor-loaded coordinates, and residue-cover proof obligations.
```

The pairwise observable is exact implication for the lonely-point predicate:
"kills this address", "bounds this component", or "supplies a witness".  With
the proof-strength gauge (exact obstruction before sampled kill rate), the tie
path is

```text
single-coordinate divisor load
  > common-gcd normalization
  > fixed-denominator atlas
  > random kill-rate evidence.
```

This quotient preserves the exact rational witness predicate and discards the
irrelevant labels of the twelve fixed core runners.  Its lesson is structural,
not numerical: primitivity is the wrong normalization for denominator
complexity.

## Deterministic audit

`04-computation/lrc14_primitive_finite_atlas_obstruction_codex_s67.py`
checks the residue margin, the explicit half-integer witnesses, finite-atlas
killing, prefix blindness, exact least denominators on representative rows,
and the fixed-address thin-component family using integer arithmetic only.

Frozen SHA-256 values:

```text
script  3ec85395d0b146b860f06c9347eedf0551a14d4b4b0d793d37e1a961d641060b
output  01a1449248749097a25499b627ad471e8ec6544c3d31829b1a382410fa7f5882
```
