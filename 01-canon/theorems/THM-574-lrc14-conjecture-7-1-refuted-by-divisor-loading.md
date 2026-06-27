---
id: THM-574
title: The paper's Conjecture 7.1 is false for k=13 as stated; divisor-loaded non-tight covering rows kill any proposed universal denominator
status: PROVED
date: 2026-06-27
session: codex-2026-06-27-S255
depends_on:
  - THM-566
related:
  - HYP-3088
  - HYP-3087
  - THM-573
  - THM-565
  - HYP-2072
results:
  - 05-knowledge/results/lrc14_conj71_refutation_normalized_arc_codex_s255.out
external:
  - arXiv:2604.23906 Conjecture 7.1
---

# THM-574: Conjecture 7.1 Refuted by Divisor Loading

## Statement

For `k=13`, Conjecture 7.1 of the paper *Eleven, twelve, and thirteen lonely
runners* is false as stated.

More explicitly, there is no constant `D` such that for every integer `d >= D`,
every non-tight primitive 13-tuple has a witness time in `(1/d)Z`.

## Construction

For an integer `B >= 6`, set

```text
L_B = lcm(1,2,...,B)
N_B = 84 L_B
S_B = {1,2,...,11,13,N_B}.
```

This is a primitive 13-tuple because it contains `1`.

## Denominator obstruction

For every `d <= B` and every numerator `a`, since `d | L_B`, we have

```text
N_B * a/d in Z.
```

Thus `||N_B a/d|| = 0 < 1/14`, so `a/d` is not a witness for `S_B`. Hence
`S_B` has no witness in `(1/d)Z` for any `d <= B`.

This is the THM-566 obstruction, specialized to the paper's exact denominator
quantifier.

## Non-tightness

The same tuple is not tight once `B >= 6`.  Since `12 | N_B`, the time

```text
t_B = 1/12 + 1/(2N_B)
```

makes the loaded runner exactly opposite the observer:

```text
N_B t_B = N_B/12 + 1/2 in Z + 1/2.
```

For the small speeds:

- if `1 <= i <= 5`, then `i t_B` is farther than `1/14` from an integer;
- if `i=6`, it is near `1/2`;
- if `7 <= i <= 11`, the critical distance to `1` is
  `(12-i)/12 - i/(2N_B)`, minimized at `i=11`, and this is `> 1/14`
  whenever `N_B > 462`;
- for `i=13`, the phase is `1/12 + 13/(2N_B)` modulo `1`, again `> 1/14`
  from an integer.

For `B >= 6`, `N_B >= 5040 > 462`, so all inequalities are strict. Thus
`S_B` is non-tight.

## Refuting the universal denominator

Given any proposed constant `D`, choose `B >= max(D,6)` and set
`d = max(D,6)`. Then `d >= D` but `d <= B`, so the denominator obstruction says
there is no witness in `(1/d)Z`, while the explicit time `t_B` says `S_B` is
non-tight. This contradicts Conjecture 7.1.

## Proof-route impact

This does not damage the LRC(14) witness route; it corrects its coordinate.
The false invariant is a uniform largest interval in the original time circle.
The viable invariant is the normalized slow/ruler-coordinate arc floor from
THM-565, after a large apex has been peeled. Divisor loading forces direct
time components to have length `O(1/N_B)`, but their apex-normalized width can
remain bounded.

This also corrects HYP-3088: the paper's Conjecture 7.1 is not literally the
project's witness route. The repaired bridge is:

```text
paper I(k,p,1) / c=2,7 lifts
  <-> THM-573 level-7 sieve + dyadic descent
  <-> normalized finite-ruler arc floor in the <=6-multiples-of-7 residual.
```

## Verification

`04-computation/lrc14_conj71_refutation_normalized_arc_codex_s255.py` checks
`B in {6,10,14,26,41,67}` exactly. It verifies:

- `gcd(S_B)=1`;
- the explicit witness margin is positive;
- no denominator `d <= B` has a grid witness in the tested rows;
- the direct largest time component shrinks with the loaded apex, e.g. for
  `B=6`, the largest component is `1/5880`.
