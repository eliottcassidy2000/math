---
source: opus-2026-06-02-S562 (remote-control)
status: METHOD/SYNTHESIS — multi-sieve implemented; three assumptions challenged with data; the pinch (pair-sum) sieve is the complete, apex-free one
tags: [LRC, n14, multi-sieve, recursive, pinch, pair-sum, apex, even-fold, S551, S554, S557, S561]
---

# Multi-sieve, recursively: the apex dissolves, and the natural moduli are pair-sums

**Prompt (user):** implement multi-sieve concepts; think recursively; challenge
assumptions.

A config is *caught* by a sieve if it finds a witness `t` with `‖v_i t‖ ≥ 1/n`
∀i. I implemented three sub-sieves and combined them, then attacked three
assumptions I had been carrying.

## The three sub-sieves (`lrc_multi_sieve_recursive_s562.py`)

- **(D) Division** — `t = a/m`, `m ∈ {2,…,M}`. The standard sieve (= my S561 / the
  literature base).
- **(P) Pinch** — `t = a/(v_i+v_j)` over pairs. By S557 the *optimal* witness sits
  at a **pair-sum** modulus, so the natural sieve moduli are pair-sums.
- **(F) Recursive even-fold** — even speeds `2u` fold to `{u}`; find a good time `s`
  for the fold (by D/P/F recursively), take a doubling-preimage `t∈{s/2, s/2+½}`
  that is also odd-safe. Recurses on the 2-adic structure (n=14 → n=7 → …).

**Composition over 329 random+loaded configs:**

| sieve | caught |
|---|---|
| div `{2..14}` | 88.8% |
| recursive fold | 86.0% |
| div `{2..14}` ∪ fold | **95.1%** |
| div `{2..28}` | 100% (this sample) |
| **PINCH (pair-sum moduli)** | **100% — COMPLETE** |

## Assumption #1 — "the apex (multiple of 2q) is a hard obstruction." FALSE.

The apex was stuck only because I used **one** modulus. Of 2366 n=14 configs
containing a multiple of 14 (apex-stuck at `m=14`):

- caught by `m=14` alone: **0%** (the multiple of 14 sits at the observer at every
  `a/14`);
- caught by some `m ∈ {2,…,13}`: **91%** (e.g. rescued at `m=9` or `m=13`).

A different modulus has a different stuck runner; a config must be stuck at **all**
sieve moduli to survive. **Multi-sieving has no apex** — the S559/S561 "apex
obstruction" is an artifact of the single-corrector (polynomial-method) mechanism,
not of loneliness. This sharpens the program: don't fight the apex, change modulus.

## Assumption #2 — "sieve moduli should be small integers `{2,…,M}`." Wrong scale.

Small-integer division is provably *incomplete* at any finite `M` (S551: loaded
configs need ever-larger `M`; the residual's min-modulus reached 23+ here and is
unbounded in general). But the **pinch sieve catches 100%** — every config at its
own optimal time `a/(v_a+v_b)`, because `M(S)=r/s` with `s` a reduced pair-sum
(S557). So:

> The natural LRC sieve moduli are the **`O(k²)` pair-sums `v_i+v_j`**, not the
> integers `{1,…,M}` and not the single `c=(k+1)` lift. The pinch sieve is
> complete with **bounded modulus *count*** (one per pair) and has **no apex**.

(Honest: "complete" here means it evaluates `M(S)` exactly, so of course it catches
everything — LRC is true. The content is the *reframe*: pair-sum moduli, apex-free,
finite count. For the residue finite-check this says: sieve at **pair-sum residues
mod p**, not at the `c=(k+1)` base.)

## Assumption #3 — "the sieve is flat." It recurses.

The even-fold sieve (F) catches 86% by reducing n=14 to the **n=7** problem on the
halved even speeds, then lifting through a doubling-preimage. Its *only* miss is the
**odd-split** residual (S554) — exactly the loaded/apex configs — which the pinch
sieve covers. The recursion bottoms out at all-odd (witness `t=½`). So the multi-
sieve is genuinely recursive: 2-adic fold (F) × pair-sum pinch (P) × division (D).

## Synthesis: a recursive multi-sieve with the right primitives

```
catch(V):
  if division{2..n} or recursive-fold(V):   return  (cheap: ~95%)
  else:                                       return pinch(V)   (pair-sum, complete, no apex)
```

The cheap recursive tiers (division + 2-adic fold) clear ~95% at low cost; the
**pinch sieve** — pair-sum moduli, the structurally-correct primitive — clears the
rest and never sees an apex. The three obstructions I had been treating as
fundamental (apex, modulus unboundedness, flatness) were all **single-primitive
artifacts** that dissolve under the multi-sieve.

## What this changes for the n=14 attack

- **Drop the apex as a target.** It is modulus-specific; the pinch/multi-sieve has
  none. (Refines S559/S561 honestly.)
- **Sieve at pair-sums, not the `c=(k+1)` lift.** Concrete proposal for the finite-
  check: replace/augment the base lift with a **pair-sum-residue sieve mod p**
  (`t ≡ a/((v_i+v_j) mod p)`). The optimal witness lives there (S557), the count is
  `O(k²)`, and there is no zero-divisor apex.
- **Recurse 2-adically.** The fold reduces n=14 to proven LRC(13) on the even half
  (S558 correction) for ~86% before any expensive step.

**Honest scope:** these are exact witness-finders measured on samples, not an
end-to-end finite-check; "pinch complete" = computes `M(S)`. The value is the
de-mystification (no apex) and the modulus reframe (pair-sums), both implementable
in the real pipeline.

**Artifacts:** `04-computation/lrc_multi_sieve_recursive_s562.py` (+`.out`).
Builds on S551 (sieve incompleteness), S554 (even-fold), S557 (pinch radius),
S559/S561 (apex, two-tier CRT). New: **HYP-2075**.
