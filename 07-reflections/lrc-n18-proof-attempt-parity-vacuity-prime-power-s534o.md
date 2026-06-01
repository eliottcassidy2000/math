---
source: oracle-2026-06-01-S534o
status: proof attempt (LRC n=18) — honest negative + new parity-degradation finding
tags:
  - lonely-runner
  - n18
  - proof-attempt
  - parity-law
  - prime-power
  - denominator-sieve
  - cascade
---

# An Honest n=18 Attempt: the Rigorous Sub-Cases, the Vacuity of Parity at n*=9, and the Wall

LRC(18) means 17 nonzero speeds, threshold `1/18` — **deeply open** (proven only for
`k ≤ 6` runners). This is a structured attempt with the full repo machinery: it
proves the rigorous sub-cases, characterizes the residual via the prime-power
generalization of the S533c three-channel parity law, runs the cascade heuristic,
and locates the wall. **n=18 is not proved.** The genuinely new result is *why* the
parity route dies at n=18.

## Rigorous: the sieve sub-cases and the q-covering necessity

n=18 has `n* = 9 = 3²`. The denominator sieve (THM-369) gives, for **every**
`q ∈ {2,…,18}`: if no speed is divisible by `q`, then `t = 1/q` is `18`-lonely
(because `1/q ≥ 1/18`). Verified numerically (`245/245` random 17-sets with no
multiple of 9 are lonely at `t = 1/9`). Hence:

> **A counterexample to LRC(18) must be `q`-covered for every `q ∈ {2,…,18}`** —
> it must contain a multiple of `9`, of `16`, and of `5,7,11,13,17`, etc.

This is real, unconditional progress, but it only carves away the sieve-reachable
boundary; the fully-covered residual is untouched.

## New: the n*=9 parity law is VACUOUS — and why (the prime-power degradation)

The S533c law generalizes: the full-support inside debt vanishes iff `0` is **not**
representable as `Σ cᵢ vᵢ (mod n*)` with each `cᵢ ∈ {1,…,n*−1}`. Computed firing
rates (fraction of primitive `(n−1)`-sets that are debt-free,
`lrc_n18_proof_attempt_s534.py`):

```
 n   n*   φ(n*)  units mod n*            parity fires
 4    2     1    {1}                     240/400   (~60%)
 6    3     2    {1,2}                    12/400   (~3%)
 8    4     2    {1,3}                     0/400
 9    9     6    {1,2,4,5,7,8}            0/400
14    7     6    {1,…,6}                  0/400
16    8     4    {1,3,5,7}                0/400
18    9     6    {1,2,4,5,7,8}            0/400   <-- our target
```

> **The parity obstruction is essentially unique to n=4** (`n* = 2`, the only case
> with a trivial unit group). It is already rare at n=6 and **vacuous for every
> primitive set at n ≥ 8** — including n=18.

**The 3-adic reason at n=18.** Mod `9`, the nonzero residues split by 3-adic
valuation: the **units** `{1,2,4,5,7,8}` (`v₃=0`) and the `v₃=1` residues `{3,6}`.
A `v₃=1` runner can already contribute `0 (mod 9)` (take `cᵢ = 3` or `6`, both
`≢ 0 mod 9`); a **unit** runner can contribute *any* nonzero residue `{1,…,8}`. So
with a single unit runner among 17, `0` is always reachable → the debt is always
present. Debt-freeness would require **all** runners divisible by 3 — which recurses
to the **n=6 mod-3 law one 3-adic level down**, and is non-primitive. This is the
prime-power recursion: `n* = 3²` stacks an `n=6`-type problem under an n=18-type one,
and the stacking is exactly what empties the single-level parity bit.

So the "free parity" that settled half of n=4 contributes **nothing** at n=18. The
difficulty is entirely in the higher-order, coupled resonance structure — the
inside debt / coupling gap (S529/S532b), now confirmed maximal.

## Heuristic: the cascade passes (margin 2.58) but is not rigorous

The cascade (S527) threshold `(n−1)·((n−2)/n)^{n−2}` is `1.12, 2.04, 2.31, 2.58`
for `n = 7, 14, 16, 18` — all `≥ 1`, and n=18 passes with the *largest* margin.
The feasible-arc carve confirms the shape: the **AP / regular 18-gon** `(1,…,17)`
carves to **exactly 0** (tight, the conjectured unique extremal), while sampled
primitive sets stay positive (`≈ 0.148`). But the per-step `(n−2)/n` contraction
needs a discrepancy bound (Erdős–Turán); its logarithmic factor is precisely what
prevents the heuristic from being a proof — the same gap that keeps the whole
conjecture open.

## Verdict

- **Proven (rigorous):** the no-multiple-of-`q` cases for all `q ≤ 18`; a
  counterexample must be `q`-covered for every `q ∈ {2,…,18}`.
- **New structural finding:** the parity obstruction degrades with the unit group of
  `n*` and is **vacuous at n=18** (potent only at n=4); the 3-adic filtration makes
  n=18's parity a stacked n=6 problem. This explains, from the parity side, why
  large/composite `n` is hard.
- **Heuristic:** cascade passes with margin 2.58; AP is the tight extremal.
- **Not proved:** the sieve-covered residual with the discrepancy/coupling gap — the
  same wall as n=14/16 and the general conjecture. **n=18 remains open.**

## What would actually close it
Either (i) a rigorous discrepancy bound turning the cascade's `(n−2)/n` shrink into a
theorem for the sieve-covered residual, or (ii) a *coupled* (not single-character)
control of the inside debt — the order-≥3 resonances the parity bit cannot see. The
parity route is now known to be insufficient beyond n=4; the live tools are the
cascade discrepancy and the coupling-gap structure.

## Artifacts
```
04-computation/lrc_n18_proof_attempt_s534.py
05-knowledge/results/lrc_n18_proof_attempt_s534.out
```
Related: S533c (mod-6 three-channel parity, HYP-2016), S527 (cascade), THM-369
(sieve), S529 (inside debt), S532b (coupling gap), S17 ("if 15 were prime"),
S515 (prior n=18 wacky attempts).
