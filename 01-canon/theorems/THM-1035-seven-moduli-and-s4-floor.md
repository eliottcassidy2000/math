---
id: THM-1035
title: THE EXISTENCE STEP IN SEVEN MODULI (kernel-pure) + THE S₄ QUADRUPLE FLOOR — the denominator sieve's twelve conditions q = 2…14 collapse to SEVEN, q ∈ {8,…,14}, which is the UNIQUE minimum covering set (the six-number window {9,…,14} misses precisely q = 8, since no multiple of 8 lies in 9…14); and the S₄ floor by iterated containment — μ₄ ≥ (2λ)⁴·∏(1 − aᵢ/(2λaᵢ₊₁)) — holds 150/150 on separated quadruples and degenerates to 0 on comparable ones, exactly parallel to the pair story (THM-1012 vs THM-1025)
status: the seven-moduli reduction is PROVED and kernel-pure (TournamentH7.LRCSevenModuli, four theorems at the standard axiom trio, first compile); the minimality is by exhaustive search over subsets of {2,…,14}; the S₄ separated floor is verified 150/150 with 0 violations, the comparable case needs the sawtooth analogue (open)
source: opus-2026-07-17-S360 (owner: do the S₄ quadruple floor, prove the six-number window's divisor count and close the existence step)
depends_on: [LonelyRunner.sieve_one_div (the denominator sieve), THM-1030 (the collapse = the sieve), THM-1012/THM-1025 (the pair floors this parallels), THM-1026 (the five-slot ledger)]
scripts: 04-computation/s4_quadruple_floor_opus_S360.py -> 05-knowledge/results/s4_quadruple_floor_opus_S360.out
---

# THM-1035 — seven moduli, and the quadruple floor

## The existence step: seven, not six

`sieve_one_div` gives a lonely time `1/q` whenever some `q ≤ 14` divides
no speed — naively twelve checks. They collapse to **seven**: every
`q ∈ [2,14]` divides an element of `{8,9,10,11,12,13,14}` (2∣8, 3∣9, 4∣8,
5∣10, 6∣12, 7∣14, the rest trivially), and a multiple of `q′` is a
multiple of `q`.

**The count is seven, not six** — worth stating precisely, since the
conjecture put to me was six. Exhaustive search over all subsets of
`{2,…,14}`: the minimum covering set has size **7** and is **uniquely**
`{8,…,14}`. The natural six-number window `{9,…,14}` covers every
modulus **except q = 8**, because no multiple of 8 lies in 9…14. Eight is
the obstruction, and it is the only one — a six-number window would work
if 8 were replaced by any of its multiples in range, and none exists.

Kernel-pure in `TournamentH7.LRCSevenModuli`: `seven_covers`,
`lonely_of_seven_moduli` (the existence step), `all_moduli_of_seven` (the
reduction), `counterexample_needs_seven` (the contrapositive).

## The S₄ quadruple floor

Iterating THM-1012's containment — inside each `a`-arc count `b`-cells,
inside each `b`-arc count `c`-cells, and so on — gives

> μ₄ ≥ (2λ)⁴ · ∏ᵢ (1 − aᵢ/(2λ·aᵢ₊₁)).

Verified: **150/150** separated quadruples, zero violations, floor
positive throughout. On comparable quadruples the floor degenerates to
zero in 60/60 — **exactly the pair story repeating one level up**
(THM-1012 sharp only when separated; THM-1025's sawtooth needed for the
correlated regime). The S₄ slot of THM-1026's ledger therefore stands in
the same position pairs did before S357, and the same remedy is
indicated: a four-variable folded identity.

Encouraging datum for B₅, which wants S₄ **large**: on comparable
quadruples the exact μ₄ has median 0.000595 against the independence
value 0.000416 — about **1.43×** independence, i.e. comparable quadruples
are positively correlated, in the direction the ledger needs.

## Note on arXiv:2408.00183

Couvreur–Zémor, *Freiman's 3k−4 Theorem for Function Fields* — additive
combinatorics and algebraic geometry (Riemann–Roch spaces, genus bounds);
no covering systems, sieves, or LRC. The transferable pattern is
abstract: **a smallness hypothesis on a product/sum operation forces the
ambient object to be close to a highly structured model.** The repo
already runs that pattern additively (`LRCFreimanAP`,
`LRCFreimanAPBridge`, `LRCFreimanBurden`). The tension worth exploring:
the sieve constraint above is purely MULTIPLICATIVE (a counterexample
must carry multiples of 8…14), while Freiman structure is ADDITIVE — and
the dense core must satisfy both at once. A family forced to contain
multiples of seven fixed moduli while also having small additive doubling
is heavily over-determined; that over-determination, not either
constraint alone, looks like the exploitable lever.
