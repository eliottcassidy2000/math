---
id: THM-1110
title: THE EXTENDED SIEVE LEMMA FOR q > 14 IS FALSE — an explicit primitive q=15 counterexample makes the classical threshold sharp. CORRECTION: the subsequent numerator count must be stratified by gcd(v,q); the advertised s<=11 conclusion at q=90 is valid only when every speed is a unit modulo 90, not merely nonzero. The primitive set {1,5,25,35} blocks every unit numerator modulo 90 although no speed is divisible by 90.
status: the q=15 counterexample, forbidden-window formula, and classical-threshold conclusion are PROVED. The original unstratified counting theorem is REFUTED; its corrected gcd-stratified union bound is proved below and exact-referee checked over 28,680 (q,v) rows.
source: opus-2026-07-17-S375 (owner: prove the extended sieve lemma for q > 14)
depends_on: [THM-1105 (which posed this as the open target and measured 34.8% failure), THM-1035/1040 (the classical lemma, whose threshold this shows sharp)]
scripts: 04-computation/extended_lemma_false_opus_S375.py, blocking_cost_opus_S375.py, extended_sieve_gcd_stratified_referee_codex_S67.py -> 05-knowledge/results/
---

# THM-1110 — the lemma is false, and exactly why

THM-1105 named this as the open target: "q divides no speed ⟹ some p/q is
lonely", for q > 14. It had already measured 34.8% failure in that regime.
This file settles it — the statement is **false** — and extracts the sharp
structure.

## Correction to the first counting corollary

The first version of this file then made a second, independent mistake.  It
counted only the dangerous residues coprime to `q` and applied that number to
every speed not divisible by `q`.  A nonunit speed can hit a larger gcd stratum
of the forbidden window.

At `q=90`, the unit group has size `24`.  The three residues `5,25,35`, all
of gcd `5` with `90`, have pairwise disjoint eight-element kill sets whose
union is the whole unit group.  Consequently

```text
{1,5,25,35}
```

is primitive, no speed is divisible by `90`, and nevertheless every reduced
numerator `p/90` is blocked.  Thus the claimed `s<=11` reach for arbitrary
nonzero speed residues was false already at four speeds.  The all-unit
specialization survives.

## The forbidden window, and why 14 is the boundary

p/q is lonely iff, for every speed v, v·p mod q ∉ W_q where

> **W_q = {r : min(r, q−r)·14 < q} = {0, ±1, …, ±⌊(q−1)/14⌋}**,
> **|W_q| = 2⌊(q−1)/14⌋ + 1**

so **|W_q| = 1 ⟺ q ≤ 14**. For q ≤ 14 the only forbidden residue is 0, so
"q divides no speed" *is* the full hypothesis — that is the classical
lemma, and the threshold is sharp rather than an artifact of its proof.
For q ≥ 15 the window widens (W₁₅ = {0,1,14}) and the hypothesis stops
being sufficient.

## The counterexample

> V = {11, 70, 77, 137, 144, 156, 175, 213, 226, 232, 246, 262, 281},
> primitive, and 15 divides no speed.

| p | blocked by | p | blocked by |
|---|---|---|---|
| 1 | v=226 (≡1) | 8 | v=77 (≡1) |
| 2 | v=232 (≡14) | 11 | v=11 (≡1) |
| 4 | v=11 (≡14) | 13 | v=232 (≡1) |
| 7 | v=77 (≡14) | 14 | v=226 (≡14) |

All eight numerators blocked, by just four speeds each killing exactly
|W₁₅|−1 = 2. The covering is tight.

## Why the coarsest union bound has no room at basic moduli

A speed v coprime to q kills the numerators p = v⁻¹w for w ∈ W_q\{0}, so

> #bad p ≤ 13·(|W_q| − 1) = 26⌊(q−1)/14⌋ ≈ 1.857 q > q > φ(q)

and this coarse union bound does not fire — for q = 15, 23, 29, 41, 71 the bad
count (26, 26, 52, 52, 130) exceeds φ(q) (8, 22, 28, 40, 70) every time.
This is the **same 13/7 > 1 obstruction** that killed S₁, the Bonferroni
ledger, and every union bound in this program, now appearing in the
purely arithmetic setting.

This rules out only the displayed sum of per-speed worst-case cardinalities.
It does not rule out overlap-sensitive, character-sum, or cross-modulus
counting arguments.

## The corrected gcd-stratified union bound

For a speed `v`, put `g=gcd(v,q)<q`.  Multiplication by `v/g` permutes
`U_(q/g)`, while reduction `U_q -> U_(q/g)` has constant fibre size
`phi(q)/phi(q/g)`.  Therefore `v` kills exactly

```text
K_(q,g) = [phi(q)/phi(q/g)]
          * #{r in W_q : gcd(r,q)=g}.                       (1)
```

Hence the correct theorem is:

> **THEOREM.** If speeds `v_i`, none divisible by `q`, satisfy
> `sum_i K_(q,gcd(v_i,q)) < phi(q)`, then some reduced `p/q` is lonely.

The unit-speed specialization has

```text
K_(q,1) = #{w in W_q\{0} : gcd(w,q)=1} = k_q,
```

so `s*k_q<phi(q)` is valid when **every speed is coprime to `q`**.

Within the original scan `q<200`, the unit-speed ratio `phi(q)/k_q` is largest
at `q=90`, where it equals `12`.  Thus eleven speeds **all coprime to 90**
are certified.  For arbitrary nonzero residues modulo 90, however, the largest
stratum has `K_(90,5)=8`, and `3*8=24=phi(90)` is sharp: the three gcd-five
residues above attain equality and cover.

The earlier statement that this alone closed the analogue for up to twelve
runners is withdrawn: an arbitrary family need not be coprime to the selected
modulus.  Formula (1), rather than the unit-only count, is the reusable result.

## The blocking cost

Without a multiple of `q`, the gcd-stratified union bound gives a blocking
cost at least

```text
ceil(phi(q) / max_(proper g|q) K_(q,g)).
```

At `q=90` this lower bound is `3` and is attained by `{5,25,35}`.  The
incoming unit-only blocking costs remain correct for unit speeds, but are not
lower bounds for arbitrary nonzero residues.  The genuinely global problem is
simultaneous blocking across many moduli, where the same integer lift must
serve all owner obligations; this is the shared-lift atlas of THM-1099.

The exact referee
`04-computation/extended_sieve_gcd_stratified_referee_codex_S67.py` checks
(1) over `28,680` pairs with `2<=q<=240` and freezes the `q=90`
counterexample.  Source/output SHA-256 values are
`3386325bfcea3fbc18d390d9b8a7b69c58ac40e2e097afa7a3f19ef8cb95bcbd` and
`b49c5237b61db4a1faf8facaebbae388bad62567f859b035a00096729c1318c8`.
