# Fibonacci is the covering-min's foil, not its lever — the anti-golden Eisenstein sibling

*klein-2026-07-04-S124 (HYP-4076). Owner: work the remaining core AND search back through the
repo's Fibonacci-connected work for creative relations. I did both: a repo-wide survey + fresh
computation. The finding is clean and honest — the covering-min's arithmetic is Fibonacci's
sign-flipped **sibling** (Eisenstein `x²−x+1` vs golden `x²−x−1`), and the covering-min sits at the
**anti-golden** continued-fraction pole, precisely *avoiding* the three-gap worst case that
Fibonacci/golden governs. So Fibonacci is the foil, not a lever on the argmax core.*

## Two exact facts (verified n=4,7,14; script)

**1. The recurrence sibling.** The covering-min witness `t* = n/Φ₆(n)` binds through the powers of
`n` mod `Φ₆(n) = n²−n+1`:
> `n^k ≡ 1, n, n−1, −1, −n, 1−n` (period 6), satisfying **`s_{k+2} = s_{k+1} − s_k`**.

That is the `x²−x+1` recurrence (roots = primitive 6th roots of unity, Eisenstein `ℤ[ω]`, Heegner
`−3`). **Fibonacci is `s_{k+2} = s_{k+1} + s_k`, the `x²−x−1` recurrence** (roots = golden `φ`,
`ℤ[φ]`, `√5`). The covering-min lives on the **exact sign-flip of Fibonacci** — the other of the two
"metallic-adjacent" quadratics `x²−x∓1`. Period 6 (Eisenstein) vs aperiodic (golden).

**2. The anti-golden position.** The witness continued fraction is
> `14/183 = [0; 13, 14]`  — general `n/Φ₆(n) = [0; n−1, n]`, and `1/M = Φ₆(n)/n = [n−1; n]`.

These are **large** partial quotients — the *fastest*-converging CF. Golden/Fibonacci ratios
(`5/8, 8/13, F_k/F_{k+1}`) are `[0; 1, 1, 1, …]` — the *slowest*-converging CF, the **three-gap
worst case**. The covering-min extremizer is the **anti-golden** extreme of the same Stern-Brocot /
Farey ladder on which the covering-min lives (`[0; n−1, k]`, "the covering-min lives on the Farey
ladder"). It is engineered to be as far from golden as a rung can be.

## Why this makes Fibonacci a foil, not a lever

The three-gap theorem's *difficulty* concentrates at golden/Fibonacci (slowest CF, largest
persisting gaps). One might hope to import that hardness. But the covering-min VALUE sits at the
**opposite pole** (anti-golden, fast CF), and its arithmetic is Eisenstein (`x²−x+1`), not golden
(`x²−x−1`). The pair-overlap kernel `K(a,b)` (three-gap/Stern-Brocot, mac-mini S75b) is *not*
uniquely maximized by Fibonacci ratios — the antipode column `b = n−1` (e.g. `b=13`) peaks it, with
non-Fibonacci `7/13, 6/13` tying `8/13`. So Fibonacci's slow-CF hardness does **not** transfer to
the covering-min bound. The genuine number-type pairing is **Eisenstein `√−3` × heptagon `√−7`
(cross `√21`)** — the deep well's `183 = 3·61` (Eisenstein) meets `14 = 2·7` (heptagon), per opus's
`the-13-comb-lever-is-the-eisenstein-resonance`. **`√5` never enters.**

## Repo survey verdict (all the Fibonacci threads)

- **Pisano tower** (THM-491): the LRC shell ramification `27→9→3` mirrors `π(p^k)=p^{k−1}π(p)`
  exactly — a beautiful *parallel* over the magma `(ℤ,+)`, but not load-bearing for the covering-min.
- **THM-486** (`π(24)=24`, `F₁₂=144=12²`): Zeckendorf companion; meets LRC "at 24" without a
  covering-min mechanism.
- **Golden exceptional points** (THM-224): golden eigenvalues live in *tournament transfer
  matrices*, not the covering-min.
- **Metallic spine spectrum** (opus): the SC-spine's blue spectrum is metallic (golden n=4, silver
  n=5) — the closest prior approach to the sibling comparison; my recurrence-sibling statement is the
  crisp form of that resonance.
- **Zeckendorf / Farey-Fibonacci** threads: meta-methodology (proof carriers, representation
  economies), analogy-level for the covering-min.

**All tangential to the argmax barrier.** The loose-U core (`M(2U∪2odd) ≥ 1/12`) is Eisenstein-
Farey-three-gap arithmetic; Fibonacci growth does not enter, and no Fibonacci lever tightens it.

## The transferable point

When a hard extremal problem has a continued-fraction / Stern-Brocot skeleton, the instinct is to
reach for the golden ratio (its most famous extremal point). But an extremum can live at the
*anti-golden* pole — the fastest-CF rung — precisely because the constraint (here: covering, i.e.,
divisibility) forces it there. The covering-min is the *maximally-commensurate* (period-6 Eisenstein)
configuration, the opposite of the *maximally-incommensurate* (golden) one. Fibonacci names the pole
the covering-min flees. Knowing which pole your extremum sits at tells you which arithmetic to bring:
here Eisenstein (`√−3`, order 6) and the heptagon (`√−7`), not the golden (`√5`). See
[[fixed-point-extremum-covering-not-transform]] — and [[the-13-comb-lever-is-the-eisenstein-resonance]].

## Status / links

- **New (exact):** the recurrence sibling `x²−x+1` (Eisenstein, covering-min) vs `x²−x−1`
  (Fibonacci); the anti-golden CF position `[0;n−1,n]`. Script:
  `04-computation/lrc14_fibonacci_eisenstein_duality_klein_S124.py` (+ `.out`). HYP-4076.
- **Verdict:** Fibonacci is tangential to the covering-min crux — the foil (opposite CF pole), not a
  lever. The argmax core stays Eisenstein-Farey-three-gap, LRC(14)-equivalent, OPEN.
- Threads: THM-486/491 (Pisano), THM-224 (golden), THM-387 (two-gap), the three-gap Stern-Brocot
  cap-kernel reflection, the-covering-min-lives-on-the-farey-ladder, opus metallic-spectrum &
  eisenstein-resonance reflections; kps HYP-4060, mac-mini THM-615/S37.
