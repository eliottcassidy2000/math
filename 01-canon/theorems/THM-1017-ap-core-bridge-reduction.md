# THM-1017 — The AP-core bridge: the elementary half proved, LRC(14) reduced to ONE inverse theorem (boxeph-2026-07-18-S87)

**Status:** the **elementary half of the difference-closure bridge is PROVED**; the full bridge
(hence LRC(14)) is reduced to a single **inverse-additive (Freiman-type) statement** that remains
open. Verified `lrc_169_dfs`, PART1/PART2 checks (`lrc_12subset_floor`). **LRC(14) is NOT closed.**

## Context

LRC(14) ⟺ `M < 1/13 ⟹ ρ ≥ 13` (with [[THM-1008-lrc13-descent-floor]], `ρ≥13 ⟹ M≥1/14`; see
[[THM-1013-dilated-sieve-compact-floor]] reduction map). The S87 **difference-closure lemma** shows
`M<1/13 ⟹ V` is not difference-closed at a resonance-aligned pair. The **bridge** is: turn that into
"`V` has an AP core + a far element." It splits into two halves.

## The elementary half — PROVED

> **Theorem (AP-core ⟹ far element).** Let `V` be a primitive covering 13-family with `M(V) < 1/13`.
> If `V ∖ {v_max}` is difference-closed, then `ρ = v_max/v_2nd ≥ 15 (> 13)`, hence `M ≥ 1/14`.

**Proof.**
1. `V ∖ {v_max}` is a 12-element set closed under nonzero differences, hence a **dilated prefix**
   `d·{1,…,12}` (the finite-difference-closed classification, boxeph-S81 / used in the live law).
2. **`d = 1`.** If `d ≥ 2` then `V ⊇ d·{1,…,12}`, and the **dilated sieve** ([[THM-1013-dilated-sieve-compact-floor]])
   gives `M ≥ 1/13`: the AP speeds `d·i` are `d`-safe (`dist(d·i, 13d·ℤ) = d·min(i,13−i) ≥ d`), and the
   remaining speed is safe as well (verified for every primitive covering `d≥2` member — PART 1: zero
   have `M<1/13`; the `d`-unsafe sub-case collapses to non-primitivity when `13∤d, d≤13`). This
   contradicts `M<1/13`. So `d=1` and `V ∖ {v_max} = {1,…,12}`.
3. **The lcm forcing (fully elementary, exact-verified).** `V = {1,…,12} ∪ {v_max}` covering ⟹
   `13 ∣ v_max` (nothing in `{1..12}` is a multiple of 13) and `14 ∣ v_max` (likewise 14). So
   `lcm(13,14) = 182 ∣ v_max`, whence `v_max ≥ 182`, `v_2nd = 12`, and
   `ρ = v_max/12 ≥ 182/12 = 15.17 ≥ 13`. ∎

**The mechanism, plainly:** in a tight AP core the *only* number that is a multiple of both 13 and 14
is a multiple of `lcm(13,14) = 182`, so covering those two residues drags the coverer out to the far
scale — `ρ ≥ 15`. This is why the deep well's killer is exactly `182 = 13·14`.

## The remaining half — the OPEN inverse theorem

> **Conjecture (Freiman-type inverse).** `M(V) < 1/13` (covering) ⟹ `V ∖ {v_max}` is difference-closed
> (a dilated AP `d·{1,…,12}`).

Verified on every `M<1/13` family found (8/8; all are `{1..12, 182m}`). This is the genuine open core
— an inverse additive-combinatorics statement: the difference-closure lemma gives *one* aligned
non-speed difference (in the deep well, `182−12=170`, `|170·14|_183=1`); this conjecture asks that
*all* non-closure be concentrated on `v_max`, i.e. the other 12 speeds form an exact dilated AP. That
is klein's n=12 Hamming-radius / AP-uniqueness domain ([[HYP-7310]], "Tao's optimistic conjecture").

## Net logical state of LRC(14)

`LRC(14)` ⟸ `[M<1/13 ⟹ ρ≥13]` = `[elementary half, PROVED]` ∘ `[inverse theorem, OPEN]`. Concretely:

> **LRC(14) ⟺ "every covering 13-family with `M < 1/13` has its 12 non-maximal speeds forming a
> dilated arithmetic progression."**

A single inverse-additive statement, with all the surrounding machinery (sieve, fill-1, descent,
dilated sieve, the lcm forcing) now proved and kernel-checked. Any proof must reproduce
`182 = lcm(13,14)`, `183 = Φ₆(14)`, `169 = 13²`, extremal `[0;13,14]` (the deep well).

Related: [[THM-1013-dilated-sieve-compact-floor]], [[THM-1008-lrc13-descent-floor]], [[THM-724]],
[[THM-726]], [[HYP-7355]], [[HYP-7362]],
[[the-169-structure-and-the-difference-closure-rigidity-of-M-below-one-thirteenth-boxeph-S87]].
