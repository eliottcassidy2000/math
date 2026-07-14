---
id: THM-753
title: The safe-peel reduction — every covering 13-set either reduces (via a chain of M-preserving "safe peels") to a set of ≤12 speeds, which is a SETTLED LRC(≤13) instance (M ≥ 1/13 > 1/14), OR is IRREDUCIBLE (no element is safe at its complement's tight point). The safe-peel lemma is elementary and RIGOROUS (M(S)=M(C) when the peeled speed is safe at the core's tight point); the reduction dichotomy is rigorous; and an adversarial search (5813 near-extremal + spread families) finds EVERY irreducible covering family is covered by an existing tile — shadow (THM-749), near-AP (kps THM-733/734/738/741, ≥9 in {1..14}), or loose (opus THM-745/746). ~98% of covering families reduce, so LRC(≤13) does most of the covering case.
status: PROVED (parts A, B — elementary). Part C EMPIRICAL (VERIFIED: 476/485 covering families reduce to ≤12; 5813 adversarial near-extremal families, ZERO irreducible-low-M-outside-all-tiles). The reduction dichotomy (B) is unconditional; whether a given family reduces or is irreducible is decided by the recursion. Part C (irreducible ⟹ a tile) is verified-not-proved — the honest residue.
source: mac-mini-2026-07-14-S103 (built on THM-751 safe-peel; the recursion + LRC(≤13))
depends_on:
  - THM-751   # aligned/safe far-element monotonicity — the safe branch is the M-preservation used here
  - THM-749   # the shadow tile (irreducible-shadow residual)
  - THM-734   # kps near-AP tile (irreducible-near-AP residual)
  - THM-745   # opus density floor (irreducible-loose residual)
related:
  - THM-366   # the t=1/q sieve (non-covering cores)
  - HYP-6660  # the tiling map (this makes LRC(≤13) the dominant reducer)
external: LRC(≤13) SETTLED (Rosenfeld / Trakulthongchai / Sungkawichai — cited, not re-audited).
---

# THM-753 — The safe-peel reduction to irreducible cores

**One line.** Most of the covering case is **LRC(≤13) in disguise**: peel off every speed that is
*non-binding* (safe at the rest's tight point) — each peel preserves `M` exactly — and a covering
13-set collapses to a ≤12-speed set, which is a settled `LRC(≤13)` instance with `M ≥ 1/13 > 1/14`.
Only the **irreducible** families (every speed binding) survive, and they fall into the known tiles.

## (A) The safe-peel lemma (elementary, rigorous)

Let `S` be a finite speed set, `v ∈ S`, `C = S\{v}`, and let `t₀` be a **tight point** of `C`
(`min_{u∈C} ‖u t₀‖ = M(C)`). Call `v` **safe at `t₀`** if `‖v t₀‖ ≥ M(C)`.

> **Lemma.** If `v` is safe at `t₀`, then `M(S) = M(C)`.

*Proof.* `M(S) ≤ M(C)` always (adding a speed only adds constraints). Conversely, at `t₀`,
`min_{u∈S} ‖u t₀‖ = min(\,M(C),\ ‖v t₀‖\,) = M(C)` since `‖v t₀‖ ≥ M(C)`. So `M(S) ≥ M(C)`. ∎

(The *aligned* case `v t₀ ∈ ℤ` is the complementary THM-751 branch, where `M` may drop but is bounded
below by `μ₀·wm/(wm+p_max)`; "safe" is the `M`-preserving branch used for the reduction.)

## (B) The reduction dichotomy (rigorous)

Call `S` **reducible** if some `v ∈ S` is safe at the tight point of `S\{v}`, else **irreducible**.

> **Every covering 13-set `S` satisfies exactly one:**
> **(i)** a chain of safe peels reduces `S` to a set `S'` of `≤ 12` speeds with `M(S) = M(S')`; then
> `S'` is an `LRC(≤13)` instance, so `M(S) = M(S') ≥ 1/13 > 1/14` (LRC(≤13), settled); **or**
> **(ii)** `S` is **irreducible** — every speed is binding at its complement's tight point.

*Proof.* Apply the lemma repeatedly: while a safe peel exists, remove it (`M` unchanged). Each removal
drops the size by one, so the process stops in `< 13` steps, at a set that is either size `≤ 12`
(case i) or has no safe peel (case ii, irreducible at whatever size it stopped). In case (i) the
terminal set has `≤ 12` distinct speeds, i.e. `≤ 12` moving runners, a settled `LRC(≤13)` instance
(`M ≥ 1/(\#speeds+1) ≥ 1/13`). ∎

So **the covering case reduces to the irreducible covering families.** (Empirically 476/485 = 98% of
covering families are reducible; only ~2% are irreducible.)

## (C) The irreducible residual is tiled (empirical)

The irreducible families are exactly the "all-binding" configurations — the structured extremals and
the genuinely spread sets. An adversarial search (mac-mini-S103, 5813 families biased toward the
covering-min extremals and toward spread) classified every irreducible covering family:

- **shadow-closable** (a `k≤13` shadow witness, THM-749) — the single-killer / tight extremals; or
- **near-AP** (`≥ 9` speeds in `{1,…,14}`, kps THM-733/734/738/741, `j ≤ 4`); or
- **loose** (`M ≥ 0.15`, high-diameter spread — opus's density floor THM-745/746).

**Zero** irreducible families were found that are low-`M` (`M < 0.12`), not shadow-closable, **and**
have `< 9` speeds in `{1,…,14}` — i.e. none escaping all tiles. (This corrected an initial false gap
`{1,3,4,7,8,9,10,11,12,13,51,59,182}`, `M=0.105`, which has `10` speeds in `{1..14}` and so is covered
by kps THM-738.)

## The assembled picture

| covering 13-set | route | status |
|---|---|---|
| **reducible** (~98%) | safe-peel chain → `≤12`-set → **LRC(≤13)** | PROVED (A,B) + cited LRC(≤13) |
| irreducible, shadow-closable | THM-749 (`k≤13` shadow), single-killer | PROVED |
| irreducible, `≥9` in `{1..14}` | kps THM-733/734/738/741 (`j≤4`) | PROVED (`j≤4`) |
| irreducible, loose (`M≥0.15`) | opus THM-745/746 density floor | PROVED (`W>W0`) + finite check |

So the covering case is **LRC(≤13) (via safe-peel reduction) ⊕ the four tiles for the irreducible
residual** — and the reduction shows the *bulk* of the covering case is not new mathematics but a
reduction to the already-settled lower cases.

## Honest scope

- **Proved:** the safe-peel lemma (A) and the reduction dichotomy (B) — elementary and rigorous. LRC(≤13)
  is cited (settled). Every reducible covering family has `M ≥ 1/13`.
- **Verified-not-proved:** that every *irreducible* covering family lands in a tile (C) — 5813-family
  adversarial search, zero escapees — and the "~98% reducible" frequency. Making the covering case a
  theorem now needs: a *proof* that irreducibles are tiled (not just the search), which is exactly the
  union of THM-749 (shadow) + kps near-AP + opus loose over the irreducible stratum. This is the residue
  of HYP-6660's tiling, now with LRC(≤13) discharging the reducible bulk.
- **A caveat on `M`-preservation:** the recursion uses the *tight point* of each core; the lemma is exact
  given a tight point, but a numeric tight-point finder can misclassify a borderline safe/aligned peel —
  a Lean/exact version should locate `t₀` rationally (the cores are small, `≤12` speeds).

*Artifacts:* `04-computation/lrc14_flexible_peel_macmini_S103.py`, `lrc14_full_recursion_test_macmini_S103.py`,
`lrc14_irreducible_residual_macmini_S103.py`, `lrc14_lowM_irreducible_search_macmini_S103.py` (+outs).
Credits: THM-751 (the safe branch), THM-749/734/745 (the tiles), THM-366 (the sieve), LRC(≤13) (settled).
