---
source: oracle-2026-06-01-S528
status: formalization milestone (machine-checked LRC cases in Lean 4 / Mathlib)
tags:
  - lonely-runner
  - lean
  - formalization
  - denominator-sieve
  - n3
  - axiom-clean
---

# Formalizing the Proven LRC Cases in Lean: All-`n` No-Multiple, `n=2`, and the `n=3` Sieve Cover

Per the request to *"formalize these and extend proofs to as high `n` as possible,"*
this session moved the small-`n` LRC results that the repo had proven on paper
(center-grid S522o; covering + mod-3 character S526) one notch toward
machine-checked truth. The honest unit of progress here is not a new theorem about
`n=14` — it is **what part of the small-`n` story a proof assistant will now certify
with no `sorry` and no project axioms.**

## What is now machine-checked (Lean 4 + Mathlib, axiom-clean)

All in `04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean`; every
declaration's `#print axioms` reports exactly `[propext, Classical.choice,
Quot.sound]` (the Mathlib foundations — no `sorry`, no LRC-specific axiom):

1. **`lonely_of_no_multiple` — LRC for *every* `n`, no-multiple case.**
   `0 < n` and `∀ i, ¬ n ∣ v i  ⟹  Lonely n v (1/n)`. This is the `t = 1/n`
   regular-`n`-gon vertex witness, a one-line corollary of the master sieve
   (`sieve_one_div`). It is the genuinely *all-`n`* unconditional fragment: for
   any number of runners, the only way to dodge the `n`-gon witness is to carry a
   speed divisible by `n` (THM-369 / the denominator sieve).

2. **`lonely_two` — the `n=2` case, fully.**
   Any nonzero speed `a` has a `2`-lonely time, witnessed by `t = 1/(2a)`
   (then `a·t = 1/2`, at distance `1/2` from every integer). Proved through the
   reusable `far_iff_fract` bridge (distance-to-every-integer ⟺ fractional part in
   the central band `[c, 1−c]`).

3. **`three_lonely_sieve_cover` — the `n=3` sieve coverage, fully.**
   If *either* no speed is divisible by `3` *or* no speed is divisible by `2`, then
   `v` has a `3`-lonely time (`t = 1/3` resp. `t = 1/2`). This is the machine-checked
   form of "the two unit witnesses settle everything except the residual."

## The honest boundary: what Lean does *not* yet certify

`three_lonely_sieve_cover` isolates — precisely, in checked code — the **residual
kernel of LRC@3**: speed families in which *some* speed is divisible by `3` **and**
*some* speed is divisible by `2`. After the standard scaling reduction to
`gcd = 1`, this is exactly the speeds-entangled-with-`6` tangle (e.g. `{2,3}`,
`{1,6}`, `{3,4}`). The analytic proofs already in canon close this:

- **S522o (center-grid):** put the larger speed at the far-center `1/2` via
  `t_k = (2k+1)/(2v₂)`; the `v₂` equally-spaced images of `v₁` have spacing
  `1/v₂ ≤ 1/3`, so one lands in the far-band `[1/3, 2/3]`.
- **S526 (harmonic covering):** `|B_a ∩ B_b| = 4/9 + (2/9)·χ(a)χ(b)/(ab)` with `χ`
  the Legendre symbol mod 3, giving `|SAFE| = 1/9 + (2/9)χ(a)χ(b)/(ab) ≥ 0`, tight
  only at `{1,2}`.

Neither argument is formalized yet: the center-grid needs the 1-D
pigeonhole/equally-spaced-points lemma; the harmonic covering needs the mod-3
Gauss/character-sum evaluation. **That gap is the whole point of recording this
milestone honestly** — the sieve layer is checked; the residual `6`-tangle (the
genuine content of even the smallest open-looking case) is still paper-only.

## Why this is the right formalization frontier

The analytic frontier (S526/S527) has shown that `n ≥ 4` and `n = 14` reduce to the
**same** object: a sign bound on a higher-resonance character sum / the discrepancy
of the cascade. Formalizing that is a research-grade undertaking. By contrast the
**denominator sieve is fully formalizable today**, and it already delivers a
non-trivial *all-`n`* statement (`lonely_of_no_multiple`) plus complete `n = 2` and
the `n = 3` reduction. So the machine-checked proof reaches exactly as far as the
elementary methodology does — every `n` in the no-multiple regime, and the two
unit-witness thresholds at `n = 3` — and stops, transparently, at the first place
the argument stops being elementary (the resonance correction). The Lean file is
now a faithful ledger of "proved by sieve" vs "needs the hard character bound."

## Verdict / next

- **Done, checked:** 11 axiom-clean declarations; the proven-cases section now
  formalizes the all-`n` no-multiple witness, `n = 2`, and the `n = 3` sieve cover.
- **Next formalization step (in increasing difficulty):** (a) the 1-D
  equally-spaced-points pigeonhole behind S522o's center-grid, which would close
  `n = 3` completely in Lean; (b) the mod-3 character evaluation behind S526; (c)
  only then the higher-resonance bound shared by `n ≥ 4` and `n = 14`.
- The Lean module remains catalogued under THM-369 (the sieve); these are its
  proven-cases extension.

## Artifacts
```
04-computation/lean/TournamentH7/TournamentH7/LonelyRunner.lean  (Cases section)
```
Related: S522o (center-grid n=3), S526 (harmonic covering n=3), THM-369 (sieve),
THM-386 (central-box grounding).
