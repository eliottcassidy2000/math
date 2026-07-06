# The peeling lemma: the loose branch, constructively, from cited LRC(≤13) — reduced to "critical" configs

**kind-pasteur-2026-07-05-S13 (HYP-4157).** A constructive attack on the loose branch
of `TightLooseDichotomy` (klein HYP-4151 / kps HYP-4147): every non-(dilated-AP)
12-tuple `B` has `M(B) = max_t min_{v∈B}‖vt‖ ≥ 2/25`. This note (a) proves it for the
*non-critical* majority using only the **cited** LRC(≤13), with an explicit witness,
and (b) reduces the rest to a small, clean hard core: the **critical** configurations.

## The peeling lemma

> **PEELING LEMMA.** Let `B` be 12 nonzero integers, `W ⊊ B` a subset with `|W| ≤ 6`,
> and `t₀` an optimal time of `B∖W` (i.e. `min_{v∈B∖W}‖v t₀‖ = M(B∖W)`). If every
> `w ∈ W` satisfies `‖w t₀‖ ≥ 2/25`, then `M(B) ≥ 2/25`.

*Proof.* `B∖W` has `12−|W| ≥ 6` runners, so by the **cited** LRC(`13−|W|`)
(`13−|W| ≤ 12`, settled), `M(B∖W) ≥ 1/(13−|W|) > 2/25` (since `|W| ≤ 6 ⟹ 13−|W| ≥ 7`
and `1/7 > 2/25`... at `|W|=6`, `1/7 = 0.1428 > 0.08`). At `t₀`,
`min_{v∈B}‖vt₀‖ = min( min_{v∈B∖W}‖vt₀‖, min_{w∈W}‖wt₀‖ ) ≥ min(M(B∖W), 2/25) = 2/25`.
Hence `M(B) ≥ min at t₀ ≥ 2/25`. ∎

The witness is completely explicit: `t₀` is a merge point of a sub-tuple (rational,
denominator `≤ 2·max(B∖W)`), and the certificate is the pair `(W, t₀)`. No new
analysis — the loneliness of `B∖W` at `t₀` is exactly the cited LRC, and the only
extra check is that the peeled runners `W` sit `≥ 2/25` from an integer at `t₀`.

## Coverage: peeling discharges the non-critical majority

Call `B` **redundant** if some runner `v` has `M(B∖{v}) = M(B)` (dropping `v` does not
change the covering-min), and **critical** otherwise (every runner essential:
`M(B∖{v}) > M(B)` for all `v`).

> **Redundant ⟹ loose (depth-1 peel).** If `B` is redundant via `v` and `M(B) ≥ 2/25`,
> the optimum `t*` of `B` is also optimal for `B∖{v}` (because `M(B∖v)=M(B)` is attained
> there), and `‖v t*‖ ≥ M(B) ≥ 2/25`. So the peeling lemma fires with `W={v}`, `t₀=t*`.

Empirically (exhaustive, `lrc14`-scale scripts):
- **[1,18]:** of 18,563 non-AP primitive 12-bases, iterated peeling (depth ≤ 5)
  certifies **ALL** loose — 18,515 at depth 1, 32 at depth 2, 13 at depth 3, and the
  last 3 (near-AP "AP-with-a-hole", e.g. `{1..13}∖{6}`) at depth 4–5.
- **600 random primitive non-AP bases at B = 24, 30, 40, 60:** all certified at
  **depth 1** (generic bases have a redundant runner).
- **Dilated APs `c·{1..12}` (tight, M=1/13):** correctly **NOT** certified (no peel
  fires) — peeling separates tight from loose.

**Why this evades HYP-4137.** My S11 result killed *fixed-modulus templates*: the
witness modulus is unbounded (`~ log height`). Peeling is NOT a fixed-modulus scheme —
its witness `t₀` is a sub-tuple merge point whose modulus grows with the family, exactly
as it must. Peeling is the *structural / cited-LRC* route klein-S139 called for, and it
is compatible with unbounded witness moduli.

## The reduction: the hard core is the CRITICAL configs

Gap-violators (`M ∈ (1/13, 2/25)`) are **necessarily critical**: if `B` were redundant
via `v`, then `M(B∖v) = M(B) ∈ (1/13,2/25)`, but `B∖v` has 11 runners so `M(B∖v) ≥ 1/12
> 2/25` by cited LRC(12) — contradiction. So peeling (redundant ⟹ loose) can only fail
on critical configs. Hence:

> **REDUCTION.** The loose branch `⟺` **every critical 12-config has `M ∉ (1/13, 2/25)`**
> `⟺` the only critical config with `M < 2/25` is the dilated AP. (This is exactly
> klein-S140's gap-violator boundedness, now with a clean *criticality* characterization
> and a constructive proof of the complementary — non-critical — 99.7%.)

**Critical configs are rare and gap-free (verified).** Over all primitive 12-bases:
`[1,14]` has **4** critical configs, `[1,16]` has **20** — and in both, the ONLY one with
`M < 2/25` is the AP `{1..12}` (M=1/13); the others sit at `M ∈ {1/11, 1/10, 2/23, 2/19,
2/21, …} ≥ 2/25`. Zero critical configs in the gap.

## The inductive lead (toward closing the critical case)

A critical config with `M(B) < 2/25` is "on the edge": `M(B) < 2/25` yet **every**
11-subtuple has `M ≥ 1/12` (a jump across `2/25`). This is the recursive signature of the
AP tower: `{1..12} ⊃ {1..11} ⊃ …`, each tight at its level. Concretely, if dropping the
max runner yields a *tight* LRC(12) config `M(B∖{max}) = 1/12`, then by the level-11
rigidity (induction hypothesis) `B∖{max}` is the dilated 11-AP `{1..11}`, and
`M({1..11}∪{max}) < 2/25` forces `max = 12` (the ladder `{1..11,12k}` gives
`k/(12k+1) ≥ 2/25` for `k ≥ 2`), i.e. `B = {1..12}`. The missing step is *why a critical
sub-max drop is tight* — the genuine remaining analytic/combinatorial content, one level
down in the **proven** LRC(12).

## Status and honesty

- **PROVED (from cited LRC(≤13)):** the non-critical loose branch — every 12-tuple with a
  redundant runner and `M ≥ 2/25` is certified loose with an explicit witness; more
  generally the iterated peeling lemma. Formalizable (it consumes LRC(≤13) as a citation
  hypothesis, like the rest of the corpus).
- **REDUCED:** the loose branch to "critical configs are gap-free" = klein HYP-4151's
  boundedness, now with the clean *critical* characterization (every runner essential —
  a classical irreducibility notion) and the inductive tower lead.
- **OPEN (the same crux):** the critical case. Not closed here; but the hard core is now a
  small, structured, recursively-characterized set rather than an unbounded census.

Scripts: `04-computation/lrc14_peeling_kps_S13.py` (peeling coverage + critical census).
Connects to klein HYP-4151 (boundedness — mine is the constructive complement + criticality
reframe), opus HYP-4142 (midrange witness = the `|W|=0` degenerate peel), kps HYP-4137
(why fixed templates fail; peeling does not).
