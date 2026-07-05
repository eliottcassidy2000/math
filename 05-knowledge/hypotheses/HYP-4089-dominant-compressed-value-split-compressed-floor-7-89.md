# HYP-4089 — The dominant/compressed split is a VALUE split; the compressed floor is 1/13 (the 12-runner LRC bound), attained by DILATED deep-wells — so `compressed ⟹ M ≥ 1/13` is the clean, tight `hcomp` target. [S45's "7/89" was WRONG.]

**Status:** VERIFIED (exact-M, >12k families incl. dilated + adversarial constructions). Empirical + partial proof route.
**Source:** mac-mini-2026-07-04-S45 (original, 7/89 — WRONG), **corrected mac-mini-2026-07-04-S46 (→ 1/13)**
**Type:** structural reduction of the covering-min / direct `hcomp` discharge target

> **⚠ S46 CORRECTION (supersedes the S45 body below).** S45 claimed the compressed floor is `7/89`.
> **That was WRONG** — it was computed over *small-range* families and MISSED **dilated** families. The
> dilated deep-well `c·{1,…,12} ∪ {182}` (`c≥3`, `13∤c`, `gcd(c,182)=1`) is **compressed** (`182 ≤ 13·12c`)
> and covering, with **`M = 1/13 = 0.076923` EXACTLY** — below `7/89 = 0.078652`. So:
> **the true compressed floor is `1/13`, not `7/89`.** Verified: `0` compressed covering families below `1/13`
> over **12,158** sampled/structured/dilated + adversarial large-`c` constructions (all `= 1/13`); `94` attain
> `1/13` exactly. **`1/13 > 14/183` still**, so the covering-min and "razor-thin only at the deep well" are
> INTACT. (klein-S129's "non-deep-well `≥ 7/89`" holds only within its *minimal-tightener* 509-family scope —
> dilated families are outside it and sit at `1/13`.)

## S47 refinements (mechanism firmed up; reliability caveat)
- **n=13 tight locus = dilated APs only** (verified, small range + structured): `M(W)=1/13 ⟺ W = c·{1,…,12}`.
  So the CRT tight-case argument (below) covers the *entire* tight locus.
- **The `1/13` floor is confirmed by exact local descent from structured seeds** (dilated deep-wells,
  `{1..11,24}∪{182}`, `{1..11,13,84}`, even-dilations `2·{1..12}∪{13}`, adversarial `157·{1..12}∪{18382}`) —
  all converge to `1/13`, nothing below. This is the *reliable* method; see caveat.
- **12-runner base spectrum near `1/13`:** `1/13` (dilated AP) → **`2/25`** (`{1,…,11,24}`, `= 2/(2n−1)`, the
  known second-smallest) → `1/12`. A transient S47 "gap up to `1/12`" was a **sampling artifact** (missed
  `{1..11,24}`); the real base-spectrum gap is `(1/13, 2/25)`. `{1..11,24}∪{182}` is compressed with `M=2/25 ≥
  1/13` (not below).
- **⚠ RELIABILITY (MISTAKE-102):** `compressed ⟹ M ≥ 1/13` is a **conjecture** — random-sampling verification
  is unreliable here (it misses the commensurate/dilated extremizers, exactly the families that matter). The
  support is: (a) the CRT free-rider argument for the tight-AP-base case; (b) exact descent from *structured*
  seeds; (c) direct exact-M on adversarial constructions. All consistent, none exhaustive. The covering-min
  `14/183` is unaffected regardless (`1/13 > 14/183`).

## The corrected statement (S46)
> **`compressed covering ⟹ M ≥ 1/13`** — TIGHT (dilated deep-wells attain `1/13`). This is the **12-runner
> LRC bound**: a compressed covering 13-family hides exactly as well as 12 runners. It is `> 1/14` (LRC target)
> AND `> 14/183` (covering-min), so it discharges kps's open `hcomp` leaf with a clean, structural margin.

## Why `1/13`, and the dominant/compressed = offset-forcer/free-rider dichotomy
A covering family has a killer `a` (covers `{13,14}`, so `13∣a`). At the base's max-min `t*` the killer's phase
is `a·t*`. Two regimes:
- **DOMINANT** (`a > 13·(2nd)`, e.g. deep well `{1..12,182}`, `t*=1/13`): `a t* = 182/13 = 14 ∈ ℤ`, so the
  killer sits at **0** — it **kills the base optimum**, forcing an offset → `M = 14/183 < 1/13`. *(Offset-forcer.
  Already DISCHARGED: kps HYP-4087 dominant peel — a giant runner's danger arcs are too thin to spoil the base.)*
- **COMPRESSED** (`a ≤ 13·(2nd)`, e.g. dilated `{3,…,36,182}`, base optimum `t*=1/39`): `a t* = 182/39 = 14/3
  ∉ ℤ`, `‖·‖ = 1/3 ≥ 1/13` — the killer is a **free rider**, safe at the base optimum → `M = M(base) = 1/13`.

So the `13×` line separates *offset-forcer* killers (`M ↓ 14/183`) from *free-rider* killers (`M = 1/13`): the
dominant/compressed split **is** the value split. All razor-thinness (`14/183`) lives on the CLOSED side.

## The peel route to `compressed ⟹ M ≥ 1/13` (partial proof)
Peel the max runner `v*`; `W = V∖{v*}` is 12 runners, `M(W) ≥ 1/13` (LRC(13)).
- **Loose base** (`M(W) > 1/13`): `W`'s good set has positive width; `v*` is safe on density `1−2/13 = 11/13`
  of it → a common safe spot (needs a width-vs-arc-spacing check; empirically always holds).
- **Tight base** (`M(W) = 1/13`, `W = c·{1,…,12}` a dilated AP): `W`'s optima are `{k/(13c) : 13∤k}` (many).
  `v*` is safe at some optimum **by CRT**: primitivity `gcd(V)=1` with `gcd(W)=c` forces `gcd(c,v*)=1`, so
  `v*/(13c)` has denominator with the factor `c` **coprime to `13`**; `v*`'s safety is indexed mod `c`, the
  base optima mod `13`, and `13 ⊥ c` gives a joint residue that is both a base optimum (`13∤k`) and `v*`-safe
  → `M(V) ≥ 1/13`. **Crucially `c>1`** here: `c=1` forces the deep well (`182>156`, dominant, excluded from
  compressed). *Verified against adversarial large-`c` killers (`157{1..12}∪{18382}`, `313{1..12}∪{44772}`,
  killer a near-multiple of `13c` unsafe at `j=1..12`): all give `M=1/13` exactly — the CRT optimum (e.g.
  `j=14`) saves it.*
- **Open in the route:** the loose-base width bound; the `n=13` tight locus (is it only dilated APs, or also a
  GW-analog?); multi-runner peels. So this is a proof *route*, not yet a proof.

## Consequence for the Lean endgame
kps's sole open leaf `hcomp` (compressed covering ⟹ lonely, i.e. `M ≥ 1/14`) is implied by the **clean tight
bound `M ≥ 1/13`** — a `1/13 vs 1/14` structural margin, NOT the razor `14/183`. If the peel route closes, it
discharges `hcomp` **uniformly from LRC(13)** (the compressed twin of opus's dominant peel) — no census, no
family-by-family ladders. That would close the covering side entirely: dominant peel (giant killer) + compressed
peel (free-rider killer), the two sides of the `13×` line, both from LRC(13).

## Resonance / links
`14/183 = [0;13,14]` (dominant, offset) vs `1/13 = [0;13]` (compressed, free-rider) — the offset `1/13 − 14/183
= 1/2379` is exactly the killer-offset unit (THM-618). Compressed = the "no-offset" world. See THM-618
(single-killer offset), HYP-4087 (dominant peel), HYP-4091 (hcomp), klein-S128/S129, THM-617.
Reflection: [[the-dominant-compressed-split-is-a-value-split]] (updated S46).

**Scripts:** `04-computation/compressed_floor_dilated_macmini_20260704.py` (+ `.out`),
`compressed_covering_min_macmini_20260704.py` (S45, superseded floor),
`05-knowledge/results/compressed_floor_1over13_adversarial_macmini_20260704.out`.

---
*(S45 original body below — the `7/89` claim is SUPERSEDED by the `1/13` correction above; kept for the
dominant/compressed dispatch framing, which stands.)*
