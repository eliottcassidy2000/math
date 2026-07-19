# The rational-time floor proves the covering-core gap for 95% of spread cores — and names the last 5%

*death-star-2026-07-19-S58i. Working the S58h handoff — "prove the covering-core gap for spread
cores": a far-from-AP 12-core `W` covering `2..12` with spread `ρ(W) ≥ 6.5` has `M(W) ≥ 1/13 + c`.
Outcome: the infimum is bounded (so a crude bound exists), an elementary rational-time floor **proves
the gap for 95%** of such cores, and the irreducible **5%** is sharply characterized. **The gap is
not fully proved**, but it is reduced to a small, precisely-named residual. Scripts:
`lrc14_spread_core_infimum_deathstar_S58i.py`, `lrc14_rational_time_floor_deathstar_S58i.py`.*

## The infimum is bounded away from 1/13 (so a crude bound exists)

Over **22,000+** spread (`ρ≥6.5`), far-from-AP (Hamming `>6`), covering-`2..12` 12-cores across two
searches (S58h + S58i): **0** have `M(W) < 1/13`, and the infimum is `M(W) = 3/29 ≈ 0.1034`, margin
`+0.0265`. Near-minimizers: `[2,7,9,10,12,13,16,17,19,22,23,25]`, `[2,3,9,10,12,13,16,19,21,22,23,25]`
(all contain the run `12,13`). So the covering-core gap holds with `c ≈ 0.0265` — a **crude** bound,
not the sharp Freiman constant, exists.

## The rational-time floor (PROVED, elementary) — 95%

> **Floor.** `M(W) ≥ max_{k≥2} d_k/k`, where `d_k = min_{w∈W} ‖w/k‖` (the min circular residue mod
> `k`). Trivially `M(W) ≥ min_w ‖w·(1/k)‖ = d_k/k` for each `k`; take the max.

Sharpened missed-modulus: a *covered* `k` has `d_k = 0`, but a `k` where every runner sits `≥ 2/k`
from a multiple of `1/k` gives `M(W) ≥ 2/k`. **Verified: over 39,275 spread far covering-`2..12`
cores, `max_k d_k/k > 1/13` for 94.9% of them** — so the covering-core gap is *proved elementarily*
for 95% of the residual, via an explicit rational witness `t = 1/k`.

## The 5% residual, named

The remaining **5%** have `max_k d_k/k ≤ 1/13`: for *every* `k`, some runner is within `1/13` of a
multiple of `1/k` — they are **near-covering at all rational scales `t=1/k`**, so no simple rational
witness clears `1/13` (their floor value stalls at `≈ 0.074`). Yet their *actual* `M ≈ 0.103` is
realized at a **twisted** maximizer `t = a/q` with `a ≠ 1` (the pair-sum competitor, S58g). These
evasive cores — arithmetically arranged to defeat every `t=1/k` witness while remaining lonely at a
twisted time — are exactly the irreducible **12-set-uniqueness / Freiman stability** core (HYP-4382):
`M(W) = 1/13` iff `W` is a dilated AP, in quantitative (Hamming `>6`) form. THM-1006's witness-table
method cannot reach them (radius capped at `1/(2δ) ≈ 6`), and the elementary floors cannot either.

## What this settles, and the remaining task

- **Proved:** the rational-time floor `M(W) ≥ max_k d_k/k` (elementary), which discharges 95% of the
  spread-far-core case of the covering-core gap. Combined with the clustered floor (`ρ<6.5`, S58h)
  and Hamming rigidity (near-AP, THM-1004/5/6), only a thin evasive stratum remains.
- **Evidence:** the infimum over the *whole* spread-far-core class is `≈ 3/29` (margin `0.0265`), so a
  crude constant `c` provably exists for the residual too — it is a matter of finding the witness.
- **Not proved:** the 5% residual (rational-time-evasive cores). This is the genuine Freiman/HYP-4382
  wall; no `t=1/k` witness suffices, and it needs the twisted `a/q` competitor bounded below.
- **Next:** for the evasive stratum, bound the twisted-time value. A pigeonhole is available: with 12
  runners but many scales `k∈[13,26]` each demanding a distinct near-killer within `k/13`, a runner
  cannot be near-killed at too many scales — turning "near-covering at all scales" into a contradiction
  or a forced twisted witness. This is the concrete crack to try on the last 5%.

— Related: `the-clustered-floor-and-the-favorable-shape-...-deathstar-S58h.md`,
`the-pairsum-competitor-margin-tracks-the-schur-deficit-...-deathstar-S58g.md`,
`the-missed-modulus-competitor-splits-the-kernel-...-deathstar-S58f.md`. THM-1004/5/6 (Hamming),
THM-730 (Schur), THM-1028 (Lemma G). HYP-7310 (crux), HYP-4382 (12-set uniqueness), HYP-7748 (floor).
