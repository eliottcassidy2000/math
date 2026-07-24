---
source: klein-2026-07-24-S424
status: NEW PROVED NECESSARY CONDITION for tightness at n=13, derived from a "forced witness" reframe.
  Cheap, structural, prunes ~96% of configurations. Also records what it does NOT do (far from sufficient)
  and a failed second-order refinement.
tags: [lrc14, open-q-108, tight-locus, forced-witnesses, unit-residues, necessary-condition, proved]
---

# Forced witnesses and the unit-residue condition

**klein-2026-07-24-S424.** Produced by a deliberately wild reframe: *stop classifying configurations, classify
the witnesses instead.*

## The reframe: the witnesses are forced
For `τ = a/14` with `gcd(a,14)=1`, and any speed `v`,
`‖v·a/14‖ = min(r, 14−r)/14` where `r = va mod 14`, which is `≥ 1/14` **iff** `14 ∤ va` **iff** `14 ∤ v`.
So:

> **Every configuration with no speed divisible by 14 has all six unit fractions
> `{a/14 : a ∈ (ℤ/14)* } = {1/14, 3/14, 5/14, 9/14, 11/14, 13/14}` as witnesses.**

Measured: the tight configurations `AP = {1..13}`, `GW = {1..11,13,24}` and `3·AP` have **exactly** these six
witnesses and no others (and `2·AP` has the dilated twelve). Loose configurations have extra witnesses forming
intervals. Hence the clean restatement:

> **TIGHT ⟺ (no multiple of 14) AND (no witnesses beyond the six forced ones).**

The first clause is trivial; all the content is in "no extra witnesses."

## The new necessary condition (PROVED)
Fix a forced witness `τ = a/14`. The speeds achieving equality `‖vτ‖ = 1/14` are exactly those with
`va ≡ ±1 (mod 14)`, i.e. `v ≡ ±a^{-1}`. Perturb `τ = a/14 + ε`:
- if `va ≡ +1`: `‖vτ‖ = 1/14 + vε` — **increases**;
- if `va ≡ −1`: `‖vτ‖ = 1/14 − vε` — **decreases** below `1/14`.

So moving right is blocked only by a speed `≡ −a^{-1}`, and moving left only by a speed `≡ +a^{-1}`. If either is
absent, the lonely set contains an interval at `a/14`, so `gap > 1/14` and the configuration is **loose**. As `a`
ranges over the units, `±a^{-1}` ranges over all units (the units are closed under negation:
`−1≡13, −3≡11, −5≡9`). Therefore:

> **THEOREM (necessary).** If `S` is tight (`gap = 1/14`) then `S` contains a speed `≡ u (mod 14)` for **every**
> unit `u ∈ (ℤ/14)* = {1,3,5,9,11,13}`.

**Verification.** `AP`, `GW`, `3·AP` all satisfy it. Of 2865 sampled configurations missing at least one unit
residue, **zero** have `gap ≤ 1/14`; the minimum observed is `gap·14 = 1.663`, comfortably loose.

**Power.** Only `4.4%` of random 13-subsets contain all six unit residues — the condition **prunes ~96% of
configurations outright**, structurally and at zero cost. It is an ideal first filter for any search or
enumeration (and for the defect ladder's finite checks).

## What it does not do (honest limits)
- **Far from sufficient.** Among 400 sampled configurations that *do* satisfy it, the gap distribution is
  `min 1.878`, `median 2.601`, `max 4.667` (units of `1/14`); **none** is within 10% of tight. The residual —
  "no lonely component *away* from the unit fractions" — is the genuinely hard content, and this condition says
  nothing about it.
- **The second-order refinement fails.** I hoped that *how many* speeds sit on unit residues would refine the
  condition. It does not: `corr(#unit-residue speeds, gap·14) = +0.249` on small samples — weak, and in the
  *opposite* direction to the naive guess. Recorded so it is not retried.

## Why 14 is the special modulus here
The argument needs the witness to be *forced*, which happens only because `‖v·a/14‖ ≥ 1/14` automatically for
`14 ∤ v`. For any other modulus `q` the analogous points are not forced, so no comparable condition follows.
This is the same "modulus-14 rigidity" opus noted, now with the precise mechanism and a usable consequence.

→ klein-S419 (defect-1 classification), klein-S422/S423 (peeling lemma, clustering frame), opus-S4 (band-width,
modulus-14 rigidity remark), OPEN-Q-108.
