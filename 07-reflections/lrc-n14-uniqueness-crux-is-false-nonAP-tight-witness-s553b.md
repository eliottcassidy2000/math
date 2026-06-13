---
source: opus-2026-06-01-S553 (remote-control, disproof hunt)
status: RIGOROUS DISPROOF of a reduction's premise (NOT a disproof of LRC); two exact methods
tags: [LRC, n14, tight-family, disproof, sporadic, non-transversal, non-unit-pair, S552, S553, HYP-2055]
---

# LRC@14: the "AP is the unique tight 13-set" crux is FALSE

**Prompt (user):** see if a disproof is possible, now that the crux "AP is the
unique tight 13-set" has been identified.

## Result

The crux is **refuted**. There is a non-AP tight 13-set at n=14:

> **`V* = {1,2,3,4,5,6,7,8,9,10,11,13,24}`**  (the AP with `12 → 24`)

verified by **two independent exact computations**:
- max-collar `M(V*) = max_t min_i ||v_i t|| = 1/14` exactly (crossing/peak
  candidate set), attained at `t = 3/14`;
- safe set at threshold `1/14` has **measure 0** with witnesses
  `{1,3,5,9,11,13}/14` (so tight = nonempty + measure zero).

`V*` is gcd-1 and not a dilate of the AP. So **the AP is NOT the unique tight
13-set**, and the clean route "AP-unique-tight ⇒ LRC@14" is dead — exactly as the
sporadic tight sets at n=5,6,8 already warned.

## Why the earlier census missed it, and what it is

The exhaustive census (HYP-2055/C1) ran over speeds `≤ 21 = 1.5n`; `V*` has speed
`24`. So C1 ("AP unique tight among speeds ≤1.5n") stands, but its **lift to all
speeds is false**. `V*` is precisely an inhabitant of the **non-unit-pair hole**
flagged in C3: mod `2n-1 = 27`, its residues **miss the non-unit antipodal pair
`{12,15}`** and **double `{3,24}`** (both 3 | 12,15,3,24). It is therefore a
**non-transversal** that oracle-S553's Link-1 witness `t = a^{-1}/27` cannot
reach (`a` must be a unit). A local search confirms `V*` is the **unique
distance-1-from-AP** non-AP tight set.

## What is NOT disproved

**LRC(14) itself stands.** `V*` is *tight* (M = 1/14), i.e. still lonely, just at
the boundary. The hunt found **no** config with `M < 1/14` across:
- all 8191 gcd-1 antipodal transversals mod 27 (speeds ≤ 26);
- 120 000 non-transversals missing a non-unit pair mod 27 (speeds ≤ ~80);
- 250 minimise-measure hill-climbs (speeds ≤ 60);
- the AP add/remove neighbourhood (added speed ≤ 80).

So a *disproof of LRC@14* was not found; what was disproved is the **uniqueness
reduction**.

## The corrected picture

The right invariant is the spectral gap, not uniqueness. LRC@14 ⟺ **every config
has `M ≥ 1/14`**. The tight boundary `M = 1/14` is **not** just the AP — it now
provably includes at least one sporadic `V*` (and likely a family of non-unit-pair
non-transversals). All known tight configs are lonely at the `j/n` lattice (C2),
so the sporadics neither help nor threaten LRC; the content is entirely "nothing
dips below `1/14`."

**Handoff:**
- to oracle (S552/S553): your "M = 1/n only for the AP-tight family" needs the
  non-transversal exception; `V* = {1..11,13,24}` is an explicit n=14 witness in
  the non-unit-pair hole. The reduced gap must be proved on **transversals ∪
  (non-unit-pair non-transversals)**, not transversals alone.
- characterise the full sporadic tight family at n=14 (search ≥ distance-2 from
  AP, speeds to ~2n); conjecture: all are non-unit-pair non-transversals, all
  lonely at `j/14`.

**Artifacts:** `04-computation/lrc_n14_disproof_hunt_s553.py` (+`.out`),
`05-knowledge/results/lrc_n14_nonAP_tight_family_s553.out`. Updates **HYP-2055**
(status → partially refuted). Builds on oracle-S552/S553, HYP-2052.
