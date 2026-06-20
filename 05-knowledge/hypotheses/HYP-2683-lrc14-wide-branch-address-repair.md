---
id: HYP-2683
title: LRC(14) wide-branch address repair via private-sector and compatibility profiles
status: OPEN; claimed for exact scout
source: codex-2026-06-20-S55
depends_on:
  - HYP-2675
  - HYP-2682
  - HYP-2681
  - HYP-2677
  - HYP-2650
related:
  - HYP-2240
  - HYP-2241
  - HYP-2530
  - HYP-2455
  - OPEN-Q-108
---

# HYP-2683 - Wide-Branch Address Repair

## Claim Being Tested

Recent work isolates HYP-2675 as the live LRC(14) sector-route crux:

```text
wide row -> p0(E) <= cap_k
```

via Weyl/decorrelation plus the exact plateau bound `Q(k-1)`.  HYP-2682 says
low-height AP resonances need finite phase/support routing before any scalar
estimate is applied.

HYP-2683 imports an older repo lesson: scalar/product shadows mix proof states
until an address coordinate is restored.  The candidate address here is a
private-sector / compatibility profile of the sector wall word:

```text
row E
  -> exact wall atoms
  -> which inner sectors are privately owned by which runners
  -> missed-state compatibility profile
  -> wide/direct-p0 risk routing
```

This is inspired by the owner-private repair in HYP-2240/HYP-2241 and the
compatibility-wall lesson in HYP-2530, but the LRC predicate preserved here is
direct `p0(E) <= cap_k`.

## Planned Exact Scout

Create `04-computation/lrc14_wide_address_repair_codex_s55.py` with output
`05-knowledge/results/lrc14_wide_address_repair_codex_s55.out`.

The scout should:

1. scan exact k=9 true-wide rows in a bounded box already touched by HYP-2675;
2. compute direct `p0`, margin, additive/Freiman diagnostics, state-word
   support, and private-sector owner mass;
3. audit repair channels by mixed-bucket width: scalar shadows, additive
   profiles, private-sector profiles, and full state-address profiles;
4. include Tournament Analysis whose vertices are proof channels, not runners;
5. report whether private-sector or compatibility addresses give a plausible
   finite-router for the HYP-2675 decorrelation/glue proof.

No LRC(14) proof is claimed.
