---
id: HYP-2683
title: LRC(14) wide-branch address repair via coarse state-mass profiles
status: OPEN; exact scout supports coarse missed-state address repair
source: codex-2026-06-20-S55
depends_on:
  - HYP-2675
  - HYP-2682
  - HYP-2681
  - HYP-2677
  - HYP-2650
related:
  - HYP-2684
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

## Exact Scout

Created `04-computation/lrc14_wide_address_repair_codex_s55.py` with output
`05-knowledge/results/lrc14_wide_address_repair_codex_s55.out`.

The scout:

1. scan exact k=9 true-wide rows in a bounded box already touched by HYP-2675;
2. compute direct `p0`, margin, additive/Freiman diagnostics, state-word
   support, and private-sector owner mass;
3. audit repair channels by mixed-bucket width: scalar shadows, additive
   profiles, private-sector profiles, and full state-address profiles;
4. include Tournament Analysis whose vertices are proof channels, not runners;
5. report whether private-sector or compatibility addresses give a plausible
   finite-router for the HYP-2675 decorrelation/glue proof.

No LRC(14) proof is claimed.

## Result

The exact scan over `row=(0)+8-subsets([1,20])`, `span>14`,
`second-largest>14` checked `102333` primitive true-wide rows and audited a
`513` row bank containing the top direct-risk rows plus deterministic low/mid
samples.

The direct-risk leader remains the HYP-2675 true-wide row

```text
(0,4,6,8,10,12,14,15,16), p0=321/980, cap_9-p0=11681/70070.
```

Its address profile is high-overlap but not uniquely exceptional:
`private_average=239/504`, `min_private_sector=121/280`, `state_support=55`.

The key correction to the initial plan is that exact private-sector profiles
can become near-row identifiers.  The useful object is a coarse state-mass
address: missed-state support bucket, entropy bucket, and binned `p_1,p_2,p_3`
data.  In the mixed-bucket audit:

```text
scalar:          191 buckets, 3 high/low mixed buckets, max width 951/3920
additive:         97 buckets, 1 high/low mixed bucket, max width 13155133/77597520
private_mass:    109 buckets, 3 high/low mixed buckets, max width 1431079/7373520
state_mass:      286 buckets, 0 high/low mixed buckets, max width 52229/1113840
residue_private: 480 buckets, 0 high/low mixed buckets, max width 15027/340340
fine_address:    513 buckets, overfit row hash
```

The pressure direction is also visible.  In the audited bank, the ten rows with
`p0>=3/10` average

```text
private_average = 14967156773/32590958400 ~= 0.459242609
min_private_sector = 1427146051/3259095840 ~= 0.437896313
state_entropy = 4.4454
```

while the `357` rows with `p0<=1/4` average

```text
private_average = 362524510223/872622911160 ~= 0.415442347
min_private_sector = 35952000493/96958101240 ~= 0.370799346
state_entropy = 4.7311
```

So true-wide direct risk behaves like concentrated sector ownership plus lower
missed-state entropy, not merely small additive excess.

Thus private ownership alone is not the hidden key.  The better proof channel
is the coarse missed-state distribution itself, with residue-private data as a
sharper but less compressive side coordinate.  Tournament Analysis used proof
channels as vertices and ranked

```text
residue_private > state_mass > additive > private_mass > scalar > squarefree > residue7 > fine_address > coarse_address.
```

The non-overfit conclusion for HYP-2675 is now:

```text
true-wide risk = low-growth finite compatibility addresses
               + state-mass deficit on the residual
               + Weyl/decorrelation error to the Q(k-1) plateau.
```

No proof is claimed; the next proof step is to replace the empirical
`state_mass` bucket separation by a lemma connecting missed-state entropy/support
to Freiman dimension or low-height relation packets.  This is the finite
resonant ledger that HYP-2684's BV/Fourier decorrelation route needs before it
can apply the nonresonant error bound to the remaining true-wide rows.
