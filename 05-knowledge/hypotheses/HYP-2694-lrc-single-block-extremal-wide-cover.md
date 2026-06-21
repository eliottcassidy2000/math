---
id: HYP-2694
title: The LRC(14) wide-residual extremal shape is the single coherent block (the apex-prime partition-function twin of max-c3=regular)
status: OPEN; coherent-block quotient PROVED by THM-557, full arbitrary-shape/multi-carrier discrepancy still open
source: kind-pasteur-2026-06-20-S21
depends_on:
  - HYP-2675   # wide => p0 <= cap (the sole LRC sector residual)
  - THM-554    # the score partition function (the tournament twin)
  - THM-555    # cut-space/cycle-space wall; max-c3=regular extremal shape
  - THM-557    # coherent-block quotient and single-block diagonal-freeze error
related:
  - HYP-2644   # far-element plateau
  - HYP-2653d  # corrected wide-bound framing
  - THM-027    # max-c3 = regular (the tournament extremal shape)
  - HYP-2606   # the relation-lattice signed sum (the cycle-space correction)
---

# HYP-2694 - The single block is the LRC wide-cover extremizer

## The claim

LRC(14) is `p0(E)=meas(S7(E)) <= cap_k` for primitive k-sets, k=8..12 (k<=7 pigeonhole; finite
check span<=14 DONE). For WIDE E (span>14), `p0` factors into a DECORRELATED part (runners as
independent particles on Z/7 = the cut-space inclusion-exclusion) plus the decorrelation error (the
cycle-space relation-lattice correction). The decorrelated cover is a partition function over the
CLUSTER SHAPE of E. CLAIM:

> **The sup over wide cluster shapes of the decorrelated cover is the SINGLE COHERENT BLOCK**
> `{M, M+1, ..., M+m-1}` (m=k-1, one all-sweeping cluster), and its decorrelated cover is
> **< cap_k with margin >= 0.19** for k=8..12. Splitting into >=2 clusters strictly LOWERS the cover.

This is the apex-prime partition-function TWIN of THM-027/555 "max c3 = regular score": maximize the
cut-space invariant and the gas condenses to its most-coherent occupancy (single block) / most-balanced
occupancy (regular score).

## Evidence (decorrelated cover = mean over (anchor phi, slow x) of the coverage indicator)

**S61 exact update / THM-557.**  The S21 floating-grid evidence has now been
replaced by exact shared-`x` Fraction arithmetic in
`04-computation/lrc14_single_block_extremality_margin_codex_s61.py`.  In the
HYP-2694 coherent-block quotient, the anchored speed `0` is fixed as its own
cluster and the `m=k-1` nonzero speeds are partitioned into far coherent
consecutive blocks.  Exact enumeration of all integer partitions of `m=7..11`
proves that the one-part block `[m]` is the decorrelated maximizer.

Single-block decorrelated cover `D_m=p0_decorr([m])`,
`E={0}U{M..M+m-1}`, vs `cap_k`:
```text
k= 8 (m=7):  D_m=283/1470     cap-D_m=1111/5880
k= 9 (m=8):  D_m=629/2058     cap-D_m=111019/588588
k=10 (m=9):  D_m=16969/41160  cap-D_m=102803/535080
k=11 (m=10): D_m=30551/61740  cap-D_m=184957/802620
k=12 (m=11): D_m=71111/123480 cap-D_m=34729/123480
```

The closest split is always `[m-1,1]`, and its exact loss from the single block
is:

```text
m=7:  1111/10290
m=8:  374/5145
m=9:  6561/96040
m=10: 42661/864360
m=11: 9047/172872
```

THM-557 also gives an elementary diagonal-freeze error bound for actual shifted
single blocks:

```text
E_M={0} union {M,...,M+m-1},
|p0(E_M)-D_m| <= 7*binom(m,2)/M.
```

Thus the single-block branch is rigorously below cap once `M` exceeds
`779,1040,1312,1367,1369` for `k=8..12`, respectively.  These cutoffs are
conservative; exact samples at `M=19` already have error magnitude at most
`0.0125` and remain below cap.

## Why this is the right object (the abstraction)

The single-block cover = the phi-shift-average of the consec_m sector pattern = the "most coherent
occupancy" on Z/7. It is the cut-space (single-particle) extremum, exactly as the regular score is the
tournament cut-space extremum. The decorrelated cover is a PARTITION FUNCTION over the cluster shape;
its sup is the most-coherent shape; and it sits comfortably below cap. See reflection
`the-apex-prime-partition-function-tournaments-and-runners-are-one-gas`.

## What remains (the cycle-space residual, now with a 0.19 budget)

1. **Single-block is the coherent-block decorrelated sup** — PROVED by
   THM-557/S61 in the exact quotient where the nonzero runners are partitioned
   into coherent consecutive blocks.  The remaining rearrangement question is
   whether arbitrary bounded cluster shapes can always be compressed to this
   coherent-block quotient without decreasing the decorrelated cover.
2. **Single-block decorrelated cover < cap_k** — PROVED exactly by THM-557 for
   `k=8..12`, with minimum cap margin `111019/588588`.
3. **Single-block finite-scale error <= margin** — PROVED by THM-557 once
   `M > 7*binom(m,2)/(cap_k-D_m)`.  The finite small-`M` part is a bounded
   exact check and is numerically much safer than the conservative cutoff.
4. **Multi-block decorrelation error <= margin** — the genuine analytic
   content still open: prove a joint carrier BV/Erdos-Turan bound for several
   separated blocks, using the exact split gap plus the cap margin as budget.

Closing the remaining arbitrary-shape compression and joint multi-block
decorrelation-error bounds would close HYP-2675, hence the LRC(14) sector
route.
