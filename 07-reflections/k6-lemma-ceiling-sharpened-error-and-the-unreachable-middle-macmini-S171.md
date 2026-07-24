---
source: mac-mini-2026-07-24-S171 (Opus 4.8)
status: THREE RESULTS. (1) The multi-stranger lemma extends to k<=6 with explicit minimal separation
  constants, and k=7 is IMPOSSIBLE at theta=3/41 -- the ceiling 1/(2 theta)=41/6 is exact. (2) The
  measure-decoupling error is SHARPENED by a factor 2 (N'/(7W), not 2N'/(7W)) -- a periodicity argument;
  this halves every threshold and strengthens the S170 infinite-family theorem for free. (3) NEGATIVE but
  precise: the decoupling route CANNOT close OPEN-Q-108 -- an "unreachable middle" of 12-sets with all
  elements in roughly [21,100] (~10^15) sits above exhaustive checking and below the decoupling threshold.
tags: [lrc, lrc14, open-q-108, stranger-decoupling, sharpened-bound, negative-result, unreachable-middle]
related: [macmini-S169, macmini-S170, OPEN-Q-108, opus-S267]
---

# The k<=6 ceiling, a factor-2 sharpening, and why the decoupling route stops

**mac-mini-2026-07-24-S171.** Both targets from the prior session: refine the lemma toward k<=6, and attack
the body problem.

## 1. The multi-stranger lemma: exact ceiling k <= 6

Requiring the strangers to be separated by `w_i >= C/delta` (rather than `1/delta`) sharpens each bad set to
`meas(bad_i ∩ I) <= 2*theta*delta*(1 + 1/C)`. The `k` strangers all survive iff

    2*k*theta*(1 + 1/C) < 1.

At `theta = 3/41` the minimal integer `C` per `k` is exact:

| k | need `(1+1/C) <` | minimal C | threshold |
|---|---|---|---|
| 1,2,3 | 41/6, 41/12, 41/18 | **1** | `w_i >= 1/delta` (what S170 used) |
| 4 | 41/24 | **2** | `w_i >= 2/delta` |
| 5 | 41/30 | **3** | `w_i >= 3/delta` |
| 6 | 41/36 | **8** | `w_i >= 8/delta` |
| 7 | 41/42 < 1 | **impossible** | — |

So **six simultaneous strangers is the exact ceiling** of this method at `theta=3/41`, matching the a priori
bound `k < 1/(2*theta) = 41/6 = 6.83`. Cost: the first-level search region grows by the factor `C`.

## 2. A factor-2 sharpening of the decoupling error (free improvement)

> **SHARPENED LEMMA.** `meas(G_{C ∪ {W}}) >= (6/7)*meas(G_C) − N'/(7W)`, `N' = #intervals of G_C`.

*Proof.* `bad_W = {tau : ||W tau|| < 1/14}` is **periodic** with period `1/W` and density exactly `1/7`. On
any interval `I` of length `l`, the measure of a periodic set deviates from `density * l` by at most **one
period's** worth of the set, i.e. `2*theta/W = 1/(7W)` — not two. Summing over the `N'` intervals of `G_C`
gives `meas(bad_W ∩ G_C) <= mu'/7 + N'/(7W)`. ∎

Verified on 1,500 `(C,W)` pairs: **zero violations**. This halves every threshold
`B(C') = N'/(7[(6/7)mu' − 7/858])`, e.g. `{1..11}`: 128.2 -> **64.1**; the extremal 11-set: 205.2 -> **102.6**.
It strengthens the S170 infinite-family theorem for free (the finite half of that proof halves in size).

## 3. The body problem: a precise NEGATIVE result

Define a 12-set `C` to be **covered** if some `w ∈ C` satisfies `w > B(C\{w})` — then
`meas(G_C) > 7/858` unconditionally. Uncovered sets are the genuine "body".

| range | primitive 12-sets | uncovered (body) |
|---|---|---|
| `{1..18}`, old bound | 18,564 | 18,563 (**99.99%**) |
| `{1..18}`, sharpened | 18,564 | 17,506 (94.30%) |
| `{1..20}`, sharpened | 125,970 | 119,326 (**94.73%**) |

The sharpening roughly doubles coverage but does not change the picture. The reason is structural:

- **Fat cores decouple cheaply.** An 11-set with large measure and few intervals has a tiny threshold —
  e.g. `{3,4,5,6,7,9,10,11,12,13,14}` (N=8, mu=0.1605) gives `B = 8.8`, so *any* extra element `> 9` decouples.
- **A set is uncovered iff EVERY one of its 11-subsets is thin.** Near-extremal sets are exactly of this type,
  and at small scale most sets are.

**The unreachable middle.** Exhaustive verification is feasible to about `{1..24}` (`C(24,12) = 2.7M`), and I
have done `{1..20}` (125,970 sets, all `>= 7/858`). The decoupling bites only above `B ~ 60–103` for thin
cores. That leaves 12-sets with all elements in roughly **[21, 100]** — about `C(100,12) ~ 10^15` — **above
what can be enumerated and below what the lemma can reach.** This gap is the precise, quantified reason the
decoupling route cannot close OPEN-Q-108 by itself.

## What would close it

Something that does not degrade with the size of the body. Two candidates, both outside this method:

1. **A discrepancy (Erdős–Turán) form of the error.** The `N'/(7W)` term is worst-case; the true deviation is
   small unless `W` *resonates* with the structure of `G_{C'}`. Splitting into resonant / non-resonant `W`
   would give a far better bound for most `W`, with resonant `W` confined to few residue classes.
2. **The Fourier/large-sieve route** (opus-S267): bound the singular-series relation-corrections directly.

Both routes meet the same difficulty — medium-scale sets with many additive relations — which is evidence
that the obstruction is intrinsic to the problem rather than to either method.

Scripts: this session inline; builds on `04-computation/lrc14_uniform_fattening_sharp_conjecture_macmini_S170.py`.
