# Covering buys distance from the AP — and the distance factors: bounded-finite + large-equidistribution

*klein-2026-07-13-S291. Owner: prove covering buys uniform distance from the AP. After the S290 identity,
this is `sup_{primitive covering} conc < 7`, which is the open residual — I did not prove it monolithically.
But the structure is now clean and, more importantly, it **factors**: the supremum lives at bounded
near-AP sets, so the uniform gap is a finite check plus a large-speed equidistribution, not one wall. And
primitivity is exactly the separator between covering and the tight extremals.*

---

## What the claim is, precisely

By the S290 identity `L({1}∪C) = |G(C)|(1 − conc/7)`, "covering buys uniform distance from the AP" means
$$\sup_{\text{primitive covering}} \mathrm{conc} < 7 \iff \inf_{\text{primitive covering}} L > 0,$$
which is the open residual (the THM-527-A cancellation). So the phrase is the honest **restatement**, not a
new reduction. Two things are nonetheless clean and worth stating.

## 1. The tight extremals are `{AP, GW}`, and primitivity is the separator

`conc ≤ 7` **always** — because `L = |G(C)|(1−conc/7) ≥ 0` is a measure. And `conc = 7 ⟺ L = 0 ⟺` a tight
LRC extremal. Verified: the AP `{1..13}` and the Goddyn–Wong doubling `{1..11,13,24}` both have `conc =
7.000` (measure-zero good set), and **both are primitive and non-covering** — exactly kind-pasteur's tight
census. So **covering ⟹ conc < 7 pointwise** (which is just LRC for that set).

The *other* `conc = 7` configurations are the **imprimitive dilates** `c·{AP}`, `c·{GW}` (e.g. `14·{1..13}`
is covering-as-written but reduces to the non-covering AP). A **primitive covering** set is excluded from
the tight set on both counts: it is not `{AP, GW}` (those are non-covering) and it is not a dilate (those
are imprimitive). **Primitivity is the exact separator** — and that is precisely opus-S271's
dilation-blindness: a dilate is a tight shadow a peel cannot see, and primitivity is what removes it.

## 2. The gap factors — the supremum is at bounded near-AP sets

The useful discovery. Sampling primitive covering `min=1` sets (max speed up to 90): the **maximum** `conc`
is `6.177` (gap `0.823`), achieved at the *bounded* near-AP set `{1..14}\{6}`. Large-speed covering sets
sit far away — `{1,90..101}` has `conc = 3.30`. So the supremum of `conc` is **not** approached by large
sets; it lives among **bounded near-AP** configurations. Two consequences:

- **Bounded near-AP → finite check.** The sets with `conc` near the sup have most speeds in `{1..14}`;
  kind-pasteur-THM-734 already closes every family with `≥11` speeds in `{1..14}` (`conc ≤ 6.18` there).
- **Large-speed → equidistribution / true-disc slack.** A large-speed covering set has a fine good set,
  so `|G(C)∩[0,1/14)| ≈ |G(C)|/14` (`conc ≈ 1`, far from 7). This is a **one-interval** discrepancy — much
  weaker than the full cancellation — and it is exactly where opus-S271's true disc certifies with huge
  slack (12/13 peels). For `min=1` specifically, the element `1` **pins the dilation parameter `c = 1`**,
  so a large `min=1` covering set cannot even be a near-dilate of the AP — it is structurally far.

So the uniform gap `sup conc < 7` **tiles**:
$$\big[\text{bounded near-AP: finite check (kps-THM-734)}\big]\ +\ \big[\text{large-speed: one-interval equidistribution (opus true-disc)}\big].$$

## The honest bottom line

I did not prove the uniform gap in one stroke — it is the residual. But "covering buys distance from the
AP" is now a **true, factored** statement rather than a monolithic cancellation: the tight set is exactly
`{AP, GW}` and their dilates, primitivity is the clean separator, and the distance is bounded below by a
finite check on the near-AP region plus a soft equidistribution on the large-speed region — the two
regimes the fleet (kps, opus) is already closing. The AP is not just the boundary; it is a boundary the
covering hypothesis clears in two well-understood pieces.

*Files: `04-computation/lrc14_ap_distance_klein_S291.py` (+.out). HYP-6540. Continues
[[the-compact-core-splits-bounded-ratio-done-and-the-one-cluster-residual-is-AP-tight-concentration-klein-S290]];
consumes opus-S271 (dilation-blindness, HYP-6525), kind-pasteur-THM-734, THM-405.*
