---
id: HYP-2623
status: CORRECTED (MISTAKE-079) — a=3,4 dips CONFIRMED (infinite families); the "a_max unbounded"
        claim is RETRACTED: the natural family collapses to the floor at a=5. a>=5 reachability OPEN.
source: kind-pasteur-2026-06-19-S9
related:
  - THM-539
  - HYP-2052
  - HYP-2084
  - HYP-2621
  - HYP-2622
---

# HYP-2623: how deep can the LRC max-min second point dip below the floor `1/(k+1)`?

## What is CONFIRMED (exact)

The second spectrum point lives on the Stern-Brocot mediants `a/(a(k+1)-1)`, and `g(k)·k^2 -> 1/a`.
The realized levels:

| level a | family / locus | k | g·k^2 | status |
|--------|----------------|---|-------|--------|
| 2 | doubled apex `{1,..,k-1,2k}` | all k | →1/2 | exact, universal |
| 3 | `F(k,3)={1,..,k-2,k,3(k-1)}` | `k≡7,13,19,25 (mod30)` | →1/3 | exact, infinite family |
| 4 | `F(k,4)` | `k≡1 (mod30)` (k=31,61,…,181) | →1/4 | exact, infinite family (codex-verified) |

So `a_max(k) >= 4` for infinitely many `k`, and the realized `g·k^2 ∈ [≈1/4, 1/2]`.

## What was RETRACTED (MISTAKE-079)

I first claimed `a_max(k)` is UNBOUNDED, with `a=5,6` at `k=211,2311` (`k-1` primorial), hence
`liminf g·k^2 = 0`. **This was wrong.** A covering test returned `M(F(211,5)) < 5/1059`, which I
misread as a *deeper dip*; in fact `M(F(211,5)) = 1/212 = the FLOOR exactly` — the family
COLLAPSES to a tight configuration (`g=0`), it does not dip to level 5. Verified by two
independent covering implementations. `F(k,5)` is TIGHT for every `k` with `2·3·5·7 | (k-1)`
(k=211,421,631,841) and above all small mediants otherwise. So the natural family TOPS OUT at
level 4.

## Adversarial-workflow confirmation (8 agents)

All four verify-claims CONFIRMED by independent exact code (denominator lemma; `M(F(31,4))=4/127`
via crossing + covering; the primorial divisibility pattern; exhaustive small-k `sigma_2`).
E2 independently computed `M(F(211,5))=1/212` (the floor) — the source of the correction.
**E3 exhaustively confirmed `sigma_2(13)=3/41`** (177 below-mediant 13-sets all collapse to the
single witness `{1,…,11,13,36}`), and `sigma_2(11)=2/23`, `sigma_2(12)=2/25` (no dip). **E1's
"beat the family" search concluded: levels top out at 4; generic 2; 3,4 need `k-1` divisible by a
primorial; level 5 collapses to the floor** — i.e. `g(k) = Theta(1/k^2)` with constant `≈ 1/4`–`1/2`
across everything tested. No construction reached `a>=5`.

## The OPEN question (sharpened)

**Is level `a>=5` realizable by ANY gcd-1 k-set, for some k?** Equivalently: is `g(k)·k^2`
bounded below by a positive constant (so `g = Theta(1/k^2)` uniformly), or can it be made
arbitrarily small?
- Evidence FOR bounded (codex S16/S17 + this correction): the F-family saturates at a=4;
  every prefix-containing set re-optimizes to level 1 (`t≈1/k`); no construction reaches a=5.
- The denominator lemma (THM-539 Lemma A, PROVED) says any a>=5 witness needs `max(S) >= q/2 ~
  5(k+1)/2`, i.e. genuinely larger speeds than the a<=4 families use — but that is necessary,
  not sufficient.

## Bonus: a new tight-locus member

`F(k,5)={1,…,k-2,k,5(k-1)}` collapsing to `M=floor` when `2·3·5·7|(k-1)` means it is a NON-AP
TIGHT configuration (max-min `=1/(k+1)`), binding pair `(1, k)` summing to `k+1`. This is a fresh
input to the lonely-measure tight-locus thread (tight locus ⊇ {AP, Goddyn–Wong, …}); worth
checking whether it is AP-equivalent under dilation or genuinely new.

See THM-539, MISTAKE-079, `07-reflections/lrc-spectral-gap-dips-along-primorials-kps.md`,
`04-computation/lrc_a5_collapse_check_kps.py`.
