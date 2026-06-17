---
id: THM-522
title: Two exact levers for the LRC(14) lonely measure — scale-invariance L(cS)=L(S), and quantization L(S)∈(1/(14·lcm S))·ℤ (so L>0 ⟹ L≥1/(14·lcm S)) — reducing inf L>0 to a compactness statement (the extremizers have bounded lcm); and a correction of the inf estimate (≤1/1260, via the minimal AP-perturbation 12→36, not the multiple-of-14 family's 0.0052)
status: PROVED (scale-invariance: τ↦cτ measure-preserving; quantization: rational arc endpoints — both elementary). VERIFIED exactly (`lrc14_exact_lonely_measure_kps.py`, `lrc14_tight_locus_and_true_inf_kps.py`): L(cS)=L(S) for c=1,2,3,5,6; L·14·lcm∈ℤ on all tested configs; the minimal single-element perturbation of the tight AP {1..13} is 12→36 giving L=1/1260≈0.000794 (vs the repo's restricted-search 0.0052), confirmed by an independent fine grid; 2-element perturbations (w≤72) found nothing below 1/1260.
source: kind-pasteur-2026-06-16-S7
depends_on:
  - THM-501   # the lonely measure L(S)=lim D(q,S)/q = meas{τ:||vτ||>1/14 ∀v}
  - THM-515   # L(S) = the lonely measure (the form these levers act on)
  - THM-518   # stranger-decoupling (the m→∞ escape direction; complementary to the lcm-bounded reduction)
related:
  - HYP-2561  # inf L = 1/1260 conjecture; the tight locus is finite
  - OPEN-Q-097 # inf L>0 (now reframed via the tight locus + quantization)
  - MISTAKE-075 # the repo's inf≈0.0052 was a restricted-search artifact (multiple-of-14 strangers)
  - reflection: the-lrc-inf-is-a-quantized-gap-from-the-tight-locus-and-the-extremizer-is-a-sporadic-perturbation-kps
external: Lonely Runner Conjecture, n=14; Csorba/Goddyn-Wong tight-config classification.
---

# THM-522 — scale-invariance and quantization of the LRC(14) lonely measure

**Object.** `L(S) = meas{ τ∈[0,1) : ||vτ|| > 1/14 ∀v∈S }` (THM-515; the lonely measure =
the singular series). `inf_S L(S) > 0 ⟹ C'(14) ⟹ LRC(14)` (THM-398/501). The danger set of a
speed `v` is `D_v = {τ: ||vτ||≤1/14} = ∪_k [(14k−1)/(14v), (14k+1)/(14v)]` (`v` arcs, total
measure `1/7`), and `L(S) = 1 − meas(∪_{v∈S} D_v)`.

## A. Scale-invariance: `L(cS) = L(S)` (PROVED)

For any nonzero integer `c`, `L(cS) = ∫_0^1 ∏_{v∈S} 1_{||cvτ||>1/14} dτ`. The map
`τ ↦ cτ mod 1` is a `|c|`-to-1 **measure-preserving** cover of `[0,1)`, so for any 1-periodic
`f`, `∫_0^1 f(cτ)dτ = ∫_0^1 f(u)du`. With `f(u)=∏_{v}1_{||vu||>1/14}`, **`L(cS)=L(S)`**. ∎

Consequences: (i) `L` depends only on the **primitive projective shape** of `S` (WLOG
`gcd(S)=1`); (ii) the "stranger scale" and the overall scale are the *same* degree of freedom —
the THM-518 family `({1..13}\{j})∪{14m}` and its dilations carry identical `L`.

## B. Quantization: `L(S) ∈ (1/(14·lcm S))·ℤ`, so `L>0 ⟹ L ≥ 1/(14·lcm S)` (PROVED)

Every arc endpoint `(14k±1)/(14v)` has denominator dividing `14·lcm(S)` (since `v ∣ lcm S`),
so all breakpoints of `∪_{v∈S} D_v` lie in `(1/(14·lcm S))·ℤ`; the union's measure (a sum of
merged-interval lengths) and hence `L = 1 − meas` lie in `(1/(14·lcm S))·ℤ`. As `L≥0`:

> **`L(S) > 0 ⟹ L(S) ≥ 1/(14·lcm S)`.**

VERIFIED: `L·14·lcm ∈ ℤ` for the AP, the extremizers, the sporadic, and Sidon sets.

## C. The compactness reduction (the reframe of the prize)

For any family `F` of primitive 13-sets with `lcm(S) ≤ M` (`∀S∈F`):
`inf{ L(S) : S∈F, L(S)>0 } ≥ 1/(14M) > 0` (B). Hence:

> **`inf_S L(S) > 0` ⟸ the `L`-minimizing configs have bounded `lcm`** (quantization B gives
> `inf ≥ 1/(14M)` over any `lcm≤M` family). Scale-invariance (A) kills dilation; THM-518
> stranger-decoupling kills the one-entry-`→∞` direction (`L → (6/7)·meas(bounded core) > 0`).

> **CORRECTION (THM-523, the prove/disprove dialectic):** the clean form "small `L` ⟹ bounded
> `lcm`" is **FALSE** — near-tight 12-cores can carry `lcm ~ 10^7` with small `L` (e.g.
> `{1,3,4,5,7,8,9,10,11,12,13,17,25}`, `lcm≈3·10^7`, `L=2/425`). The quantization `L≥1/(14 lcm)`
> still holds, but the correct compactness bounds the **perturbing ELEMENTS** (via the
> decoupling floor `L(C∪{w})≥(6/7)meas(G_C)−r/(7w)`, THM-523 B), not the lcm. The true reduction
> of `inf L>0` is THM-523 D: a **uniform lower bound on `meas(G_C)` over all 12-subsets** (≡
> finiteness of the tight locus), NOT "bounded lcm of the extremizers."

## D. Correction of the inf estimate (MISTAKE-075): `inf L ≤ 1/1260`, not `0.0052`

The tight locus (`L=0`) is **not just the AP** `{1..13}`: it includes sporadics, e.g.
**`{1..11,13,24}`** (`L=0`, the HYP-2055 sporadic, verified). The repo's inf `≈0.0052`
(THM-518) came from restricting extremizers to **multiple-of-14 strangers** `({1..13}\{j})∪{14m}`.
But the **minimal single-element perturbation of the tight AP** is `12→36`:

> **`{1,2,…,11,13,36}` has `L = 1/1260 ≈ 0.000794`** (≈6.7× below the restricted 0.0052),
> verified by exact arc-sweep AND an independent fine grid. (`36` is not a multiple of 14;
> both the AP and the sporadic `{1..11,13,24}` reach this config in one move.)

So **`inf_S L(S) ≤ 1/1260`** over primitive 13-sets, and the true extremizers are minimal
perturbations of the **full tight locus** (AP + sporadics), not the multiple-of-14 family.
Exact search found **nothing below `1/1260`** (1-drops `w<400`, 2-drops `w≤119`; 2-drop floor
`1/980`). And the **tight locus is sparse**: among 1-element replacements of `{1..13}` with
`w≤120`, the *only* tight (`L=0`) config besides the AP is the sporadic `{1..11,13,24}`; zero
2-element tights with `w≤79` — strong evidence the primitive tight locus is finite (HYP-2561).
This is the MISTAKE-073 lesson recurring: sweep the full orbit; the extremizer is never the
obvious family member.

## Scope / honesty

PROVED: A (scale-invariance), B (quantization + the `1/(14 lcm)` lower bound), C (the
compactness reduction — a logical equivalence, with the open content named). VERIFIED exactly:
the `1/1260` minimal AP-perturbation and the sporadic tight `{1..11,13,24}`. The claim
`inf L = 1/1260` (vs merely `≤`) is CONJECTURAL (HYP-2561) — extensive exact search (1-drops
`w<400`, 2-drops `w≤72`) found nothing below, but 3-drops / larger / other tight-locus
perturbations are unsearched, and whether `inf L > 0` at all still hinges on the tight locus
being finite/bounded-lcm (C). The `0.0052→1/1260` correction is a genuine sharpening, not a
disproof of LRC(14) (every config found is loose, `L>0`). Cross-links: THM-501/515/518,
HYP-2055 (the sporadic), HYP-2561, OPEN-Q-097, MISTAKE-075.
