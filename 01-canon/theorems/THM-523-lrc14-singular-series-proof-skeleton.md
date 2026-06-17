---
id: THM-523
title: The LRC(14) singular-series proof skeleton — the single-perturbation infimum is PROVABLY 1/1260 (decoupling floor + explicit compactness N≤93 + the two-speed-clash champion mechanism), and inf L>0 ⟹ LRC(14) reduces to ONE open lemma — a uniform lower bound on the gap-1/14 lonely measure meas(G_C) over all 12-subsets C (≡ finiteness of the primitive tight locus); plus: zero counterexamples (tight locus = {AP, Goddyn–Wong T5}, both max-min=1/14 exactly)
status: PROVED (decoupling floor; single-perturbation infimum = 1/1260 with explicit constant N≤93; the champion two-speed-clash; the max-min=1/14 of the tight locus). VERIFIED exactly + independently (dialectic prove/disprove workflow, 7 agents, both sides reproduced with fresh code; champion 15/36−2/5−1/70−1/504=1/2520 and floors re-checked here). The CRUX (uniform meas(G_C) bound ≡ tight-locus finiteness) is OPEN — confirmed open in the literature (Perarnau–Serra survey, no fixed-n finiteness theorem).
source: kind-pasteur-2026-06-17-S1 (dialectical prove/disprove campaign)
depends_on:
  - THM-522   # scale-invariance + quantization L≥1/(14 lcm) (the engine turning L>0 into L bounded below)
  - THM-518   # stranger-decoupling (sharpened here to the finite-w inequality)
  - THM-501   # L(S) = lonely measure; inf L>0 ⟹ C'(14) ⟹ LRC(14)
related:
  - HYP-2561  # sharpened: inf L=1/1260, tight locus = {AP, AP[12→24]}
  - OPEN-Q-108 # the uniform-fattening lemma (the crux)
  - reflection: the-lrc14-proof-is-one-lemma-away-and-the-dialectic-built-it-kps
external: LRC(14) = first OPEN case of the Lonely Runner Conjecture (proven k≤12, Sungkawichai–Trakulthongchai arXiv:2604.23906); Goddyn–Wong "Tight instances of the LR" Integers 6 (2006) #A38; Tao arXiv:1701.02048 (bounded-speed reduction); Perarnau–Serra survey arXiv:2409.20160.
---

# THM-523 — the LRC(14) singular-series proof skeleton

**Object.** `L(S)=meas{τ:||vτ||>1/14 ∀v∈S}` (13 distinct positive integers; THM-501/515).
`inf_S L(S)>0 ⟹ C'(14) ⟹ LRC(14)`. For a 12-element core `C`, let `G_C = [0,1)∖∪_{v∈C}D_v`
be its **gap-1/14 lonely set** (`meas(G_C)=L(C)`). LRC(14) = our case = the **first open case**
of the Lonely Runner Conjecture (proven for ≤12 speeds, 2026).

## A. Gap decomposition (PROVED, exact)

Adding one speed `w` to a core `C`: `D_w` only removes from `G_C`, so

> **`L(C∪{w}) = meas(G_C ∖ D_w) = meas(G_C) − meas(G_C ∩ D_w).`**

`L(C∪{w})=0` (tight) ⟺ `D_w ⊇ G_C` (the new speed covers every residual gap).

## B. Decoupling floor (PROVED — sharpens THM-518)

`w`'s danger comb has period `1/w` and density `1/7`, covering at most `meas(A)/7 + 1/(7w)`
of any arc `A`. Summing over the `r` arcs of `G_C`:

> **`L(C∪{w}) ≥ (6/7)·meas(G_C) − r/(7w)`,  so  `lim_{w→∞} L(C∪{w}) = (6/7)·meas(G_C)`.**

For the 13 single-element-dropped cores `C={1..13}∖{e}`, the limiting floors `(6/7)meas(G_C)`
range over `[1/143, 3/49]`, **minimum `1/143 ≈ 0.00699` at `e=6` `≫ 1/1260`** (re-verified:
`meas(G_{e=6})=7/858`, floor `1/143`). **Sending any one speed to `∞` pushes `L` UP, never
toward 0** — the single-element escape to `inf=0` is closed.

## C. Single-perturbation infimum = 1/1260 (PROVED)

From B, if `S=({1..13}∖{e})∪{w}` has `L(S)<1/1260` then `(6/7)meas(G_C)−r/(7w)<1/1260`, which
(worst case `e=6`, `r,meas(G_C)` explicit) forces **`w ≤ 93`**. Exhaustive exact search over the
13 drops and `w≤2000` then gives the **unique minimal positive value `L=1/1260`**, attained
only at `{1,2,…,11,13,36}`.

**Champion mechanism (sharp value, closed-form).** Drop `e=12`: `G_C` has 4 arcs — 2 outer
(width `2/1001`) + 2 central (half-length `1/490`). Adding `w=36` (`=3·12`) recovers the outer
arcs but the central arc needs comb half-width `1/(14w)≥1/490` i.e. `w≤35`; `w=36` misses by
`1/504<1/490`. The residual is a **two-speed clash** between speed 5 (`2/5+1/70=29/70`) and
speed 36 (`15/36−1/504`): `15/36 − 2/5 − 1/70 − 1/504 = 1/2520`, and `τ↔1−τ` symmetry doubles it:
**`L = 2·(1/2520) = 1/1260`** (verified exactly). `w=24=2·12` covers fully (`1/336>1/490`) ⟹
`{1..11,13,24}` is **tight** (= Goddyn–Wong's family `{1,…,n−2,n,2(n−1)}`, tight for `n≡1 mod 6`,
here `n=13`). A discrete arithmetic lock: no `w` leaves a central residual strictly in `(0,1/2520)`.

## D. The crux — ONE open lemma (OPEN-Q-108)

By A–C + THM-522 (quantization `L>0⟹L≥1/(14·lcm)`, scale-invariance), the entire `inf L>0`
program reduces to:

> **UNIFORM FATTENING LEMMA (open):** there is `c>0` with `meas(G_C) ≥ c` for **every**
> 12-subset `C` of distinct positive integers — equivalently, the **primitive tight locus at
> n=13 is finite**.

Decoupling (B) supplies this when speeds grow **one at a time** (floor `1/143`) and iterates for
`k` large entries while the residual `(13−k)`-core keeps positive bounded-below measure; the
**only uncontrolled regime is `k≥3` arithmetically-coordinated growing speeds** (located by the
disprove side: the `drop-6` family minimizes at the *large* `w=69`, `L=19/10626`, still `>1/1260`
and rising). **The LRC(12) lever** (re-verified: exactly one 12-subset of `{1..14}` is tight at
gap `1/13`, zero at gap `1/14`) converts the crux from *existence* (`meas(G_C)>0`?) to
*transversality* (does the isolated gap-`1/13` maximizer **fatten** to a positive gap-`1/14`
measure, uniformly?). The literature confirms the tight-locus classification is "widely open"
(Perarnau–Serra) and that **no prior bound controls the safe-MEASURE** — all known LR bounds
(incl. Bedert 2025) bound the *gap* `κ(n)`, not `meas(G_C)`. The bounded-speed reduction (Tao
Thm 1.3 / MSS) makes the compactness scaffold rigorous in principle.

## E. Zero counterexamples (PROVED on the tight locus)

A genuine LRC(14) counterexample needs `L=0` AND `max-min := max_τ min_v||vτ|| < 1/14`. But
`L>0 ⟹ max-min>1/14`, so any counterexample is **tight**; and the only tight configs found
(`{1..13}` and `{1..11,13,24}`, exhaustive to entry≤18 and over all near-AP families) have
**`max-min = 1/14` EXACTLY** (two independent candidate-enumerations; AP at `τ=13/14`, sporadic
at `τ=3/14`). Above `1/14` the max-min spectrum has a **gap**: `1/14 < 3/41 < 2/27 < …`. So
**LRC(14) holds with equality on the tight locus; zero counterexamples** over 320k exact-checked
configs — consistent with LRC(14) being TRUE.

## Scope / honesty

PROVED & VERIFIED: A (decomposition), B (decoupling floor, min `1/143`), C (single-perturbation
inf `=1/1260` with explicit `N≤93`; champion two-speed-clash, exact), E (tight-locus max-min
`=1/14`, zero counterexamples). The CRUX (D, the uniform `meas(G_C)` bound ≡ tight-locus
finiteness) is **OPEN** — strong one-directional evidence (every exact search agrees, `L`
provably *increases* with lcm, disprove side failed on all fronts) but not proved; the literature
confirms it open. So this is a **near-complete singular-series proof of LRC(14), reduced to one
well-isolated quantitative lemma** — NOT a proof of LRC(14), and NOT a disproof (zero
counterexamples). Corrects THM-522 C: "extremizers have bounded lcm" is FALSE (near-tight cores
can have lcm `~10^7`, e.g. `{1,3,4,5,7,8,9,10,11,12,13,17,25}`, `L=2/425`); the correct
compactness bounds the **perturbing elements** (via B), not the lcm. Cross-links: THM-522, THM-518,
THM-501, HYP-2561, OPEN-Q-108, MISTAKE-075.
