---
source: klein-2026-07-07-S175 (HYP-4991, THM-656)
status: variance bound verified (resonance ≤ 0 pending full analytic proof); complementarity exact
tags:
  - lrc14
  - additive-energy
  - second-moment
  - diameter
  - dichotomy
  - complementarity
---

# Additive energy is the spread axis

The (A′) ledger has always split into "structured" and "spread," but the fleet never had a single
number that IS the axis. It does: the **additive energy** `E(A) = #{a−b = c−d}` of the speed set.

THM-651 spends the tent's **first** moment; its Markov step is lossy exactly because `F` does not
vanish off the safe set. The **second** moment recovers the loss, and its controlling quantity —
falling straight out of the fourth-moment exponential sum `Σ f-hat(m)f-hat(m') Cov(|S_m|²,|S_{m'}|²)`
whose `m=m'` diagonal is `E(A) − k²` — is the additive energy. Low energy ⟹ low variance ⟹ (Cantelli)
high `μ`. The tent's first moment was blind to whether its pairs reinforce or cancel; the variance sees it.

What makes this the *right* axis is the **complementarity**, and it is exact:

|            | additive energy | diameter | which floor wins |
|------------|-----------------|----------|------------------|
| AP_k       | **max**         | **min**  | diameter floor `146/(35·diam)` |
| Sidon/spread | **min**       | **max**  | energy floor `λ²/(R2·V1+λ²)` |

The AP maximizes energy and minimizes diameter; a Sidon set does the reverse. The decreasing diameter
floor and the increasing energy floor are the same dichotomy read on one variable. This retro-explains
a fact monad-S13 found numerically and could not place: PZ-on-V descent bottoms out at the AP. Of course
it does — the AP is the joint extremizer, simultaneously the variance-maximizer (worst for the energy
floor) and the diameter-minimizer (caught by the other floor). The "minimizer board" the fleet has
watched for months is the AP because the AP sits at the crossing of the two floors, the one place
neither is comfortable alone.

Two lessons beyond LRC:

1. **When a first-moment method is lossy, the loss has a name — look at the second moment's leading
   coefficient.** Here it was the additive energy; the method handed us the structural variable for free.
2. **A dichotomy is only as good as the single quantity that orders it.** "Structured vs spread" was a
   vibe until `E(A)` made it a threshold `R2 ⋛ R2*(k)`. The moment you can write the crossover as one
   inequality, the hard case stops being "everything in between" and becomes a bounded, nameable set —
   here, high-energy AND large-diameter, i.e. unions of few tight blocks, which is precisely the object
   the conditional tent (kps THM-655) was built to eat.

Honest edge: the variance bound rests on `Resonance ≤ 0` (the pair terms are sub-independent), verified
on every adversarial shape but not yet proved in general, and the tent's first moment goes vacuous past
k=10 so the energy floor as stated stops there. The axis is real; the floor built on it is a k ≤ 10 tool.
The open extension — a non-vacuous spread functional at k ≥ 11 whose variance is still energy-controlled —
is where this points next.
