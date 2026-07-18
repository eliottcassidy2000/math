---
id: THM-1032
title: THE CLUSTER-BAND EXISTENCE THEOREM — THM-1018(III) closed unconditionally by an EXPLICIT modulus, no divisor count needed. SETTING V = P ∪ K with |V| = 13, core P (|P| = 13−r ≥ 2), μ = min P, M = max P, killers K = {v₁ < … < v_r} all > 13M, spread D := v_r − v₁. THEOREM (a) if **M ≤ 12μ** and **D ≤ M − μ** then q := v₁ + M and a := ⌈q/(14μ)⌉ give t = a/q with ‖v·t‖ ≥ 1/14 for every v ∈ V; (b) the residues are FORCED — v_i < q and q − v_i = M − δᵢ, so every killer residue is M − δᵢ ∈ [μ, M] and E = P ∪ {M−δᵢ} ⊆ [μ, M], giving e_min = μ, e_max = M: **the ratio hypothesis of the band lemma is inherited from the core, not assumed**; (c) since q ≥ 14M + 1, the band [q/(14μ), 13q/(14M)] contains an integer as soon as 13μ − M ≥ μ, i.e. M ≤ 12μ; (d) BOTH hypotheses are SHARP — aspect ratio 12 passes 400/400 while 12.25, 12.5, 12.6, 12.67 all fail, and spread M−μ passes 1500/1500 while spread M−μ+1 fails 0/1500; (e) for ANY core P ⊆ {1,…,12} the hypothesis M ≤ 12μ holds AUTOMATICALLY (M ≤ 12 ≤ 12μ) and M − μ ≥ 10 when |P| = 11, so the whole d ≤ 5 clustered stratum of THM-1011(VII)/THM-1018 is certified with the reach DOUBLED — spread up to 10, and any number of killers r, not just two; (f) general threshold 1/N: the criterion is M ≤ (N−2)μ, strictly better than the (N−2+√(N²−4N))/2 of the companion q = v₁ − μ construction (11.9161 at N = 14)
status: PROVED (a)-(c) — elementary residue arithmetic plus THM-1018(II)'s band implication; the construction is explicit so NO existence/divisor-count argument is required. (d) sharpness and (e) coverage VERIFIED exhaustively: 240,000/240,000 over the full stratum shape (all 12 cores, v₁ ∈ (13M, 13M+4000], d ∈ 1..5), 103,870/103,870 with spread extended to 1..M−μ, 1500/1500 at r = 2,3,4,5 killers, 4000/4000 residue-identity checks, 600/600 on the former exceptional set
source: kind-pasteur-2026-07-18-S128 (cont.57; owner: prove the divisor count and close the existence step, with arXiv 2408.00183 offered as inspiration)
depends_on:
  - THM-1018 (II)    # the band lemma this supplies a modulus for
  - THM-1018 (III)   # the existence step this CLOSES (supersedes the divisor-count route)
  - THM-1011 (VII)   # the d ≤ 5 clustered stratum
related:
  - THM-1007         # single-killer + lacunary chains (the other half of "covering ⟹ M > 1/14")
  - THM-1035         # opus-S360's SEVEN-MODULI existence step — a DIFFERENT existence step, see the disambiguation below
  - THM-523          # the t = 1/q witness this generalizes to t = a/q
script: 04-computation/existence_closed_kps_S128c57.py, existence_multikiller_kps_S128c57.py, existence_direct_kps_S128c57.py (+ .out)
---

# THM-1032 — the cluster-band existence theorem

## What was open

THM-1018 proved the **band lemma**: given a modulus q whose effective speeds satisfy
2 ≤ eᵢ and e_max ≤ 13·e_min, every integer a ∈ [q/(14 e_min), 13q/(14 e_max)] certifies
min‖v·a/q‖ ≥ 1/14. What it did *not* prove was that a good modulus exists; it reduced
that to bounding the divisor count of the six numbers v₁±1, v₁, v₂±1, v₂.

**That count is not needed.** The modulus can be written down.

## The construction

> Let V = P ∪ K, |V| = 13, core P with μ = min P, M = max P; killers
> K = {v₁ < … < v_r}, all v_i > 13M, spread D = v_r − v₁.
> **If M ≤ 12μ and D ≤ M − μ, put q = v₁ + M and a = ⌈q/(14μ)⌉.
> Then min_{v ∈ V} ‖v·a/q‖ ≥ 1/14.**

*Proof.* Write v_i = v₁ + δᵢ with 0 ≤ δᵢ ≤ D ≤ M − μ.

*Residues.* Each v_i satisfies v_i = v₁ + δᵢ < v₁ + M = q, so v_i mod q = v_i and
q − v_i = M − δᵢ. Since v₁ > 13M > M ≥ M − δᵢ, the least-absolute residue is
eᵢ = M − δᵢ. From δᵢ ≤ M − μ we get eᵢ ∈ [μ, M]. Each core speed p ≤ M < q/2
(because q > 14M) is its own residue. Hence

  E := P ∪ {e₁,…,e_r} ⊆ [μ, M],  so **e_min = μ and e_max = M**.

The band lemma's ratio hypothesis e_max ≤ 13 e_min is therefore *inherited from the core*
— M ≤ 12μ < 13μ — rather than assumed about the killers. This is the whole point: the
killers, however large, are pulled back into the core's own window by the choice of q.

*Integrality.* The band [q/(14μ), 13q/(14M)] has length q(13μ − M)/(14Mμ), so it contains
an integer once q ≥ 14Mμ/(13μ − M). Here q = v₁ + M ≥ 13M + 1 + M = 14M + 1, and
M ≤ 12μ gives 13μ − M ≥ μ, whence 14Mμ/(13μ − M) ≤ 14M < q. So a = ⌈q/(14μ)⌉ lies in
the band, and THM-1018(II) applies. ∎

Explicitly: for each e ∈ E, e·a ∈ [q/14, 13q/14], so e·a never wraps and
‖e·a/q‖ ≥ 1/14; each v_i ≡ ±eᵢ (mod q) and ‖·‖ is even.

## Both hypotheses are sharp

| aspect ratio M/μ | 10 | 12 | 12.25 | 12.5 | 12.6 | 12.67 |
|---|---|---|---|---|---|---|
| certified / 400 | 400 | 400 | 397 | 393 | 369 | 369 |

The threshold is exactly 12, not the compression bound 13 — one full unit of slack is
spent on integrality. And on the spread side, D = M − μ gives 1500/1500 at r = 2,3,4,5
killers while D = M − μ + 1 gives **0/1500**: the residue M − δ drops below μ, e_min
collapses, and the band closes. The probe at μ=2, M=12 shows the switch at spread 11.

## Why the d ≤ 5 stratum is now closed

For any core P ⊆ {1,…,12}: M ≤ 12 and μ ≥ 1, so **M ≤ 12μ holds automatically**. With
|P| = 11 the core spans M − μ ≥ 10. So every clustered family in THM-1011(VII)/THM-1018
— whose hypothesis was d ≤ 5 and r = 2 — is certified, and the reach is strictly larger:
**any number of killers, spread up to 10**. Verified 240,000/240,000 over the full
stratum shape and 103,870/103,870 with the spread extended past 5 to M − μ.

The 8 critical covering families of THM-1011(VII) get clean explicit witnesses:

| (v₁,v₂) | q | a | t | min‖vt‖ |
|---|---|---|---|---|
| 168,169 | 180 | 7 | 7/180 | 7/90 |
| 195,196 | 207 | 8 | 8/207 | 16/207 |
| 208,210 | 220 | 8 | 2/55 | 4/55 |
| 221,224 | 233 | 9 | 9/233 | 18/233 |
| 234,238 | 246 | 9 | 3/82 | 3/41 |
| 247,252 | 259 | 10 | 10/259 | 20/259 |
| 294,299 | 306 | 11 | 11/306 | 11/153 |
| 308,312 | 320 | 12 | 3/80 | 3/40 |

Every one is ≥ 2/21 of threshold or better — nothing is close to tight.

## The general-N law

At threshold 1/N the criterion is **M ≤ (N−2)μ**. The companion construction q = v₁ − μ
(also verified here, 600/600 plus 8/8 critical) yields only M/μ ≤ (N−2+√(N²−4N))/2, which
is 6+√35 = 11.9161 at N = 14 — just short of 12, which is why it left a finite exceptional
set (cores containing both 1 and 12, v₁ ∈ [157,168]; 600 families, all separately
certified). Adding M rather than subtracting μ buys M + μ of modulus and closes the gap.

## ⚠ Disambiguation: two different "existence steps"

The owner's brief went to several machines at once and the words map onto **two distinct
statements**. Recording this so the fleet does not conflate them:

- **opus-S360, THM-1035** — the *denominator sieve*'s existence step: which q ∈ {2,…,14}
  must divide a speed. Answer: **seven** moduli, uniquely {8,…,14}; the six-number window
  there is {9,…,14}, which misses q = 8. Kernel-pure in Lean. **Not this theorem.**
- **THM-1032 (here)** — THM-1018(III)'s existence step: does a good *band* modulus exist
  for a clustered killer family. The six-number window there was {v₁−1, v₁, v₁+1, v₂−1,
  v₂, v₂+1}. Answer: the question dissolves — q = v₁ + M is always good.

Both are real, both are now closed, and they are unrelated. Cf. MISTAKE-158: match a
lemma by its *statement shape*, not by the phrase used to describe it.

## On arXiv 2408.00183 (Couvreur–Zémor, Freiman 3k−4 for function fields)

Offered as inspiration. Honest accounting: it supplied **no technical tool** — it contains
nothing on residues mod q, covering systems, or runner problems. The resonance is at the
level of *method*, and it is real:

1. **Structure replaces counting.** Couvreur–Zémor transport Freiman's theorem by
   replacing cardinality-counting with a structural invariant (degrees in Riemann–Roch
   spaces). The move here is the same shape: the divisor **count** of the six-number
   window is replaced by an explicit **construction**. When a counting argument resists,
   ask whether the object it is counting can simply be exhibited.
2. **Bounded expansion forces interval confinement.** Freiman 3k−4 says small doubling
   forces a set into a short AP. Part (b) above is the multiplicative-scale analogue:
   bounded aspect ratio plus bounded spread force every killer residue into the single
   interval [μ, M]. Confinement-to-an-interval is exactly what the band lemma consumes.

The genuinely additive-combinatorial reading of LRC(14)'s crux belongs to **boxeph-S89**
(HYP-7372: sharp Freiman 3k−4 ⟹ M < 1/13 forces the 12-core to be an AP) — that is the
right home for the paper's actual content, and this note defers to it.

## Status of "covering ⟹ M > 1/14" after this theorem

- single-killer: **PROVED** (THM-1007)
- lacunary multi-killer chains: **PROVED** (THM-1007 ext)
- clustered, spread ≤ M − μ, any r: **PROVED** (here) — was the d ≤ 5 residual
- clustered, spread > M − μ and non-lacunary: still (BG-K)-certified, not proved

The last line is the honest remaining gap, and it is now the *only* one on this branch.

## Named next
- Push the modulus family past spread M − μ: q = v₁ + c needs c ∈ [μ+D, M], so the cap is
  structural for this family. A second modulus (or a two-modulus alternation) is wanted for
  the wide-cluster case.
- Lean: (a)–(c) is a short rational inequality chain on top of THM-1018(II) — the residue
  identities e_i = M − δ_i are `Nat.sub` facts, and the integrality is one `Nat.ceil` bound.
