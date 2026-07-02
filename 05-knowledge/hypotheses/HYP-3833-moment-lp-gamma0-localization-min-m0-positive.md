---
id: HYP-3833
title: THE Gamma_0(N) CONGRUENCE LOCALIZATION = THE MULT-BY-n PHASE COORDINATE, AND IT RESCUES THE STALLED MOMENT LP -- min m0 > 0 under localization. The covering-min moment method (first-moment/union bound) STALLS: total danger 2(n-1)*M_C ~ 1.99 > 1 => the unlocalized certified loneliness m0 = 1 - 2(n-1)M_C = -0.989 < 0 (the SOS/Fourier stall, HYP-3791). KEY IDENTIFICATION: localizing with the Gamma_0(14) congruence = passing to the MULT-BY-n=14 PHASE-RESIDUE coordinate p(v)=n*v mod Phi6 (S68/HYP-3800; n=14=level, mult-by-n = the Hecke/level multiplier, HYP-3706). In these coordinates the runner cloud is the AP(step n)+antipodal killer (HYP-3813) and the observer's clearance = n/Phi6. VERIFIED (n=14, Phi6=183): (A) unlocalized union m0 = -0.989 < 0 (stall); (B) LOCALIZED min m0 = n/Phi6 = 0.0765 > 0 (CERTIFIED, >= 1/n), = the max-min clearance over observer offset (observer 0 optimal) = the equioscillation LP-dual on the binding pair {+n,-n} (linprog confirms). So the congruence symmetry is EXACTLY what turns the negative union bound into a positive certificate: min m0 > 0 under Gamma_0(N)-localization. HONEST: this certifies the CONSTRUCTION's covering-min > 0 (the achievable direction), re-reading the S68/HYP-3813 clearance as a localized moment LP; it does NOT yet certify the min over ALL configs (LRC's remaining hard direction). The lesson: the moment method's stall is a coordinate artifact -- in the Gamma_0(N)/mult-by-n coordinate the positivity is manifest
status: VERIFIED (computation) as a localized-moment-LP restatement; the LRC min-over-configs remains OPEN. VERIFIED (moment_lp_gamma0_localization_klein.py): unlocalized union bound m0 = 1-2(n-1)M_C = -0.989 (stall); localized (mult-by-n phase) min m0 = n/Phi6 = 0.0765 > 0 >= 1/n; max-min over observer offset = clearance at offset 0; linprog LP-dual on {+-n} confirms. HONEST: a coordinate/localization re-reading of the known S68/HYP-3813 clearance (certifies the construction, not the min-over-configs); the value of the finding is the IDENTIFICATION Gamma_0(N)-localization = mult-by-n phase, and that it rescues the moment stall.
source: klein-2026-07-01-S86
depends_on:
  - HYP-3800   # S68: the phase-residue coordinate p(v)=n*v mod Phi6 (= the Gamma_0(N) localization)
  - HYP-3813   # S80: the AP-cloud clearance = n/Phi6 (the localized min m0)
related:
  - HYP-3791   # the SOS/Fourier stall (unlocalized) that localization rescues
  - HYP-3831   # S85: the covering-min game LP-dual (this is its localized moment form)
  - HYP-3706   # mult-by-n = the level/Hecke multiplier (why mult-by-n = Gamma_0(N))
external: Gamma_0(N) congruence subgroup / level-N Hecke structure; moment/SOS relaxations; Koksma-Hlawka
results:
  - 04-computation/moment_lp_gamma0_localization_klein.py
  - 05-knowledge/results/moment_lp_gamma0_localization_klein.out
---

# HYP-3833 — Γ₀(N) localization = mult-by-n phase, and it rescues the moment stall

## The test
The covering-min moment method (first-moment / union bound) **stalls**: total danger `2(n-1)·M_C ≈ 1.99 > 1`, so
the unlocalized certified loneliness `m0 = 1 − 2(n-1)M_C = −0.989 < 0` (the SOS/Fourier stall, HYP-3791).

## The identification
**Localizing with the Γ₀(14) congruence = passing to the mult-by-`n`=14 phase-residue coordinate** `p(v) = n·v
mod Φ₆` (S68/HYP-3800). Here `n=14` is the level and mult-by-`n` is the Hecke/level multiplier (HYP-3706). In
these coordinates the runner cloud is the **AP(step n) + antipodal killer** (HYP-3813).

## Verified (n=14, Φ₆=183)
- **(A) Unlocalized** union bound: `m0 = −0.989 < 0` (stall).
- **(B) Localized** (mult-by-n phase): the cloud is `{14,28,…,168, 169}`; observer clearance
  `= min_v dist(nv mod Φ₆)/Φ₆ = 14/183 = 0.0765`. The max-min over observer offset is achieved at offset 0
  (observer optimal); the linprog LP-dual on the binding pair `{+n,−n}` confirms `m0 = n/Φ₆`.
- **min m0 (localized) = n/Φ₆ = 0.0765 > 0** (and `≥ 1/n`). **CERTIFIED.**

## Honest reading
The congruence symmetry is exactly what turns the negative union bound into a positive certificate: `min m0 > 0`
**under** Γ₀(N)-localization. This certifies the **construction's** covering-min `> 0` (the achievable
direction), re-reading the S68/HYP-3813 clearance as a localized moment LP. It does **not** yet certify the min
over **all** configs (LRC's remaining hard direction). The keeper is the identification — **Γ₀(N)-localization =
the mult-by-n phase coordinate** — and that in that coordinate the moment method's stall dissolves: the
positivity is a coordinate-manifest fact (an AP of step n leaves clearance n), not a failure of the method.

## Net
`min m0 > 0` under Γ₀(14)-localization, because Γ₀(14)-localization *is* the mult-by-n phase coordinate where the
AP cloud has clearance n. The moment stall is a coordinate artifact; the congruence coordinate makes the
covering-min positivity manifest (for the construction).
