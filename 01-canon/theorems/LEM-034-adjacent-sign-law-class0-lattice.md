---
id: LEM-034
title: THE ADJACENT-CLASS SIGN LAW AND THE CLASS-0 LATTICE (LEM-033's N₀ named check RESOLVED — structural, by ATTRIBUTION, not full blocking). (A) SIGN LAW (universal, any cluster, any s): an endpoint whose attributed owner has boundary class c ≡ s (mod 7) is ALWAYS an exit (sign −1: just after it that owner sits in section s), and c ≡ s+1 is ALWAYS an entry (just before, the owner sits in section s) ⟹ per owner N_s ≤ 0 ≤ N_{s+1} signed. Referee: 1080 endpoints, 6 clusters × 7 sections, 159 class-s endpoints all exits + 159 class-(s+1) all entries, 0 violations. COROLLARY (reflection bijection): x ↦ 1−x maps class-s exits of R_s to class-(s'+1) entries of R_{s'}, s' = 6−s (indices j ↦ 7e−j, sections σ ↦ 6−σ) — hence the exact global count match 159 = 159. (B) THE COINCIDENCE LATTICE + ATTRIBUTION LEMMA: at a class-0 crossing x = k/e of a 7-full owner e = 7m, the co-boundary owners are EXACTLY {f ∈ E : m | fk}; every f = mt ∈ E is co-boundary at EVERY k (index tk). The endpoint machinery attributes min(owners), so if E contains mt with t ≤ 6 (< e), owner e can NEVER receive a class-0 endpoint: N₀^(e) = 0 IDENTICALLY (balanced: 21←{12,15}, 28←{12,20}, 35←{15,20,30}; near-AP: 14←{8,10,12}); with no small co-multiple (two-large 56/84, family70 70) attribution permits N₀ < 0 (measured −4, −6, −6). (C) BLOCKING LEMMAS + BOUNDARY CENSUS (convention-independent layer): (B1) a co-lattice runner with index fk/m ≡ 1 (mod 7) sits in section 0 just before x — both sides dead; (B2) at 7 | k with < 5 non-multiples of m, occupancy of {1..5} is impossible — dead. Census (exact): most crossings genuinely dead (e=21: 10 B1 + 1 other-zero + 9 occupancy of 21; e=28: all 28 dead; e=35: 33/35; e=14: 13/14), but SURVIVORS EXIST on balanced/near-AP (21@k=11, 35@k=22, 14@k=12) — true R₀ exits, re-attributed to smaller co-lattice owners. Survivor ⟺ endpoint referee exact; attribution counts reproduce N₀ exactly (56: 5 survivors/4 attributed; 84: 8/6; 70: 6/6, at k ≡ 1 mod 10, k ≠ 1 — observed pattern)
status: PROVED ((A) 2 lines + reflection 3 lines; (B) one-line integrality + min-convention; (C) B1/B2 short) + REFEREED EXACT (6 clusters, all sections, all crossings of all 7-full owners; survivors cross-checked against the endpoint list; N₀ reproduced from attribution counts)
source: boxeph-2026-07-17-S63 (owner directive: the N₀ structural check; keep proving little statements)
depends_on: [LEM-033 (the named check), S25 endpoints machinery (the min-attribution convention), S26 owner_data]
related: [LEM-030 (σ_e owner imbalance — now decomposed with forced-sign classes: family70's owner 70 has N₁ = +15 inside σ = +4), THM-886 (N_c comb amplitudes)]
script: 04-computation/lrc14_class0_blocking_boxeph_S63.py -> 05-knowledge/results/lrc14_class0_blocking_boxeph_S63.out
---

# LEM-034 — the adjacent-class sign law and the class-0 lattice

## (A) The sign law (universal)

If the attributed owner's boundary index at an endpoint has class c ≡ s, the
owner enters section s just after — so the right side is outside R_s and the
endpoint is an exit (−1). If c ≡ s+1, the owner sits in section s just before
— the left side is out and the endpoint is an entry (+1). ∎ Hence
**N_s^(e) ≤ 0 ≤ N_{s+1}^(e)** for every owner, cluster, and s.

**Reflection corollary.** x ↦ 1−x sends {ex} to 1−{ex}: sections σ ↦ 6−σ,
so R_s ↦ R_{6−s}; boundary indices j ↦ 7e−j (classes c ↦ −c mod 7); entries
and exits swap. A class-s exit of R_s (c = s) maps to a class-(−s) entry of
R_{6−s}, and −s ≡ (6−s)+1 (mod 7): exactly the adjacent-entry class. ∎
The global referee count — 159 class-s exits, 159 class-(s+1) entries across
6 clusters × 7 sections — is this bijection made visible.

## (B) The coincidence lattice and the attribution lemma

At x = k/e, e = 7m: runner f is at its own section boundary iff 7f·x ∈ Z iff
**m | fk** (referee: exact). In particular f = mt ∈ E is co-boundary at every
k, with index tk. The machinery attributes endpoints to min(owners), so:

> **If E contains a multiple mt of m with t ≤ 6, then N₀^(e) = 0 identically**
> — every surviving class-0 boundary of e is re-attributed to a smaller
> co-lattice owner. This is the whole mechanism of the S62 observation.

Discriminator across the six clusters: balanced and near-AP have small
co-multiples under every 7-full owner (forced 0); two-large (m = 8, 12) and
family70 (m = 10) have none — and indeed N₀ = −4, −6, −6 there (negative, as
(A) forces). CONVENTION CAVEAT: this layer is a fact about the min-attribution
convention (which S26/S62 data inherit consistently); the layer below is
intrinsic.

## (C) The boundary census (intrinsic layer)

Is x⁻ ∈ R₀ at all? Blocking lemmas: **(B1)** a co-lattice runner with index
fk/m ≡ 1 (mod 7) sits in section 0 just left of x — both sides dead (right
side is always dead: e enters section 0). **(B2)** at 7 | k every co-lattice
runner sits at the top of section 6, so with < 5 non-multiples of m the
sections {1..5} cannot all be covered — dead. Census (exact, all crossings):

| owner | crossings | B1 | other-zero | occupancy | survivors |
|---|---|---|---|---|---|
| balanced 21 | 21 | 10 | 1 | 9 | **1** (k=11) |
| balanced 28 | 28 | 11 | 2 | 15 | 0 |
| balanced 35 | 35 | 16 | 1 | 17 | **1** (k=22) |
| near-AP 14 | 14 | 9 | 1 | 3 | **1** (k=12) |
| two-large 56 | 56 | 9 | 20 | 22 | **5** |
| two-large 84 | 84 | 11 | 31 | 34 | **8** |
| family70 70 | 70 | 13 | 27 | 24 | **6** (k ≡ 1 mod 10, k ≠ 1) |

Every survivor is a genuine endpoint (exit) of R₀ — the balanced/near-AP ones
exist but belong to smaller owners by (B); the sparse-lattice ones stay with
the 7-full owner, and the attribution counts reproduce N₀ exactly (5→4, 8→6,
6→6).

## Resolution of the named check

**N₀ = 0 for 7-full owners on balanced/near-AP is STRUCTURAL — forced by the
coincidence lattice through the attribution convention — while the intrinsic
boundary phenomenon is "mostly blocked, rarely surviving" (B1 + occupancy +
cluster-specific section-0 hits), with survivors re-labeled.** The correct
intrinsic invariant for downstream use is the survivor count, not N₀.

## Evidence log
- [x] sign law: 1080 endpoints, 0 violations; N_s ≤ 0 ≤ N_{s+1} all owners
- [x] reflection bijection: 159 = 159 exact
- [x] lattice rule own = {f : m | fk} exact; forced owners all N₀ = 0
- [x] census complete (255 crossings, 7 owners, 4 geometries); survivors ⟺ endpoints
- [x] named check RESOLVED (LEM-035, S64): the multiplication-permutation law — clean columns 6r < M, rescued column r = M/6, sporadics (M−2,5) at M = 11,12; M = 11 prediction k = 64 confirmed
- [ ] named: σ_e forced-class decomposition as a LEM-030 refinement (N₁ = +15
      inside σ = +4 at family70/owner 70 — the imbalance lives in the forced classes)
