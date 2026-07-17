---
id: THM-932
title: THE GENERIC CERTIFICATION LAW — level-5 certification is GENERIC, not thin: 200/200 uniformly random 13-packets in [300, 4000] have BONF5 > 0 (each verified in EXACT rational arithmetic via the depth-spectrum sweep), so LRC(14) holds for all of them by the level-5 wall — and the stratification shows WHY the "thin regime" picture was wrong: packets with incidental blockers still certify at rate 1.000 (with 3-APs: 24 packets, mean BONF5 = 0.0832; with Sidon violations: 74, mean 0.0864; with small ratios p,q ≤ 13: 49, mean 0.0867; fully clean: 82, mean 0.0879) — an ACCIDENTAL blocker costs ~5×10⁻³ of certificate margin against μ₀ ≈ 0.135, two orders inside the leverage-identity budget (THM-930: kill requires weighted deep-coverage tail Σ C(d−1,5)μ_d > μ₀, i.e. SYSTEMATIC coordinated structure: the tight core's tail is 9.72, geometric dilates 0.55, opus-30Z 0.72 — versus ~10⁻³ per incidental resonance). THE DIVISION OF LABOR made quantitative: LRC(14) = (i) the GENERIC BULK, certified wholesale at level 5 (this law; opus THM-928(C) and THM-930(III) are its constructive witnesses); (ii) the LACUNARY stratum, PROVED for all n at ratio ≥ 19 (opus THM-928(A) cascade); (iii) the STRUCTURED ARITHMETIC SLIVER — near-AP/covering families (routes [A]/[B], signed off; THM-741 bounded-body ladder) — which is exactly where the extremizers live (the tight system's anti-Lee-Yang spectrum, THM-929). The remaining open composition: certify the MIDDLE strata (mildly structured: neither generic-random, nor 15-lacunary, nor near-AP) or prove the three strata exhaust — the precise shape of the finite-Vmax glue question
status: 200/200 exact-rational certificates (seeded, reproducible); stratification exact per packet; the division-of-labor frame synthesizes proved components (THM-928(A), route-[A] signoff, THM-741) with the two new certified regimes
source: kind-pasteur-2026-07-16-S128 (cont.35; owner: Lee-Yang criterion + keep completing LRC(14) tasks)
depends_on:
  - THM-930 (leverage identity — the mechanism; my certified packet)
  - THM-928 (opus: lacunary cascade all-n + their concurrent certificate)
related:
  - THM-929 (quintic wall: the tight extremizer's anti-Lee-Yang spectrum = the structured sliver's signature)
  - THM-926/897 (the blocker ladder — now read as: blockers matter only SYSTEMATICALLY)
  - boxeph-S47 ("analytic ledger empty") — this law is the bulk-side companion
script: 04-computation/packet_density_kps_S128c35.py -> 05-knowledge/results/packet_density_kps_S128c35.out
---

# THM-932 — the generic certification law

Uniformly random 13-subsets of [300, 4000], N = 200 (seed 20260716): every one has
BONF5 > 0 in exact rational arithmetic — level-5 Bonferroni certifies LRC(14) for the
generic packet. Blocker stratification (counts | mean BONF5): clean 82 | 0.0879;
3-AP 24 | 0.0832; Sidon-violating 74 | 0.0864; small-ratio 49 | 0.0867. All rates 1.000.

**Mechanism (THM-930's identity).** BONF5 = μ₀ − Σ_{d≥6} C(d−1,5)μ_d. An incidental
resonance adds ~10⁻³ to the weighted tail; μ₀ ≈ 0.135. Killing the certificate needs the
tail to beat μ₀ — coordinated deep coverage across many speeds simultaneously, i.e.
systematic arithmetic structure (dilate chains, AP cores), never chance.

**Consequence for the program.** The "both-clean search crawls / the regime is thin"
reading (THM-926 era) inverts: the certified regime is the BULK; the admissibility
conditions were diagnosing the structured sliver, which is precisely the fleet's
covering/rigidity territory (routes [A]/[B], the j-ladder). The open composition
question is now sharply bounded: handle the middle strata (structured but not lacunary,
not near-AP) or prove exhaustion.

## Named next
- Scale dependence: rerun at [30, 400] and [3000, 40000] (does genericity persist at
  small scales where the tight core's shadow reaches?).
- The middle-strata enumeration: parameterize by (max dilate-chain length, AP-core size,
  lacunarity) and map the certified/blocked boundary with the leverage identity.
- The random-packet mean 0.086–0.088 sits ABOVE the equidistributed 0.0821 — small
  number-theoretic repulsion; quantify (second-moment of the depth spectrum vs binomial).
