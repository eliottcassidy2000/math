---
id: THM-946
title: THE EXACT SUPPORT DECOMPOSITION AND THE ABSOLUTE-BUDGET SHARPNESS LAW — (I) THE M_s MASSES ARE EXACT RATIONALS, NOT ESTIMATES: define recursively over subsets A (|A| = 2..5) M(A) := μ(∩_{v∈A}B_v) − (2λ)^{|A|} − Σ_{S⊊A,|S|≥2} M(S)·(2λ)^{|A|−|S|} (each μ(∩) an exact integer-grid sweep on 14·lcm(A)); then M_s = Σ_{|A|=s}M(A) and THM-935's identity B₅ = 2052/16807 + Σ_s (−1)^s E_s M_s holds as EXACT RATIONAL EQUALITY against the direct depth-spectrum B₅ — refereed on three real 13-packets (certified {307..3310}, opus-30Z, random): identity EXACT all three. For CONCRETE packets this replaces all tail estimation (THM-944's T_s remain only for the UNIVERSAL statement) — the relation side of the theory is now exactly computable end to end; (II) THE SHARPNESS LAW (the answer to the signed-vs-absolute question at the E_s level): the slack factor |Σ_s|E_s||M_s|| / |Σ_s(−1)^s E_s M_s| = 1.00 on ALL three packets — because the masses are SYSTEMATICALLY NONNEGATIVE at s ≥ 3 (M₂ alone can be slightly negative: −6×10⁻⁵ on the random packet) and the sign pattern makes every term a debit except the tiny +E₂M₂ credit. NO signed cancellation exists at the support level: the absolute budget (codex's LRCB5RelationBudget frame) is ALREADY SHARP there, and any sharpening (boxeph's near-orthogonality) must operate INSIDE the M_s — per-subset — where the top masses are also all positive. STRUCTURAL READING + NAMED HYP: relation resonances systematically INCREASE deep coverage (the arcs' overlaps are positively-associated) — conjecture M(A) ≥ 0 for |A| ≥ 3 (FKG/positive-association flavor; the s = 2 exception is real: pair repulsion exists); (III) THE EXACT AUTOPSY INSTRUMENT: per-subset masses name failure culprits with exact fees — opus-30Z's B₅ = −0.6016 decomposes as M₃ = 0.2318 + M₄ = 0.3650 + M₅ = 0.5059 led by the SCHUR TRIPLE {420,450,870} (420+450 = 870) at +0.01250 and the AP clusters {450,510,570}, {450,570,690} at +0.0075 each — the exact fee data death-star's deviation ledger (THM-938/939/942/943/945) prices axiomatically, now measured; the certified packet's largest subset mass is 0.0022 (16× smaller), with M₂ = +0.00031 (its ratio content near-neutral)
status: (I) exact identity referee 3/3 packets (rational equality; script cited); (II) sharpness measured exactly (slack 1.000-1.001); the M(A) ≥ 0 (s ≥ 3) conjecture named with the observed exception structure at s = 2; (III) autopsies exact
source: kind-pasteur-2026-07-17-S128 (cont.39; owner: highest-leverage angles toward finishing LRC(14))
depends_on:
  - THM-935 (the E_s identity this makes exactly computable)
  - THM-930 (depth-spectrum sweep = the direct side of the referee)
related:
  - THM-944 (universal tails — now needed ONLY for the universal statement; concrete packets are exact)
  - death-star THM-938/939/942/943/945 (the fee ledger — (III) is its measured input)
  - codex LRCB5RelationBudget/LRCB5NormalizedBridge (the absolute frame — (II) says it is sharp at E_s level)
  - boxeph LEM-031 (character-mass factorization — the inside-M_s signed structure lives there)
script: 04-computation/exact_support_decomposition_kps_S128c39.py -> 05-knowledge/results/exact_support_decomposition_kps_S128c39.out
---

# THM-946 — exact support decomposition; the sharpness law

## (I) The computation

M(A) is the excess overlap of the subset A after removing the independent prediction and
all lower-support excesses — Möbius inversion over the subset lattice, with every
ingredient an exact rational sweep. The identity referee (all rational equalities):

| packet | B₅ direct | B₅eq + signed | M₂ | M₃ | M₄ | M₅ |
|---|---|---|---|---|---|---|
| certified | +0.082053 | EXACT ✓ | +0.000311 | +0.006682 | +0.023545 | +0.030061 |
| opus-30Z | −0.601603 | EXACT ✓ | +0.000098 | +0.231821 | +0.365015 | +0.505866 |
| random | +0.088381 | EXACT ✓ | −0.000060 | +0.002007 | +0.019141 | +0.027255 |

## (II) Sharpness

Signed total vs absolute bound: 0.040039 vs 0.040083 (certified), 0.72370 vs 0.72371
(30Z), 0.033711 vs 0.033711 (random). The absolute budget loses nothing at this level.
Conjecture (named): **M(A) ≥ 0 for |A| ≥ 3** — positive association of arc overlaps;
the pair level genuinely admits repulsion (M₂ < 0 observed).

## (III) The autopsy

The instrument outputs, per packet, the ranked exact subset fees. The 30Z verdict:
its own gcd-30 dilate structure contributes nearly nothing at support 2 (+0.0001) — the
kill is the ADDITIVE CLUSTER: one Schur triple and two APs inside {420..870} carry the
lead fees, and the mass cascades through supports 4–5 as those clusters embed. This is
the measured form of THM-926's diagnosis, and the exact input for the fee-ledger route.

## Named next
- The M(A) ≥ 0 (s ≥ 3) conjecture: prove (positive-association/FKG on arc systems) or
  find the counterexample census.
- Wire the exact fees into death-star's ledger (replace axiomatized fee floors by
  measured values where packets are concrete).
- The inside-M_s signed structure via boxeph's character factorization (the only
  remaining sharpening room).
