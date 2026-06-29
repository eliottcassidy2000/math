---
id: HYP-3548
title: The proof/disproof boundary of LRC(14) has TWO razor-thin lines at different scales -- the GAP line (covering sets cluster at M>=7/89, +10%, NOT thin; quantized by binding-pair denominators) and the MEASURE-FLOOR line (R'>0, thin only for UNBOUNDED covering sets); a disproof must be a primitive, covering, UNBOUNDED, even-heavy 13-set with a resonant 'lonely-poor' descended level; the apex-7/units extremal {1..13} is NON-COVERING and OFF the critical path; my recent metagraph/Ky-Fan/perfect-number work is sigma-ODD WITNESS-side, ORTHOGONAL to the sigma-EVEN floor gatekeeper, except as the apex-7 BASE of the 2-adic descent (Heegner SOS)
status: SYNTHESIS + boundary computation. VERIFIED: M({1..11,13,84})=7/89 exactly; single-large family min at L=84; covering sets cluster >>1/14. Integrates Explore-agent map (THM-523/580, HYP-3129/3415, OPEN-Q-108). Reorganizes proof targets; does NOT close the gatekeeper.
source: mac-mini-2026-06-29-S15
related:
  - THM-523    # covering reduction: LRC(14) <=> M>=1/14 for primitive COVERING sets
  - THM-580    # 2-adic parity descent (peels even speeds -> per-level rho_j)
  - HYP-3415   # the covering floor R'>0 = the critical path / OPEN-Q-108
  - HYP-3129   # certified R'>=0.642 (bounded); exact-low + Parseval-tail technique
  - HYP-3547   # the three pillars = the three stages of ONE descent (apex-7 base)
  - THM-581    # floor is even-category; Borsuk-Ulam is witness-side (orthogonal)
  - HYP-3535   # S75e Heegner Fejer-Bochner SOS (the base-level certificate)
results:
  - 04-computation/lrc14_disproof_boundary_margins_macmini_20260629.py
  - 04-computation/lrc14_covering_set_boundary_macmini_20260629.py
  - 05-knowledge/results/lrc14_covering_set_boundary_macmini_20260629.out
---

# HYP-3548 -- the two razor-thin lines and the reorganized proof targets

## A disproof, precisely
A primitive COVERING 13-set `S` (gcd 1, a multiple of every `q in 2..14`) with `M(S) < 1/14`
strictly -- the lonely set empty, the danger combs (each measure `1/7`) covering `[0,1)`. Non-covering
sets are closed by the `t=1/q` witness (THM-523), so a disproof MUST be covering.

## Two razor-thin lines (the key reframe)
- **GAP line (M):** covering sets cluster at `M >= 7/89 = 0.07865` (**+10.1%** above `1/14`). VERIFIED
  `M({1..11,13,84}) = 7/89` exactly (`84 = lcm(12,14)`); the single-large family `{1..11,13,L}` is
  tightest at `L=84` and only loosens for larger `L`. Quantization: `M = (v_a q - v_b p)/(v_a +- v_b)`,
  denominator a binding-pair sum/difference -> a discrete lattice the covering constraint keeps clear
  of `1/14`. **The disproof is FAR in gap-value.**
- **MEASURE-FLOOR line (R'):** `meas(L) = R'*(m_R m_Q)`, `R' > 0 => M >= 1/14`. Certified `R' >= 0.642`
  for BOUNDED covering families (60% margin); **UNBOUNDED open**. This is the only thin line; the single
  gatekeeper (OPEN-Q-108) is `rho_j >= c > 0` on the 2-sheet 2-adic descent (THM-580).

## Necessary conditions for a counterexample (a short list)
Primitive + covering + **UNBOUNDED** speeds (bounded certified `R'>=0.642`) + **even-heavy** (the binding
obstruction is 2-adic, S259) + a **resonant/lonely-poor descended level** (where the per-level
Cauchy-Schwarz goes non-positive). A disproof lives only in the un-certified tail of the measure floor.

## Honest placement of the tournament thread
The floor is `sigma`-EVEN (existence invariant under `t->1-t`, SOS-provable without Borsuk-Ulam). My
recent metagraph / Ky-Fan / signed-cycle-index / perfect-number work is `sigma`-ODD WITNESS-side
(Redei, odd index, klein's R-odd Betti) and is **orthogonal to the floor gatekeeper**. Its one
critical-path contact: the 2-adic descent bottoms out on an ALL-ODD face = the apex-7 cyclotomic object,
where the gentlest-cyclotomic Heegner SOS (`F_7`, `Q(sqrt-7)` h=1) is the base-level `rho` certificate.
The three pillars (HYP-3547) are NOT three attacks but the three STAGES of one descent:
**2-adic peel (Mersenne) -> per-level SOS (Heegner) -> Borsuk-Ulam witness backstop (3 mod 4).**

## Reorganized proof targets
1. **PRIMARY (gatekeeper):** `rho_j >= c > 0` uniformly via HYP-3129's exact-low + Parseval-tail PER
   2-sheet descent level (simpler than the monolithic 14-sheet: 2 sheets, smaller sets).
2. **Use the 10% gap as slack:** covering sets sit at `M >= 7/89`, so target `R' > 0`, not `R' >= 0.642`;
   the gap margin absorbs Cauchy-Schwarz loss.
3. **Fragile point has a silver lining:** lonely-poor descended sets = the structured/resonant ones =
   exactly where cyclotomic SOS is STRONGEST. Deploy Cauchy-Schwarz on generic strata and SOS on the
   resonant ones (complementary, not uniform).
4. **Gap-quantization route (alternative):** bound the binding-pair denominator `D = v_a +- v_b` for
   covering sets so `M = j/D` with `14j >= D` forced -- converting the analytic floor into a
   number-theoretic statement about speed differences (= tournament arcs), where the Redei/OCF thread
   could finally touch the critical path.

## What this buys
A precise, short target: ONE even-category inequality `rho_j >= c` on 2-sheet descents, with a 10% gap
cushion and an apex-7/Heegner base case. And clarity on what is NOT on the path (the entire sigma-odd
tournament scaffolding), so effort can be aimed at the gatekeeper.
