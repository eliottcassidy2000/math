# Unknotting number one: exact status + in-repo decider design (klein-S691, 2026-07-28)

**Prompt (owner puzzle bundle / Epoch FrontierMath open problem):** an
algorithm taking a PD diagram (≤ 100 crossings), returning True iff the
represented knot has unknotting number 1, under an hour on a typical
machine. (https://epoch.ai/frontiermath/open-problems/unknotting-number)

## Status ledger (checked 2026-07-28)

- **Decidability of u(K)=1 is OPEN.** Not known decidable, per the problem
  page itself and the literature; a provably-correct always-terminating
  decider would be a genuine theorem, not an engineering artifact.
- **Known theory:** u=0 (unknot recognition) decidable, in NP ∩ co-NP
  (Lackenby for co-NP); *diagrammatic* unknotting variants are NP-hard
  (de Mesmay–Rieck–Sedgwick–Tancer, "the unbearable hardness of
  unknotting"); **genus-one knots: u=1 IS decidable** (Coward–Lackenby,
  "Unknotting genus one knots"); Lackenby et al. 2025 study practical
  unknotting via RL (Experimental Mathematics).
- **Montesinos:** u(K)=1 ⟹ the double branched cover Σ(K) is ±d/2-surgery
  on a knot (d = det K) ⟹ H₁(Σ) cyclic; **Lickorish linking-form
  obstruction:** the linking form must represent ±2/d; **Murasugi:**
  |σ(K)| ≤ 2u(K); Ozsváth–Szabó d-invariant refinements exist (heavy).
- **Consequence for the benchmark prompt:** any honest submission is a
  three-valued engine (True-certified / False-certified / UNKNOWN) plus an
  argument that UNKNOWN is empty on the benchmark's test distribution —
  or a new theorem. State this loudly; never present the pipeline as a
  proven decider.

## The in-repo pipeline (04-computation/unknot1_decider.py)

1. **PD parsing + faces** (rotation system) → checkerboard graph.
2. **Exact invariants:** Goeritz matrix; signature via Gordon–Litherland
   (integer congruence diagonalization, no floats); det; Alexander
   polynomial (Seifert matrix route; evaluation-interpolation with exact
   rationals).
3. **False-certificates:** |σ| ≥ 4 (Murasugi); Lickorish linking form
   ±2/d test via Smith normal form of the Goeritz matrix (convention
   fixed against the classic u(7₄)=2 example before it is allowed to
   decide); optional Alexander/Nakanishi module checks.
4. **True-certificates:** for each crossing change in the given diagram:
   det=1 screen → greedy R1/R2 reduction to the empty diagram = certified
   unknot ⟹ u=1 with an explicit certificate (crossing + move sequence).
5. **Escalation lane (not implemented day one):** bounded Reidemeister
   exploration before crossing changes (the u=1 change need not be visible
   in the input diagram — the classic completeness gap); d-invariant
   obstructions; SnapPy/Regina outsourcing where available.

## Repo-native remark

THM-2176/2191 (Gordian continuation cocycle, catalytic group length,
closure of 9_10) already live one street over: the continuation-kernel
machinery is the repo's own language for crossing-change moves. If a
session ever pushes the u=1 lane seriously, start from THM-2191's stable
envelope rather than raw diagram search, and respect MISTAKE-230–235
(no syntax-only bridges between knot moves and other carriers).
