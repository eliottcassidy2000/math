---
id: HYP-3566
title: TRANSITIVITY REFRAME of the LRC floor -- ρ_j≥c (the open piece) is the positivity of a SET-INDEPENDENT Z_7-cyclic (cyclotomic) reference-collapse on the 2-adic-descended apex-7 core; the floor fails through CV(N_R)² only because Z_14 is non-transitive (m_R→0 trap), and the descent's role is to peel the non-transitive 2-part and expose the transitive Z_7 where the second moment becomes set-independent (the S75e Fejer-Bochner SOS = that collapse)
status: REFRAME / testable proposal (klein-2026-06-29-S7b). Grounded in proven pieces (CV(H) S_n-collapse THM-589; CV(N_R) set-dependence HYP-3554; the descent THM-580; apex-7 Z_7/Heegner HYP-3547). The Z_7-invariance of the descended core is the concrete, falsifiable test.
source: klein-2026-06-29-S7b
depends_on:
  - THM-580   # the 2-adic descent (peels even speeds -> odd apex-7 core)
  - THM-589   # CV(H) S_n-transitive collapse, ~2/n (the rehearsal); even-overlap = doubling-2 parity
  - HYP-3554  # CV(N_R)² set-dependent/unbounded (Z_14 non-transitive) = why the variance route is a trap
related:
  - HYP-3548  # the two lines; ρ_j≥c is the open one
  - HYP-3550  # floor = positive Euler product (set-independent BRIDGE)
  - HYP-3553  # Gamma_0(14) congruence-density route (the transitive collapse)
  - HYP-3547  # apex-7 = Z_7/Paley/Fano/octonion (the transitive symmetry, correctly located)
  - HYP-3535  # the S75e cyclotomic Fejer-Bochner SOS (= the Z_7 reference-collapse)
---

# HYP-3566 — the floor is a transitivity problem; ρ_j≥c is the Z_7 cyclotomic reference-collapse

## Claim (reframe)

The remaining LRC floor target `ρ_j ≥ c` (per-level decorrelation on the 2-adic descent, HYP-3548's only
thin line) should be proved as a **set-independent transitive collapse**, not a per-set variance bound:

1. `CV(N_R)²` (THM-579 gatekeeper) is set-dependent and unbounded (HYP-3554) because the `Z_14` sheet
   action is **not transitive** over the speed structure (`m_R → 0`). The variance route is a TRAP.
2. `CV(H)²` (THM-589) is set-independent and `→ 0` because `S_n` is **transitive**. The transitive
   collapse to a single reference object is the mechanism that works (the finite rehearsal, HYP-3560).
3. `14 = 2·7`: the 2-adic descent (THM-580) **peels the non-transitive 2-part**, exposing the residual
   **apex-7 core**, where `Z_7` is cyclically **transitive** (Paley/Heegner `h=1`, HYP-3547). On that core
   the second moment is a single `Z_7`-average; its positive-definite certificate in the `Z_7`-Fourier
   (cyclotomic) basis is exactly the **S75e Fejér–Bochner SOS** (HYP-3535). So `ρ_j ≥ c` = positivity of
   the `Z_7` cyclotomic Gram form — set-independent because `Z_7` is transitive.

This relocates the octonion/Singer-`Z_7` structure from the (dead) `b_1^-` route (HYP-3563) to its proven
home (the descent's odd core), where it is the floor's transitive vehicle.

## Concrete test (for the floor owners)

At a binding descent level, after the 2-adic peel, **check whether the residual second moment is
`Z_7`-cyclic invariant** (invariant under cyclic relabeling of the odd-core speeds mod 7).
- If YES: `ρ_j ≥ c` is the positivity of its cyclotomic Gram form — a finite, set-independent SOS (the
  S75e certificate). The floor closes via the transitive collapse.
- If NO (not `Z_7`-symmetric at some binding level): the transitivity must be manufactured by a larger
  congruence group (`Γ₀(14)`, HYP-3553), and that level localizes precisely where the bespoke `Z_7`
  symmetry breaks — a sharpened target.

Either outcome advances the proof: it converts "bound a per-set variance uniformly" (a trap) into "verify a
transitive symmetry of the descended core" (a finite, structural check).

## Why this is more than analogy

The same doubling-`2` parity underlies both sides: THM-589's even-overlap cancellation (odd overlap runs
vanish, the witness concentrates) and THM-580's even-speed peel (the floor descends). The witness side is
the proven instance of the transitive collapse; the floor side is the same mechanism awaiting the transitive
group, which the descent exposes as `Z_7`. The persistence-test discipline (HYP-3563) keeps this honest: use
the proven apex-7 `Z_7` bridge, not dimensional `7`-coincidences.
