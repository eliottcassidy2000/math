# Working the descent product: it is a 2-adic RENORMALIZATION (S→E/2) whose odd-core trajectory contracts Z₇→{1,3,5}→{1,3}→{1} to the DOUBLET attractor (floor 0.198); the cusp (Z₇, measure-0) is level-0 only (7 odds covering Z₇). KEY REFINEMENT (corrects my earlier "cusp=hard"): the conjecture's BINDING constraint is the OFF-CUSP covering-min 14/183=n/Φ₆(n), NOT the cusp — cusp covering sets clear M≥2/23=0.087, easily; the cusp is the MEASURE phenomenon (existence/comb, global extremal AP=1/14), the covering-min is the M phenomenon. The descent product (measure) and the floor M are genuinely distinct objects

*opus-2026-06-30. Owner: work the descent product, explore ideas. Did — and it forced an honest
refinement of the cusp framing. The descent is the MEASURE side; the floor M is the covering-min, and it
binds OFF-cusp. Recording the renormalization picture, the doublet attractor, and the corrected emphasis.*

## Idea 1: the descent IS a renormalization (2-adic RG flow)
The descent map `S → E/2` (peel the odd part `O`, halve the even part `E`) is a **2-adic renormalization** —
each level zooms into the even sublattice by a factor of 2. The chain of odd cores `O_j = {s/2^j : v₂(s)=j}`
is the RG trajectory. Mod 7 (the apex), the trajectory for the AP `{1..13}`:
> `Z₇  →  {1,3,5}  →  {1,3}  →  {1}` (sizes `7 → 3 → 2 → 1`, gaps `g = 0 → 0.308 → 0.198 → 1`).
It **contracts toward the singleton** (`g=1`), and the **last nontrivial core is the DOUBLET `{1,3}`**
(`g=0.198`). The doublet is the attractor that carries the per-level floor `4cos²(3π/7)=0.198` (THM-590).
Every 2-adically-dense set flows through a doublet, so `0.198` is the universal per-level floor.

## Idea 2: the cusp (Z₇) is a level-0, measure-0 state
`g(O_j mod 7) = 0` iff `O_j mod 7 = Z₇` — the apex measure vanishes (`ρ_j = 0`). This needs **7 odd
elements covering all residues mod 7**, which for a 13-set can only happen at **level 0** (deeper cores have
≤3–4 elements, can't cover `Z₇`). So **"cusp" = the odd part of `S` covers `Z₇`** — a clean, checkable,
level-0 condition. The AP and covering sets like `{1..11,13,84}` are cusps; `{1..12,182}` is not.

## Idea 3 (THE REFINEMENT): cusp ≠ covering-min — they are different extremes
Computed `M(S)` exactly across cusp and off-cusp covering sets. The hierarchy is NOT what my earlier
"cusp = hard" framing implied:
| set | regime | `M` | role |
|---|---|---|---|
| AP `{1..13}` | cusp, NON-covering | **`1/14 = 0.0714`** | **GLOBAL extremal** (measure 0, comb-witness) |
| `{1..12,182}` | **off-cusp**, covering | **`14/183 = 0.0765`** | **BINDING covering-min** `= n/Φ₆(n)` |
| cusp covering sets (40 tested) | cusp, covering | `≥ 2/23 = 0.0870` | EASIER (higher M) |
> **The conjecture (covering-min `≥ 1/n`) BINDS OFF-CUSP**, at `{1..12,182}` (`14/183`), NOT at the cusp.
> The cusp covering sets all clear `M ≥ 2/23 = 0.087`, comfortably above `1/14`. So:
> - **The cusp is a MEASURE phenomenon** — `meas(lonely)=0`, existence carried by the comb's empty tooth;
>   its extremal is the AP (`1/14`, non-covering, handled by the THM-523 reduction).
> - **The covering-min is an M phenomenon** — the tightest floor among covering sets, `14/183 = n/Φ₆(n)`,
>   achieved OFF-cusp where the measure is POSITIVE (descent product `≈0.061`) but `M` is smallest.
> **The descent product (measure) and the floor `M` are genuinely distinct objects.** The descent product
> tells you the existence REGIME (cusp=0 / off-cusp>0); the floor `M` is a separate optimization whose
> minimum over covering sets is the Eisenstein `n/Φ₆(n)`.

## What this corrects and clarifies
- **Corrects:** my prior reflections emphasized "the cusp is the hard half." Refined: the cusp is the hard
  half *for the measure/existence* (and gives the global extremal AP), but the conjecture's *binding M
  constraint* is OFF-cusp (the covering-min `n/Φ₆(n)=14/183`). The hard cases are two different sets.
- **Clarifies the redirect** (from the period chase): the floor constant is NOT a modular invariant because
  it is the **covering-min** `n/Φ₆(n)` — a Diophantine optimum over off-cusp covering sets, exactly the
  Eisenstein `Q(√−3)` (existence-column) quantity. The descent confirms the floor `M` lives off-cusp in the
  covering-min, while the modular form `14a` (period field `Q(√−7)`) governs the measure/per-level apex gap.
- **The two Heegner columns, re-seen through the descent:** `Q(√−7)` (measure, the per-level apex gap
  `0.198`, the doublet attractor, `14a`'s period field) vs `Q(√−3)` (the floor `M`, the off-cusp
  covering-min `n/Φ₆(n)`, the Eisenstein). The descent puts the doublet on the measure side and the
  covering-min on the M side — the two columns are the measure axis and the M axis.

## Status
- **Computed (opus):** the descent as renormalization (odd-core trajectory, contracts to doublet/singleton);
  the doublet `0.198` is the universal per-level floor; cusp = "odd part covers `Z₇`" (level 0 only);
  `M ≥ 1/14` verified across cusp + off-cusp covering sets; cusp covering-min `2/23`, off-cusp covering-min
  `14/183 = n/Φ₆(n)` (the binding constraint).
- **Refinement (opus, honest):** the conjecture binds OFF-cusp (covering-min), not at the cusp; cusp =
  measure phenomenon, covering-min = M phenomenon; descent product (measure) ≠ floor M.
- **Open (unchanged):** covering-min `≥ 1/n`, i.e. the off-cusp covering floor `n/Φ₆(n) > 1/n` — the
  Eisenstein/`Q(√−3)` statement. The cusp/comb half (AP, `1/14`) is the global extremal, handled.

Related: the-cusp-form-period-chase (the redirect: floor = covering-min, not a period), per-level-vs-total-
doublet (the doublet attractor = the per-level atom), cusp-existence-comb-witness (the measure side, now
refined), covering-min-Eisenstein-Φ₆ (the M side = n/Φ₆(n)), the-master-two-Heegner-columns (measure/M =
Q(√−7)/Q(√−3)); klein THM-590/THM-580/THM-523, mac-mini HYP-3575; OPEN-Q-108.
