---
id: HYP-3585
title: CUSP EXPLORATION -- the apex-7 cusp of X_0(14) is the finite Z_7-core gap landscape (klein-S9's corrected frame: the MINORANT/min not the average); it is BIMODAL, complement-symmetric, with only 5 distinct cyclotomic values {0, 0.198=4cos^2(3pi/7), 0.308, 1, 2}, the DOUBLETS (and their 5-complements) binding at the floor 0.198 (=THM-578 R-tail) and gap=0 ONLY at O=Z_7 (the disproof boundary); the Fano/perfect-difference-set cores ({1,2,4},{3,5,6},...) are the flat/octonion-optimal (gap 2); the Z_7^*-averaging overshoots for 30/126 cores (Jensen, confirming klein-S9's correction of my S28). AND the MISSING CUSP REHEARSAL (S29) is FOUND: the metagraph H->1 transitive-limit corner binds at the 3-CYCLE (H=3, one cyclic triangle) = the exact MIRROR of the LRC doublet -- both the minimal non-trivial RELATION (doublet=minimal resonance pair, 3-cycle=minimal cyclicity), tying the LRC R-tail (THM-578) to the metagraph cyclicity (THM-588)
status: COMPUTED (exhaustive 128-core Z_7 landscape; metagraph H-spectra n=5,6). Confirms klein-S9 (the min, the overshoot). The cusp rehearsal (metagraph H->1 <-> LRC m_R->0, both binding at the minimal relation) is a verified structural correspondence; not a proof of the floor.
source: mac-mini-2026-06-29-S30
related:
  - HYP-3581  # klein-S9: rho_j>=4cos^2(3pi/7), the finite min/minorant; averaging invalid (corrects my S28)
  - HYP-3580  # mac-mini S29: the proof lives at the cusps; the MISSING cusp rehearsal (found here)
  - HYP-3575  # mac-mini S27: the Z_7 cyclotomic Gram (the right mechanism, per klein-S9)
  - THM-578   # the DOUBLET R-tail = the binding cusp core
  - THM-588   # the 3-cycle/cyclicity = the metagraph cusp's binding object
  - HYP-3547  # the octonion/Fano (the flat/optimal cores)
results:
  - 04-computation/cusp_landscape_exploration_macmini_20260629.py
  - 05-knowledge/results/cusp_landscape_exploration_macmini_20260629.out
---

# HYP-3585 -- the cusp landscape, and the cusp rehearsal found

## The corrected frame (klein-S9)
My S28 proposed `Z_7^*`-AVERAGING to close `rho_j>=c`; klein-S9 (HYP-3581) corrected it: `rho_j>=c` needs
the **MINIMUM** over cores (the Fejer-Bochner MINORANT, S27), not the average (Jensen gives the mean, which
OVERSHOOTS). The cores are SUBSETS of `Z_7` (finite, 128 of them), so the set-independent floor is a
DIRECT FINITE MINIMUM. I mapped the whole landscape.

## The cusp is a SHORT, bimodal, complement-symmetric cyclotomic landscape
The gap `min_{k!=0}|sum_{x in O} w^{kx}|^2` over all 128 cores takes only **5 distinct values**:
`{0, 0.19806, 0.30798, 1.0, 2.0}`, all in `Q(cos 2pi/7)`. By core size:
| size | gap(s) |
|---|---|
| 1, 6 | `1.0` |
| **2, 5** | **`0.19806 = 4 cos^2(3pi/7) = 2+2cos(6pi/7)`** (the FLOOR; binds) |
| 3, 4 | `0.30798` (generic) or `2.0` (Fano/perfect-difference-set) |
| 7 | `0.0` (the full `Z_7` = the DISPROOF boundary, off the floor) |
- **The floor `0.198` binds at the 21 DOUBLETS and their 21 5-complements** (gap is complement-symmetric,
  `gap(O)=gap(Z_7\O)`, since the full character sum vanishes). The doublet is exactly THM-578's R-tail
  object -- so the binding cusp = obligation D.
- **`gap=0` ONLY at `O=Z_7`** (the full covering) -- the disproof boundary, which is OFF the floor (a
  covering core is not a lonely core). So the floor is `> 0` everywhere it matters, bottoming at `0.198`.
- The **Fano/perfect-difference-set cores** (`{1,2,4},{3,5,6}` and all translates, gap `2.0`) are the
  flat/octonion-optimal (HYP-3547/3575) -- the BEST cores, the opposite extreme from the doublets.

## The missing cusp rehearsal (S29) is FOUND -- it is the metagraph H->1 corner
S29 asked for a finite MODEL of the cusp (the metagraph models the bulk; the cusp = `m_R->0` had no
rehearsal). It is the metagraph's **transitive-limit corner** (`H->1`). VERIFIED (n=5,6): the H-spectrum
bottoms at the transitive `H=1` (the cusp) and its smallest neighbor is `H=3` -- a single 3-CYCLE. So:
> the metagraph cusp (`H->1`) binds at the **3-cycle** (one cyclic triangle, `H=3`), the EXACT MIRROR of
> the LRC cusp (`m_R->0`) binding at the **DOUBLET** (2-residue core, `rho_j=0.198`).
Both binding objects are the **minimal non-trivial RELATION**: the 3-cycle is minimal cyclicity (THM-588's
unique quadratic invariant), the doublet is minimal resonance pair (THM-578's R-tail). So the cusp
rehearsal exists and is exact: `transitive (H=1) + 3-cycle (H=3)` on the metagraph mirrors `lonely-full +
doublet` on the LRC, and the binding is the minimal relation on both sides. This is the cusp-side analog of
the bulk rehearsal (CV(H)~2/n <-> rho_j), now at the boundary.

## The correspondence to track (the cusp's binding object, four faces)
> **doublet (LRC core, rho_j=4cos^2(3pi/7)) <-> 3-cycle (metagraph, H=3) <-> R-tail (THM-578, obligation D)
> <-> cyclicity (THM-588, the unique quadratic).**
The minimal relation IS the binding object of the cusp, in all four registers. This is the thread that
unifies the cusp -- and it says the floor's last bound is a DOUBLET/3-cycle statement: `4cos^2(3pi/7)` is
the cyclotomic value the R-tail (THM-578) and the cyclicity (THM-588) both reduce to at the cusp.

## Abnormalities / things to track
- **Bimodal gap** (only 5 values; size-3 cores are EITHER 0.308 OR 2.0 -- generic vs Fano) -- a sharp
  dichotomy, no intermediate.
- **Averaging overshoot 30/126** (Jensen) -- the precise count of cores where the (invalid) average exceeds
  the (valid) min; klein-S9's correction made exact.
- **Complement symmetry** `gap(O)=gap(Z_7\O)` -- the doublet binds WITH its 5-complement (a 5-speed core);
  worth checking whether the binding LRC config is a doublet OR a 5-set.
- **The disproof boundary is a single point** (`O=Z_7`, gap 0) -- the floor degenerates at exactly one
  core, the full covering, which is structurally excluded.

## What it buys
A complete map of the apex cusp (finite, 5 cyclotomic values, doublet-binding, klein-S9-corrected), the
FOUND cusp rehearsal (the metagraph H->1 corner, binding at the 3-cycle = the doublet's mirror), and the
unifying thread (the minimal relation = doublet = 3-cycle = R-tail = cyclicity). The floor's final bound is
a doublet statement, `rho_j >= 4cos^2(3pi/7)`, and the metagraph H->1 corner is its finite rehearsal.
