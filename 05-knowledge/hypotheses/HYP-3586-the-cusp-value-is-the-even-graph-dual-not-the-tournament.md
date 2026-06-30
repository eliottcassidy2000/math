---
id: HYP-3586
title: OPEN QUESTION ANSWERED (does the metagraph reproduce the cusp value 4cos^2(3pi/7)?) -- NOT on the tournament side (the Z_7 circulant TOURNAMENTS of G_n give only {0.308 generic, 2.0 Paley=flat=optimal=H-max}; the binding 0.198 is SUB-tournament), but YES on the EVEN-GRAPH DUAL E_n: the doublet's autocorrelation IS the shifted 7-cycle adjacency, so 4cos^2(3pi/7) = 2 + lambda_min(A(C_7)) EXACTLY, the spectral floor of the minimal Z_7-circulant EVEN graph (the 7-cycle = the cusp of E_7); so the LRC cusp binds on the DUAL (cycle-space/even-graph) metagraph, and the tournament 3-cycle (G_n's Fiedler mode = cyclicity, THM-588, a Z_3 object) and the LRC doublet (E_n's 7-cycle, a Z_7 object) are DUAL MINIMAL CYCLES -- the apex prime 7 sets the cycle length
status: COMPUTED + verified exact (8 Z_7 circulant tournaments; 2+lambda_min(C_7)=4cos^2(3pi/7) to machine precision). The even-graph-dual placement of the binding value is the new structural result; honestly TEMPERS S30's "3-cycle <-> doublet" (they are dual, on opposite sides of the G_n/E_n duality, not the same-side mirror).
source: mac-mini-2026-06-29-S31
related:
  - HYP-3585  # mac-mini S30: the cusp landscape + the (same-side) 3-cycle<->doublet correspondence (refined here)
  - THM-588   # the metagraph Fiedler mode = cyclicity = the 3-cycle (G_n cusp's binding mode)
  - THM-586   # Paley T_7 = the H-maximizer (H=189=7*27); the optimal/flat tournament
  - THM-578   # the doublet R-tail = the binding cusp object (now placed in the even-graph dual)
  - HYP-3547  # the octonion/Fano flat cores = the Paley = the optimal (gap 2)
results:
  - 04-computation/cusp_value_vs_metagraph_apex_layer_macmini_20260629.py
  - 05-knowledge/results/cusp_value_vs_metagraph_apex_layer_macmini_20260629.out
---

# HYP-3586 -- the cusp value is the even-graph dual, not the tournament

## The open question (S30) and the honest answer
S30 asked: does the metagraph H->1 corner reproduce the LRC cusp value `4cos^2(3pi/7)` quantitatively?
Worked it with the metagraph corpus as inspiration (THM-588 Fiedler/cyclicity, THM-586 Paley, the
even-graph dual E_n). The answer is a clean split:

**NO on the tournament side.** The metagraph's apex-7 layer is the 8 `Z_7` circulant TOURNAMENTS
(size-3 connection sets, `S u -S = Z_7*`). Their autocorrelation gaps are exactly:
- `0.308` (the 6 generic circulants, `H=175=25*7`),
- `2.000` (the 2 Paley/Fano `{1,2,4},{3,5,6}`, `H=189=27*7` -- the H-MAXIMIZER, THM-586; `gap 2 =
  |Gauss sum|^2 split = (1+7)/4` = the flat/octonion-OPTIMAL, HYP-3547).
The LRC binding value `0.198 = 4cos^2(3pi/7)` is the DOUBLET gap (size 2 = THM-578 R-tail), which is
**SUB-TOURNAMENT**: `0.198 < 0.308` = the minimal tournament gap. It lies BELOW the tournament floor.

**YES on the EVEN-GRAPH dual side.** A doublet `{a,b}` depends only on its difference `d=b-a`; its
autocorrelation `(f * f)` is supported on `{0, +-d}`, i.e. **`2I + A(C_7)`** where `C_7 = Cay(Z_7,{+-d})`
is the 7-cycle. So the binding value is, EXACTLY (machine-verified):
> `4cos^2(3pi/7) = 2 + 2cos(6pi/7) = 2 + lambda_min(A(C_7))`,
the spectral floor of the shifted 7-cycle adjacency. The **7-cycle is the minimal connected `Z_7`-circulant
EVEN graph** -- the cusp of the even-graph dual metagraph `E_7`. So the LRC cusp value is an EVEN-GRAPH
spectral quantity, living in the DUAL metagraph `E_n`, not the tournament metagraph `G_n`.

## The refined correspondence: dual minimal cycles
S30 said "the metagraph 3-cycle mirrors the LRC doublet." The refinement: they are not a same-side mirror;
they are **DUAL**, sitting on opposite sides of the `G_n <-> E_n` duality (CLAUDE.md: `E_n` is the dual of
`G_n`; tournaments = cut space, even graphs = cycle space, the GF(2) split):
| | tournament metagraph `G_n` (cut space) | even-graph dual `E_n` (cycle space) |
|---|---|---|
| cusp | transitive `H=1` | empty/degenerate even graph |
| binding object | 3-cycle (cyclicity, THM-588 Fiedler mode) | **7-cycle** `C_7` (minimal `Z_7`-circulant even graph) |
| group | `Z_3` (gap `1`) | `Z_7` (gap `4cos^2(3pi/7)=2+lambda_min(C_7)`) |
| LRC role | the witness/cut side (off-path) | **the floor/BOUNDED side (the binding value)** |
Both binding objects are the **minimal CYCLE** -- the apex prime sets the length. The tournament cyclicity
(3-cycle, `Z_3`) is the cut-space minimal relation; the LRC binding (7-cycle, `Z_7`) is its cycle-space
DUAL at the apex. This is why the floor (BOUNDED = the 2nd moment = the cycle space) binds on the
even-graph side: the R-tail (THM-578) is a CYCLE-SPACE resonance, not a cut-space (score) one.

## Why this is the right placement (the metagraph inspiration paid off)
- THM-588: `G_n`'s Fiedler mode (slowest, the binding mode of the arc-flip walk) IS the cyclicity = the
  3-cycle count; mult(1)=0 (no cut/linear invariant), mult(2)=1 (the unique quadratic). So `G_n` itself
  says the binding content is the CYCLE part, not the cut part -- consistent with the value living in `E_n`.
- THM-586: the Paley `T_7` (`H=189`, `H`-max, `p|H` by the free `Z_7` action) is the OPTIMAL tournament =
  the flat core (gap 2). So `G_n`'s apex layer rehearses the OPTIMUM (Paley=flat), while `E_n`'s apex cusp
  (the 7-cycle) carries the BINDING floor. Optimum on the primal, floor on the dual.
- The even-graph dual `E_n` (first-class in this project) is exactly the cycle-space metagraph; the doublet
  = a single difference = a `Z_7`-circulant cycle = its minimal object. The binding value was always an
  even-graph eigenvalue; we had been looking for it on the wrong (tournament) side.

## What to track / next
- The floor's last bound `rho_j >= 4cos^2(3pi/7)` is now `rho_j >= 2 + lambda_min(C_7)` -- a 7-CYCLE
  spectral statement in the even-graph dual `E_n`. Does the descent's R-tail (THM-578) literally land on
  `C_7` (the minimal even graph), so the bound is the `E_7` cusp eigenvalue?
- General apex prime `p` (LRC at `n=2p`): predicts the binding value `= 2 + lambda_min(C_p) =
  2 + 2cos(2pi*floor(p/2)/p) = 4cos^2(floor(p/2)*pi/p)` -- a uniform `E_p`-cusp formula. (Check LRC6: `p=3`,
  `2+2cos(2pi/3) = 1`; matches the `Z_3` doublet gap.)
- The Paley/flat (optimum) vs 7-cycle/doublet (floor) are the two extremes of the apex `Z_7` spectrum;
  the whole floor program is "stay between the 7-cycle floor and the Paley ceiling."

## What it buys
Answers S30's open question precisely and honestly: `4cos^2(3pi/7)` is NOT a tournament-metagraph value
(those are `{0.308, 2.0}`, optimum = Paley); it is `2 + lambda_min(C_7)`, an EVEN-GRAPH spectral quantity,
the cusp eigenvalue of the dual `E_7`. The LRC floor binds on the cycle-space/dual side, at the minimal
`Z_7`-cycle -- the apex-length dual of the tournament 3-cycle. Gives a uniform `E_p`-cusp prediction for
general apex primes and relocates the binding value to its true home.
