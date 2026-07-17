# THM-941 — The block-split (GLevel) lift: the generic engine and the quad core (death-star-2026-07-16-S35)

**Status:** PROVED (Lean, kernel-pure — `TournamentH7/LRCBlockSplitLift.lean`,
axioms propext/Classical.choice/Quot.sound; verify the build report in the session
log). Lifts THM-938 through kps-S22's fat-block window step. Source: HYP-7172.

## Statement

1. `lonely_of_block_split` (**the generic composition**): cite k ≤ 12 runners ≤ B at
   gap 1/(k+1); ONE fat block `ws` whose kps-S22 mass fee fits the transported
   window L = 2δ at width ε
       Σ_{w∈ws} [L/7 + 3/(7w) + ε(wL + 3)] < L − ε;
   an S19 singles `ChainOK` tail INSIDE the ε-window; every runner cited/block/tail
   ⟹ lonely. The engines compose because `block_window_step` certifies the WHOLE
   [x, x+ε] window and `good_chain` finds its point inside it. The block size is
   implicit in the fee (union density kills it at 7 — the 7-wall).
2. `lonely_or_quadCore` (**the lifted dichotomy**): on THM-938's TripleDenseCore,
   the dense TRIPLE {w(j), w(j+1), w(j+2)} enters the block with tail width
   ε = 3/(2·w(j+3)) (B = w(j−1), or 1 at j = 0): fee passes ⟹ lonely; else
   **QuadDenseCore** = TripleDenseCore + (explicit triple-fee failure ∨ the deferred
   top-corner j ≥ 10).
3. `residualObligation_of_quadCore` / `lrc14_of_quadCore` — the wire, strictly
   sharper than THM-938's TripleCoreObligation.

## Honest measures

The position-free fat-mass ledger is ≈ 4.5× costlier than the window-aware S19 step
(a G-single enters at w·L > 34/5 vs 3/2), so few random dense cores pass the triple
fee (referee: 5 closures in a 108k dense-leaning mix — structural value, the fee is
now formal). The j ≥ 10 top-corner is deferred into the core (one session, one
family shape). The dichotomy ladder now reads: singles (THM-937) → pair (THM-938) →
triple block (this) → the 7-wall.

## Referee

Same script as THM-940: lifted-dichotomy exhaustiveness PASS on 108k families
(singles | pair | triple-block | quad-core partition, no leaks).
