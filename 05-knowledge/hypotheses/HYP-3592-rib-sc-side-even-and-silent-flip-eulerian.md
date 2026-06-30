---
id: HYP-3592
title: COMPLEMENT to klein-S12/codex-S675b on the even-graph & Eulerian structure -- THREE distinctive results (1) the equinumerosity is LABELED-only (tournaments = 2^(n-1) x even graphs via Cut(+)Cycle; the iso fan-out A000568/A002854 = 1,1.7,3.5 GROWS, the per-fiber mechanism behind klein's Royle sandwich; clean decomposition only at ODD n, Cut∩Cycle=dim n-2 at even n); (2) NEW PARITY THEOREM -- the d=1 RIBS (SC-NS section) have EVEN degree on the SC side at ALL n (the complement involution R fixes each SC class and pairs its NS-neighbors into complement-pairs), so the d=1 ribs are EULERIAN at n=4,5 (NS side also even) and the NS side first breaks parity at n=6 -- a DIFFERENT Eulerian section from codex's d=m black/blue; (3) the prior Eulerian thread is THM-413 Lemma A (a SILENT FLIP = a metagraph SELF-LOOP exists iff the value-multigraph G_x is Eulerian; prime-3, x=(2n-1)/3; n=14: C=27=3^3, x=9), so the metagraph's own degree parity (C(n,2) - #selfloops) is governed by the silent-flip/Eulerian mechanism. Verified: full arc-flip metagraph is NOT Eulerian (odd degrees n=4,5,6); d=1 edge types at n=5 = 16 SC-NS + 12 SC-SC + 2 NS-NS = 30 = E(G_5)
status: COMPUTED (n=4,5,6 exhaustive iso metagraph; SC-side-even verified n<=6 + proof sketch via R; fan-out + labeled factor exact). COMPLEMENTS klein-S12 (HYP-3591, the sandwich) + codex-S675b (black/blue d=m Eulerian sections); the d=1 rib SC-side-even theorem + the silent-flip/self-loop connection are the new pieces.
source: mac-mini-2026-06-29-S32
related:
  - HYP-3591  # klein-S12: the Royle sandwich (Eulerian<=tournaments<=all); the three Burnside counts; surplus=obstruction
  - merged-line-parity-even-odd-s675b  # codex: black=Eulerian (cycle-space), blue=odd boundary (d=m lift)
  - HYP-3590  # mac-mini S31: the floor binds on the even-graph dual E_n (the cusp = 7-cycle)
  - THM-584   # complement R = the antipodal involution; SC = R-fixed classes (the rib theorem's engine)
  - THM-578   # the doublet R-tail = the SC-NS/rib object (obligation D)
results:
  - 04-computation/equinumerosity_and_eulerian_metagraph_macmini_20260629.py
  - 05-knowledge/results/equinumerosity_and_eulerian_metagraph_macmini_20260629.out
---

# HYP-3592 -- the rib SC-side-even theorem + silent-flip Eulerian (complement to klein-S12/codex-S675b)

This session worked the owner's two asks (tournament/even-graph equinumerosity; "sections of the merged
metagraph as Eulerian"). klein-S12 (HYP-3591) and codex-S675b independently nailed the SANDWICH framing
and the d=m black/blue Eulerian sections concurrently; below are the three pieces NOT in their accounts.

## (1) Equinumerosity is LABELED-only -- the per-fiber mechanism behind the sandwich
- **LABELED**: tournaments `= 2^(n-1) x` even graphs. Verified ratio `4,8,16,32` (n=3..6) `= 2^(n-1)` = the
  **CUT/score factor**: the GF(2) split `tournament = CUT (dim n-1, scores) (+) CYCLE (dim C(n-1,2), even
  graph)`. Each even graph is the cycle-shadow of exactly `2^(n-1)` labeled tournaments (the scores).
- **ISO**: `A000568 (2,4,12,56)` vs `A002854 (2,3,7,16)` diverge -- this is klein's sandwich `Eulerian <=
  tournaments`. The MECHANISM: the iso **fan-out** `A000568/A002854 = 1, 1.7, 3.5` (n=4,5,6) GROWS -- each
  even-graph class hosts an increasing number of tournament classes, because the cut (score) adds genuine
  iso TYPES (not just "copies" -- this REFUTES the old thermodynamic intuition in
  `tournaments-and-even-graphs.md`). The cut & cycle do not decouple under `S_n` (different stabilizers).
- **PARTICULARITY**: the cycle-shadow `T_cycle=(I+L)T mod 2` is canonical only at **ODD n** (`Cut∩Cycle=0`);
  at **EVEN n**, `dim(Cut∩Cycle)=n-2` so the even-graph projection is ambiguous (the floor case n=14 is even
  -- the apex `Z_7` shadow lives one descent level down, where the odd core is reached).

## (2) NEW PARITY THEOREM: the d=1 ribs are SC-side-even (and Eulerian at n=4,5)
The d=1 (wiggly, single-arc-flip) iso-class metagraph splits by SC/NS type. At n=5: `30 = 16 SC-NS (ribs) +
12 SC-SC (spine) + 2 NS-NS (sea)` (matches `E(G_5)=30`). The **full metagraph is NOT Eulerian** (odd-degree
classes at n=4,5,6). But the RIBS have a clean parity:
> **THEOREM (verified n<=6; proof via R): every SC class has EVEN SC-NS (rib) degree.**
> *Proof.* `C` self-complementary means the complement involution `R` fixes `C`, represented by a relabeling
> `sigma` on `C`'s arcs. Flipping arc `e` lands in class `C'`; applying `R`, flipping arc `sigma^{-1}(e)`
> lands in `R(C')=comp(C')`. If `C'` is NS (`C' != comp(C')`) then `e != sigma^{-1}(e)`, so the rib-arcs
> pair up `e <-> sigma^{-1}(e)`, giving the complement-pair of NS-neighbors `{C', comp(C')}` -- TWO distinct
> classes. Hence the SC-NS degree of `C` is even. (Fixed arcs of `sigma` give SC-neighbors = spine, not
> ribs.) ∎
- VERIFIED: SC-side rib degrees all-even at n=4 `{2}`, n=5 `{2,4}`, n=6 `{2,4,8,10}`. The **NS side** is also
  even at n=4 `{2}`, n=5 `{2,6}` (so the ribs are **fully EULERIAN at n=4,5**), but first breaks at n=6
  `{2,3,4,6}` (odd `3`). So the bipartite SC-NS ribs are "SC-side-Eulerian" universally, fully Eulerian only
  for small n. This is a `d=1` section result, DISTINCT from codex-S675b's `d=m` black(Eulerian)/blue(odd).
- This ties to the LRC: the ribs ARE the SC-NS edges = the R-tail/doublet object (THM-578, obligation D);
  the SC-side-even parity is a clean certificate on the obligation-D section.

## (3) The prior Eulerian thread: silent flips = self-loops (THM-413 Lemma A)
The merged metagraph's degree parity is `C(n,2) - #(self-loops)`, and a **self-loop = a SILENT FLIP**
(class-preserving arc flip). THM-413 Lemma A: **a silent flip of runner `x` exists iff the value-multigraph
`G_x` is Eulerian** (a silent flip = a balanced edge 2-coloring, exists iff all degrees even); the
odd-degree defects of `G_x` are exactly `{x, rho(2x)}`, merging iff `x=rho(2x) iff 3x≡0 iff x=(2n-1)/3` --
the **order-3 torsion** of `Z/(2n-1)`. For n=14: `C=27=3^3`, silent runner `x=9`. So the project's existing
"Eulerian" content is the silent-flip/self-loop layer (prime-3), and the metagraph's (non-)Eulerian-ness is
governed by this mechanism -- the self-loop count whose parity decides each class's degree parity.

## Unification (with klein-S12/codex)
Both asks point to the **CYCLE SPACE / EVEN GRAPHS**: equinumerosity = `Cut (+) CYCLE`; Eulerian = even
graph = cycle space (codex black); silent flip <=> `G_x` even (Eulerian); the even-graph DUAL `E_n` is the
metagraph OF Eulerian graphs (where the LRC floor binds, HYP-3590). The complement involution `R` (= the
2-adic involution, the cut/cycle reference) is the common engine: it fixes SC (ribs SC-side-even), splits
even/odd boundary (codex black/blue), and its `-1`-eigenspace is the odd/cusp obstruction (klein-S12). The
cycle space governs the count, the Eulerian structure, and the floor at once.
