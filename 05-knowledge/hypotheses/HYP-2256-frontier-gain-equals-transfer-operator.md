---
id: HYP-2256
title: Frontier-gain = transfer operator = incremental partition function — and the residue-profile proof of range-stability
status: PARTLY PROVED (residue-profile reduction is rigorous); the unification is a frame
source: claudebox-2026-06-03-S629
related:
  - HYP-2245  # partition functions everywhere (this is its algorithmic shadow)
  - THM-411   # tight-family finiteness / range-stability (now EXPLAINED)
  - THM-415   # quantitative C'(n) (min M)
  - HYP-2230  # unit-distance bridge
---

# HYP-2256 — frontier-gain tables are transfer operators

**Principle (from the user's unit-distance optimization).** A global objective evaluated over a
search that builds configurations element-by-element should be made **state-local**: carry the value
and update it by the **marginal gain** of the added element, read from a **frontier-gain table** keyed
by the boundary state. This is exactly a **transfer operator** — the incremental form of a partition
function. So this is the *algorithmic shadow* of "partition functions everywhere" ([[HYP-2245]]):
**every partition function we have identified yields a state-local incremental algorithm**, and every
such speedup is a transfer-operator factorization.

## Three instances (measured, `frontier_gain_s629.out`)

1. **Unit distance (the user's example).** On the triangular (= Eisenstein CM) lattice the
   frontier-gain table is the 6-neighbor rule: adding a point `p`, gain `= #{neighbors of p already
   placed}` (O(1)), vs a full distance count (`C(k,2)` checks). At the **n=22 frontier the
   per-extension reduction is ≈ C(22,2) = 231×** (≈ the reported 211×); averaged over the beam's
   growth, ×83. Best edge-count identical (49 on this lattice beam).
2. **LRC.** The frontier is the **shell profile**: the running min residue-distance `md[(m,a)]` over
   the `(ℤ/m)*` witnesses (`m ≤ 2n-1`). Adding a speed `v` updates each `md` by one `min` — O(#witnesses)
   — and `gap_shells = max md/m`. In a DFS enumeration, clone the parent profile and add one speed
   (O(#witnesses)) instead of recomputing (O(#witnesses·depth)): a depth-fold win, on top of the ×56
   `gap_shells` win (S628). This is the transfer recursion `Z_{S∪{v}} = Z_S + (z-1)·local` (S626) in
   operational form.
3. **Tournaments (noted).** Vertex insertion `n→n+1` (canon "Mode A") is a frontier recursion: scores
   / `H` / the depth update by the new vertex's arc-gains — the metagraph's transfer operator.

## The breakthrough: the residue-profile reduction (rigorous)

`gap_shells(S)` reads each speed `v` only through `v mod m` for `m ≤ 2n-1`, hence (CRT) only through
**`v mod L`, `L = lcm(2,…,2n-1)`**. Verified: `[1,3,4,5]` and `[2521,3,4,5]` (≡ mod `L=2520` at n=5)
have identical `gap_shells`. Consequences:

> **The frontier of the whole LRC tight/min-M problem is the residue profile (the multiset of speeds
> mod `L`), of which there are FINITELY many — independent of the box size `R`.**

- This **proves the range-stability** observed empirically in S621 (tight count constant as `R` grows,
  up to `8n`): the tight property is a function of the residue profile, and raising `R` adds no new
  profiles. (Modulo the high-shell residual where `gap_shells < M`; on the shell-decided part it is a
  theorem.)
- It gives a **box-free algorithm** for `min M` (THM-415) and the tight count: enumerate residue
  profiles, not raw configs — collapsing an unbounded search to a finite one.
- It is the LRC echo of the disproof's "bounded root discriminant": the relevant arithmetic lives at
  a bounded level (`mod L`), so the problem is uniformly finite there.

## To do
1. Implement the residue-profile enumerator for `min M` / tight count and confirm it reproduces S621
   (and pushes to larger `n` cheaply).
2. Bound `|profiles|` and turn "range-stable ⟹ finite" into the clean finiteness proof for THM-411
   (currently VERIFIED-only) on the shell-decided part.
3. Frontier-gain DP over residue profiles as a route to the uniform quantitative-`C'` bound
   ([[HYP-2240]]): a transfer matrix on `(ℤ/L)`-profiles whose value is `min M`.
