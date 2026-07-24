---
source: kind-pasteur-2026-07-23-S135 (Opus 4.8)
status: POSITIVE residue from the S134 retraction. Worked the surviving concrete sub-task (classify the k=13
  tight instances). Results: only TWO tight instances up to dilation survive an extensive search; both attain
  the optimum at modulus EXACTLY 14 (never 14m, m>=2); and tightness has a clean covering characterization
  that lands the residual squarely in the repo's covering machinery. All exact-rational, no floats.
tags: [lrc, lonely-runner, tight-instances, classification, covering, arithmetic, residual]
related: [kps-S134 (retraction), kps-S133, THM-518, THM-515, Kravitz 1912.06034]
---

# k=13 tight instances: two of them, modulus exactly 14, and the covering reformulation

Following the S134 retraction (SGC false), the surviving concrete target was: **classify the tight instances**
(gap exactly `1/14`) — the genuine residual of any decomposition. Progress:

## 1. Only TWO tight instances survive an extensive search (up to dilation)
```
T1 = {1,2,3,4,5,6,7,8,9,10,11,12,13}         (the canonical extremiser)
T2 = {1,2,3,4,5,6,7,8,9,10,11,13,24}         (12 -> 24;  NOT a dilate of T1)
```
- BFS on the tight-instance graph (single-speed moves, replacement values ≤ 90, depth 2): **1 new node (T2), then
  nothing**.
- Extended single-replacement search from BOTH T1 and T2 with values ≤ 300: **zero new tight instances**.
So the tight set is *sparse* — nothing like a thick family in this reach. (Literature asserts an infinite family
across `n` with largest speed `~2^{n-2}`; that is one instance *per n*, consistent with `k=13` having few.)

## 2. Sharp arithmetic fact: the optimum sits at modulus EXACTLY 14
`gap = 1/14` forces the optimal time to be `τ = a/(14m)` with all speeds satisfying `v·a mod 14m ∈ [m, 13m]`.
Scanning `m = 1..8` for both tight instances:
> **Only `m = 1` ever occurs**, with `a ∈ {1,3,5,9,11,13}` = exactly the **units mod 14**. No `m ≥ 2` solution.

So a tight instance's witness time is `τ = a/14`, `gcd(a,14)=1`; the necessary speed condition collapses to
`14 ∤ v`. That condition is *weak* (most configs satisfy it) — **all the content of tightness is "no better τ."**

## 3. The covering reformulation (where the residual actually lives)
Let `D_v = {τ : ‖vτ‖ ≤ 1/14}`, 13 closed arcs each of measure `2/14 = 1/7` (total `13/7 ≈ 1.857`).
- `gap(S) ≤ 1/14  ⟺  ⋃_v D_v = [0,1)`  (the closed arcs cover the circle);
- **`S` is tight ⟺ the CLOSED arcs cover but the OPEN arcs `{‖vτ‖ < 1/14}` do not** (a boundary/tangential cover);
- **LRC(14) ⟺ 13 OPEN arcs of width `1/7` can never cover the circle.**

This is the classical covering/view-obstruction face of the problem, but stated at the exact boundary it says
precisely what the residual is: **tight instances = the tangential (measure-critical) coverings.** That places the
residual inside the repo's covering-core / singular-series machinery (THM-515/518), *not* the analytic gap-bound
route that S134 closed off.

## 4. Why this matters after the retraction
S134 killed the "buffer" premise: the Kravitz family `s/(13s+1)` leaves lossy methods no margin. What survives is
that the *hard set is tiny and arithmetically rigid*: two tight instances, both witnessed at modulus exactly 14 by
units mod 14, characterised as tangential coverings. Any proof must be **exact at these**, and there appear to be
very few to be exact about — which is encouraging for a decomposition whose residual is handled algebraically.

## 5. Open / next
- **Is {T1, T2} complete for k=13?** Search was single-replacement (≤300) + depth-2 BFS. Needs multi-speed moves
  and larger speeds before any completeness claim. **Do not claim completeness yet** (I over-generalised
  support-limited evidence once already — see S134).
- Prove the modulus-14 rigidity: *every* `k=13` tight instance is witnessed at `τ = a/14`. (Empirically true for
  both; would be a clean lemma.)
- Characterise tangential coverings by 13 arcs of width `1/7` — a finite-type combinatorial/arithmetic problem,
  and the right home for THM-518's resonance analysis.

Files: `/tmp/{tight,tight2,nbhd,isolated}.py`. Follows kps-S134 (retraction of SGC).
