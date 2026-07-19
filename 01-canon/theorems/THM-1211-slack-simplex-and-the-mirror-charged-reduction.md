---
id: THM-1211
title: THE SLACK SIMPLEX — B is six 3-simplices with |B| = 7⁻³ exactly, the four constraint slacks sum identically to 1/7, and that single identity gives a general single-run bound that collapses death-star's two-case maximiser proof into ONE inequality. (I) THE SLACK IDENTITY, which is the engine: writing the four defining slacks of death-star-S58b's bad set as f₁ = 2/7 − g_min, f₂ = 2/7 − (g_mid−g_min), f₃ = 2/7 − (g_max−g_mid), f₄ = g_max − 5/7, one has **f₁+f₂+f₃+f₄ ≡ 1/7**, identically. So each piece of B is a 3-SIMPLEX, not a box, and its incentre — all four slacks equal 1/28 — is exactly the permutation of (1/4,1/2,3/4). That is where death-star's 1/28 comes from. (II) |B| IN CLOSED FORM: in gap coordinates each piece is the corner of the cube [0,2/7]³ cut by x+y+z ≥ 5/7, i.e. a simplex with leg 1/7, so **|B| = 6·(1/7)³/6 = (1/7)³ = 1/343** exactly. This CORRECTS my own THM-1150/1151 figure |B| = 0.003367 (a coarse-grid artifact); the true concentration ratio at (1,2,3) is (2/21)/(1/343) = **98/3 = 32.67**, not 28.28. (III) A GENERAL SINGLE-RUN BOUND WITH NO CENTRE-HIT HYPOTHESIS: f₄ is continuous, decreasing, and — because the max coordinate can only ever switch to a SLOWER one — has slopes of decreasing magnitude along a run, so d_exit·L ≤ Σ_cells d_c L_c = (total drop of f₄) ≤ 1/7, giving **L ≤ (1/7)/d_exit** for every run. death-star's run = (1/28)(1/c + 1/max(a,b−a,c−b)) is the exact value AT a centre hit; this is the same bound everywhere. (IV) MIRROR-CHARGING: by death-star's mirror symmetry the mirror run r* has the same length and its exit coordinate is r's BOTTOM coordinate, so **L ≤ (1/7)/max(d_exit(r), d_exit(r*)) =: (1/7)/ρ(r)**. (V) THE COLLAPSE: summing, **sojourn ≤ (1/7)·W''** with W'' = Σ_runs 1/ρ(r), so the whole Kakeya sup bound sojourn ≤ 2/21 follows from the single inequality **W'' ≤ 2/3** — which never mentions the number of mirror pairs, and so would close death-star's 1-pair class and their ≥2-pair residue together, with no case split. VERIFIED for all 11,557 primitive directions in [1,24]³, tight at (1,2,3). (VI) AND A STRUCTURAL REASON: ρ ≥ 3P (P = number of mirror pairs) holds with 0 violations and implies W'' ≤ 2/3 immediately, being tight at P=1 (ρ=3, (1,2,3)) and P=2 (ρ=6, (5,6,23))
status: (I) PROVED — one line of algebra, 0 violations on an exhaustive lattice check. (II) PROVED — unimodular gap coordinates plus a corner-of-cube substitution, with both boundary checks strict so there is no clipping or wraparound; Kronecker sampling converges to it (+0.33%, +0.05%, −0.07% at 2·10⁵, 10⁶, 4·10⁶ points). It is also a CORRECTION to my own earlier measured |B|. (III) PROVED, including across cell boundaries, via the switching lemma (a max-coordinate switch at g_i = g_j forces d_j < d_i). (IV) PROVED given death-star's mirror symmetry; the mirror of every merged run was verified to be a merged run (0 misses). (V) the REDUCTION is proved; the closing inequality **W'' ≤ 2/3 is VERIFIED, NOT PROVED** — exhaustively for primitive d ∈ [1,24]³ only. (VI) ρ ≥ 3P is likewise VERIFIED, not proved. **Uniform r=5 remains OPEN**, and this is a sub-lemma of it, itself a rung below the r=6 sharp horn and two below the n=12 inverse theorem
source: kind-pasteur-2026-07-18-S128 (cont.83; owner: build on death-star's six-ray sojourn-max, think Kakeya needles)
depends_on:
  - THM-1154    # the six-ray/X-ray-transform reading this sharpens
  - THM-1150    # whose "six boxes, |B| = 0.003367" is corrected here to six simplices, |B| = 1/343
related: [THM-1151, THM-1152, MISTAKE-173]
script: 04-computation/simplex_slack_kakeya_kps_S128c83b.py (+ .out)
credit: death-star-S58b (HYP-7735/7736) — the bad set B, its mirror symmetry, the exact 2/21, and the centre-hit run formula whose 1/28 is explained here
---

# THM-1211 — the slack simplex, and one inequality instead of two cases

## Setting

death-star-S58b's bad set, verified there to equal `bad_from_g`, is: in the ordering
region g₍₁₎ ≤ g₍₂₎ ≤ g₍₃₎,

> g₍₁₎ ≤ 2/7,  g₍₂₎ − g₍₁₎ ≤ 2/7,  g₍₃₎ − g₍₂₎ ≤ 2/7,  g₍₃₎ ≥ 5/7.

## (I) The slack identity

Put f₁ = 2/7 − g₍₁₎, f₂ = 2/7 − (g₍₂₎−g₍₁₎), f₃ = 2/7 − (g₍₃₎−g₍₂₎), f₄ = g₍₃₎ − 5/7. Then

> **f₁ + f₂ + f₃ + f₄ = 6/7 − g₍₃₎ + g₍₃₎ − 5/7 = 1/7,  identically.**

So each piece of B is a **3-simplex** {f ≥ 0, Σf = 1/7} — not a box, as I had it in
THM-1150. Its incentre has all four slacks equal to (1/7)/4 = **1/28**, and that point is
exactly g = (1/4, 1/2, 3/4). **This is where death-star's 1/28 comes from.**

## (II) |B| = 1/343, exactly

In gap coordinates x = g₍₁₎, y = g₍₂₎−g₍₁₎, z = g₍₃₎−g₍₂₎ (unimodular, so
volume-preserving) the constraints are x,y,z ≤ 2/7 and x+y+z ≥ 5/7. Substituting
u = 2/7−x, v = 2/7−y, w = 2/7−z ≥ 0 turns this into

> **{u, v, w ≥ 0,  u + v + w ≤ 1/7}**,

the corner of the cube [0,2/7]³ — a simplex with leg 1/7, volume (1/7)³/6. The bound
u ≤ 2/7 is automatic since 1/7 < 2/7; and both boundary checks are strict
(x = 2/7−u ≥ 1/7 > 0, and g₍₃₎ = x+y+z ≤ 6/7 < 1), so there is no clipping and no
wraparound. Over six ordering regions:

> **|B| = 6 · (1/7)³/6 = (1/7)³ = 1/343 = 0.0029155.**

**This corrects my own THM-1150/1151 value |B| = 0.003367**, which was a coarse-grid
artifact — a midpoint grid has error O(surface · h) and on this set that is tens of
percent even at 180³. Kronecker sampling converges to 1/343 (+0.33%, +0.05%, −0.07% at
2·10⁵, 10⁶, 4·10⁶ points). Consequently the concentration ratio at (1,2,3) is

> (2/21)/(1/343) = 686/21 = **98/3 = 32.67**, not the 28.28 I have been quoting.

## (III) A general single-run bound — no centre hit needed

Within a cell f₄ is affine and decreasing at rate d_top. Across a cell boundary the max
coordinate can only switch to a **slower** one: a switch at g_i = g_j with both ≥ 5/7
means g_j is larger just after, so d_j < d_i. Hence along a merged run f₄ is continuous,
decreasing, with slopes of decreasing magnitude, and since every rate on the run is at
least the final rate d_exit,

> d_exit · L ≤ Σ_cells d_c L_c = (total drop of f₄) ≤ 1/7,  so **L ≤ (1/7)/d_exit.**

death-star proved run = (1/28)(1/c + 1/max(a, b−a, c−b)) **at a centre hit**. This is the
same bound for **every** run, hit or not. Verified: 0 violations, cell-runs and merged runs
alike, over all 11,557 primitive directions in [1,24]³.

## (IV) Mirror-charging

B is invariant under g ↦ 1−g and g(1−u) = 1−g(u) (death-star). So the mirror run
r* = [1−hi, 1−lo] has the **same length**, and its max coordinate is r's **minimum**
coordinate. Applying (III) to both:

> **L ≤ (1/7)/max(d_exit(r), d_exit(r*))  =:  (1/7)/ρ(r).**

(The mirror of every merged run was verified to be a merged run — 0 misses.)

## (V) The collapse: one inequality, no case split

Summing (IV) over runs:

> **sojourn(d) ≤ (1/7) · W''(d),  W''(d) := Σ_runs 1/ρ(r).**

So the entire Kakeya sup bound follows from

> **W''(d) ≤ 2/3.**

| aggregation | max over d ∈ [1,24]³ | at | verdict |
|---|---|---|---|
| W = Σ N_i/d_i, per **cell**-run | 1003/420 = 2.388 | (1,20,21) | too lossy (double-counts cells) |
| W' = Σ M_i/d_i, per merged run | 1003/420 = 2.388 | (1,20,21) | still too lossy |
| **W'' = Σ 1/ρ, mirror-charged** | **2/3** | **(1,2,3)** | **holds, and tight** |

The point is that **W'' ≤ 2/3 never mentions the number of mirror pairs.** death-star
proved the 2/21 ceiling on the 1-mirror-pair class and left the ≥2-pair class as a
*verified* residue (max 4/105, gap 2.5×). Under this reduction there is nothing to split:
one inequality covers both. For the record, W'' on the ≥2-pair class maxes at 4/7 = 0.571,
comfortably under 2/3.

## (VI) Why W'' ≤ 2/3 should be true

With P = number of mirror pairs, **ρ ≥ 3P** holds with 0 violations, and it implies the
target immediately: Σ_pairs 1/ρ ≤ P/(3P) = 1/3, so W'' = 2·Σ_pairs 1/ρ ≤ 2/3.

| P | min ρ observed | witness | 3P |
|---|---|---|---|
| 1 | 3 | (1,2,3) | 3 |
| 2 | 6 | (5,6,23) | 6 |
| 3 | 14 | (1,13,14) | 9 |
| 4 | 17 | (1,16,17) | 12 |
| 5 | 21 | (1,20,21) | 15 |

Tight at P = 1 and P = 2; slack above. The P = 1 case, ρ ≥ 3, is exactly death-star's
c ≥ 3 — which they get from a, b−a, c−b being positive integers summing to c, but only
**at a centre hit**. Off a hit that argument is unavailable, and ρ ≥ 3 is open.

## Honest status

(I)–(IV) are proved. **(V) is a reduction, not a closure: W'' ≤ 2/3 is verified for
primitive d ∈ [1,24]³ and not proved**, and (VI) is likewise verified only. What the
session buys is that the residue is now *one* inequality about a purely combinatorial
weight, rather than a proved class plus an enumerated class. **Uniform r=5 remains open**,
and it is itself a rung below the r=6 sharp horn and two below the n=12 inverse theorem —
the fleet's actual wall.

## Named next
- **Prove ρ ≥ 3.** This is the P = 1 base case and would re-derive death-star's result by
  a route that does not assume a centre hit. Concretely: rule out a run whose top and
  bottom coordinates both have rate ≤ 2.
- **Prove ρ ≥ 3P**, or directly Σ_pairs 1/ρ ≤ 1/3. The mechanism to look for is why more
  mirror pairs force larger exit rates — more runs means more crossings of g = 5/7, which
  is a counting constraint on the d_i.
- Push the verification of W'' ≤ 2/3 past [1,24]³; the check is cheap and linear in N.
