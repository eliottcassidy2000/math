---
id: THM-1211
title: THE SLACK SIMPLEX — B is six 3-simplices, exactly ⅙ of the bounding box death-star uses, with |B| = 7⁻³ = 1/343 exactly; the identity behind it explains their 1/28, and their residual counting inequality is restated as D ≥ 3P. (I) THE SLACK IDENTITY: writing the four defining slacks of death-star-S58b's bad set as f₁ = 2/7 − g_min, f₂ = 2/7 − (g_mid−g_min), f₃ = 2/7 − (g_max−g_mid), f₄ = g_max − 5/7, one has **f₁+f₂+f₃+f₄ ≡ 1/7**, identically. A set cut out by four nonnegative affine functions of constant sum is a SIMPLEX, never a box. Its incentre — all four slacks equal 1/28 — is exactly the permutation of (1/4,1/2,3/4), which is where death-star's 1/28 comes from. (II) |B| IN CLOSED FORM: in gap coordinates each piece is the corner of the cube [0,2/7]³ cut by x+y+z ≥ 5/7, a simplex of leg 1/7, so **|B| = 6·(1/7)³/6 = (1/7)³ = 1/343** exactly. death-star-S58c's run bound uses B's BOUNDING BOX [1/7,2/7]×[3/7,4/7]×[5/7,6/7]; B is the inscribed simplex, **exactly 1/6 of it**. This also CORRECTS my own THM-1150/1151 figure |B| = 0.003367 (a coarse-grid artifact, MISTAKE-176): the concentration at (1,2,3) is (2/21)/(1/343) = **98/3 = 32.67**, not 28.28. (III) PRIORITY, STATED PLAINLY: the single-run bound is **death-star-S58c's** (HYP-7737, pushed first). Their box-exit argument gives run ≤ 1/(7D) with D = max_i d_i global; the slack identity independently gives only run ≤ (1/7)/d_exit and, with mirror symmetry, run ≤ (1/7)/ρ where ρ = max(d_exit(r), d_exit(r*)) ≤ D. **Their bound is strictly stronger and their reduction is the one to use.** (IV) WHAT THIS ADDS TO THEIR RESIDUE: #runs is **always even** — 0 odd counts and 0 self-mirror runs over all 11,094 primitive directions with components ≤ 30 — so mirror symmetry gives #runs = 2P and their inequality **#runs ≤ 2D/3 is exactly D ≥ 3P**, P = number of mirror pairs. Verified with 0 violations; tight only at P = 1, with minimal witnesses (1,2,3), (1,6,7), (1,13,14), (1,16,17), (1,20,21), (1,27,28) — all of the form (1, D−1, D). My independent mirror-charged route lands on ρ ≥ 3P, which implies D ≥ 3P since ρ ≤ D; two independent reductions converging on the same constant 3
status: (I) PROVED — one line of algebra, 0 violations on an exhaustive lattice check. (II) PROVED — unimodular gap coordinates plus a corner-of-cube substitution, both boundary checks strict so there is no clipping or wraparound; Kronecker sampling converges to it (+0.33%, +0.05%, −0.07% at 2·10⁵, 10⁶, 4·10⁶ points). Also a CORRECTION to my own earlier measured |B|. (III) is an attribution note, not a claim: **death-star-S58c owns the single-run bound and the reduction.** (IV) the parity fact is VERIFIED (components ≤ 30, 0 exceptions) and makes the restatement D ≥ 3P rigorous; **D ≥ 3P itself is VERIFIED, NOT PROVED**, and so is ρ ≥ 3P. **Uniform r=5 remains OPEN**, and it is a rung below the r=6 sharp horn and two below the n=12 inverse theorem — the fleet's actual wall
source: kind-pasteur-2026-07-18-S128 (cont.83; owner: build on death-star's six-ray sojourn-max, think Kakeya needles)
depends_on:
  - THM-1154    # the six-ray / X-ray-transform reading
  - THM-1150    # whose "six boxes, |B| = 0.003367" is corrected here to six simplices, |B| = 1/343
related: [THM-1151, THM-1152, MISTAKE-173, MISTAKE-176]
script: 04-computation/simplex_slack_kakeya_kps_S128c83b.py, runs_parity_D3P_kps_S128c83c.py (+ .out)
credit: death-star-S58b/S58c (HYP-7735/7736/7737) — the bad set B, its mirror symmetry, the exact 2/21, the centre-hit run formula, and the box-exit single-run bound run ≤ 1/(7D) with its reduction to #runs ≤ 2D/3
---

# THM-1211 — the slack simplex

## Priority note, first

death-star-S58c (HYP-7737) proved the single-run bound **run ≤ 1/(7D)**, D = max_i d_i, by
a box-exit argument, and reduced the maximiser ceiling to the counting inequality
**#runs ≤ 2D/3**. That was pushed before this file and it is the stronger route: the
slack identity below independently yields only run ≤ (1/7)/d_exit, and with mirror
symmetry run ≤ (1/7)/ρ with ρ ≤ D. **Use their reduction.** What this file contributes is
the geometry underneath it, an exact |B|, a correction to my own earlier figure, and one
restatement of their residue.

## (I) The slack identity

In the ordering region g₍₁₎ ≤ g₍₂₎ ≤ g₍₃₎, B is cut out by g₍₁₎ ≤ 2/7, g₍₂₎−g₍₁₎ ≤ 2/7,
g₍₃₎−g₍₂₎ ≤ 2/7, g₍₃₎ ≥ 5/7. Put f₁ = 2/7 − g₍₁₎, f₂ = 2/7 − (g₍₂₎−g₍₁₎),
f₃ = 2/7 − (g₍₃₎−g₍₂₎), f₄ = g₍₃₎ − 5/7. Then

> **f₁ + f₂ + f₃ + f₄ = 6/7 − g₍₃₎ + g₍₃₎ − 5/7 = 1/7,  identically.**

Four nonnegative affine functions of constant sum cut out a **simplex**, never a box. Its
incentre has all four slacks equal to (1/7)/4 = **1/28**, and that point is exactly
g = (1/4, 1/2, 3/4) — **this is where death-star's 1/28 comes from**, and why the six
centres are the permutations of (1/4,1/2,3/4).

## (II) |B| = 1/343, exactly — and B is ⅙ of the box

Gap coordinates x = g₍₁₎, y = g₍₂₎−g₍₁₎, z = g₍₃₎−g₍₂₎ are unimodular, so
volume-preserving, and the constraints become x,y,z ≤ 2/7 with x+y+z ≥ 5/7. Substituting
u = 2/7−x, v = 2/7−y, w = 2/7−z ≥ 0 gives

> **{u, v, w ≥ 0,  u + v + w ≤ 1/7}**,

the corner of the cube [0,2/7]³ — a simplex of leg 1/7 and volume (1/7)³/6. The bound
u ≤ 2/7 is automatic since 1/7 < 2/7, and both boundary checks are strict
(x = 2/7−u ≥ 1/7 > 0 and g₍₃₎ = x+y+z ≤ 6/7 < 1), so nothing clips and nothing wraps.
Over six ordering regions:

> **|B| = 6 · (1/7)³/6 = (1/7)³ = 1/343 = 0.0029155.**

death-star's box-exit argument uses B's bounding box [1/7,2/7]×[3/7,4/7]×[5/7,6/7], of
side 1/7 — correct, and enough for the run bound. B is the **inscribed simplex, exactly
1/6 of that box**. Since run ≤ 1/(7D) is already tight at (1,2,3), the simplex cannot
improve the run bound; where it might bite is the residual counting inequality.

**Correction to my own canon (MISTAKE-176).** THM-1150/1151 report |B| = 0.003367. That
was a midpoint-grid artifact — grid error on a polytope is O(surface·h), tens of percent
even at 180³, and the estimates 0.00333 / 0.00283 / 0.00267 at 60³ / 120³ / 180³ do not
even converge monotonically. Kronecker sampling converges to 1/343. The sojourn values are
unaffected (2/21 is exact affine-cell arithmetic), but every "28.28× concentration" I have
published should read **98/3 = 32.67×**.

## (III) The slack route to a run bound, for the record

Within a cell f₄ is affine and decreasing at rate d_top; across a cell boundary the max
coordinate can only switch to a **slower** one (a switch at g_i = g_j with both ≥ 5/7
means g_j is larger just after, so d_j < d_i). So along a merged run f₄ is continuous,
decreasing, with slopes of decreasing magnitude, and every rate on the run is at least the
final one:

> d_exit · L ≤ Σ_cells d_c L_c = (total drop of f₄) ≤ 1/7,  so **L ≤ (1/7)/d_exit**,

and by mirror symmetry (the mirror run has the same length and its exit coordinate is this
run's bottom coordinate) **L ≤ (1/7)/ρ**, ρ = max(d_exit(r), d_exit(r*)). Verified: 0
violations over 11,557 primitive directions in [1,24]³. This is **weaker than
death-star's 1/(7D)** because ρ ≤ D. Recorded because the derivation is the slack identity
rather than the box, and because the ρ form is what produces (IV).

## (IV) What this adds to their residue: #runs ≤ 2D/3 **is** D ≥ 3P

Over all 11,094 primitive directions with components ≤ 30 that have a nonempty bad set:

> **0 directions have an odd number of runs, and there are 0 self-mirror runs.**

So mirror symmetry pairs the runs perfectly, #runs = 2P with P the number of mirror pairs,
and death-star's inequality is exactly

> **D ≥ 3P.**

Verified with 0 violations. Minimal witnesses:

| P | min D | witness | 3P |
|---|---|---|---|
| 1 | 3 | (1,2,3) | 3 |
| 2 | 7 | (1,6,7) | 6 |
| 3 | 14 | (1,13,14) | 9 |
| 4 | 17 | (1,16,17) | 12 |
| 5 | 21 | (1,20,21) | 15 |
| 6 | 28 | (1,27,28) | 18 |

Tight only at P = 1; the slack grows after that, and **every minimal witness has the form
(1, D−1, D)**. My mirror-charged route independently lands on **ρ ≥ 3P**, which implies
D ≥ 3P since ρ ≤ D — two reductions, arrived at separately, converging on the same
constant 3.

## Honest status

(I) and (II) are proved. (III) is a weaker rederivation of a bound death-star already owns.
(IV)'s parity fact is verified to components ≤ 30 and makes the restatement rigorous, but
**D ≥ 3P itself is verified, not proved**. **Uniform r=5 remains open**, and it sits a rung
below the r=6 sharp horn and two below the n=12 inverse theorem.

## Named next
- **Prove D ≥ 3P.** In the restated form this is: a rational needle that enters B in P
  mirror pairs must have max rate at least 3P. The (1, D−1, D) shape of every minimal
  witness is the obvious place to start — those are the directions with two nearly-equal
  fast coordinates and one slow one.
- The P = 1 base is death-star's c ≥ 3, which they get from a, b−a, c−b being positive
  integers summing to c — but only **at a centre hit**. Off a hit that argument is
  unavailable, and even P = 1 is open in general.
- The simplex is 1/6 of the box. The run bound cannot see the difference, but a counting
  bound should: entering a simplex is a strictly stronger event than entering its box.
