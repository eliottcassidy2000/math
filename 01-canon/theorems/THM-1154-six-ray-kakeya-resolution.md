---
id: THM-1154
title: THE SIX-RAY SOJOURN-MAX CONFIRMED AND EXPLAINED — (1,2,3) is maximal because it is the unique INCREASING representative of a six-fold symmetric orbit, not because it is uniquely resonant. (I) THE KAKEYA READING: the sojourn of a geodesic of direction d in the bad set B is exactly the **X-ray transform** of 1_B in direction d, so "which d maximises sojourn" is a Kakeya-type question — how much of a set a line can see. B is six boxes centred at the six permutations of (1/4,1/2,3/4), and a box centred at σ(1,2,3)/4 is threaded **lengthwise** by the direction σ(1,2,3). Hence exactly **six maximising rays**, which is death-star-S58b's six-ray sojourn-max. (II) ALL SIX TIE EXACTLY: sojourn = **0.095230 = 2/21** for every one of (1,2,3), (1,3,2), (2,1,3), (2,3,1), (3,1,2), (3,2,1) — 28.28×|B| each, identical to six digits. (III) EACH RAY HITS ITS OWN CENTRE, exactly and at the same parameter: direction σ(1,2,3) passes through σ(1,2,3)/4 at **u = 3/4** for all six σ. That is the head-on entry that makes the chord maximal. (IV) THE ORDERING CONSTRAINT SELECTS ONE: killers satisfy k₂ < k₃ < k₄, which forces d₂ < d₃ < d₄, and **exactly one of the six permutations is increasing** — (1,2,3). So the maximiser is picked out by a symmetry constraint on the problem, not by any special arithmetic of the vector. (V) SCALINGS PRESERVE IT: (2,4,6), (4,2,6), (6,4,2) all give 28.29×, so the whole orbit ℝ₊·{σ(1,2,3)} ties and the increasing-representative argument applies at every scale. This closes, structurally, the "why (1,2,3)" question that my beat, tube, standoff and centre-hitting framings all failed to answer
status: (I) the X-ray identification is exact. (II),(V) MEASURED at grid resolution — 0.095230 against 2/21 = 0.0952381. (III) PROVED — the centre-hitting at u = 3/4 is an exact rational identity for all six permutations. (IV) is the resolution, and it is structural rather than computational. **What remains unproved is the Kakeya sup bound itself**: that no direction outside the six-ray orbit exceeds 2/21. THM-1152 already showed 96–117 other directions hit a centre, reaching only 4–10×, so the evidence is strong, but a sup bound over all d is not established. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.83; owner: build on death-star's six-ray sojourn-max, think Kakeya needles)
depends_on:
  - THM-1152    # which enumerated the compatible relations and left the invariant unidentified
  - THM-1150    # the six-box geometry
related: [MISTAKE-173, THM-1148, THM-1211]
script: 04-computation/six_ray_kakeya_kps_S128c83.py (+ .out)
---

> **CORRECTION (same session, THM-1211).** The "×|B|" column below uses the measured
> |B| = 0.003367 I carried from THM-1150/1151. That figure is a coarse-grid artifact.
> B is six **3-simplices** (not boxes) and **|B| = (1/7)³ = 1/343 = 0.0029155 exactly**,
> so every "28.28×" here should read **98/3 = 32.67×**. The sojourn values themselves
> (0.095230 = 2/21) are unaffected — only the ratio to |B| was wrong.
>
> **CREDIT/PRIORITY.** death-star-S58b proved the exact 2/21 tie on the (1,2,3)-orbit by
> exact affine-cell arithmetic, and proved the ceiling on the 1-mirror-pair class, in
> parallel with this. Sections (II) and (V) below are independent confirmation, not
> priority. What is net-new here is (IV), the ordering-selects-the-representative reading.

# THM-1154 — six rays, and why (1,2,3) wins

## (I) The Kakeya reading

The sojourn of a geodesic of direction d in B is the **X-ray transform** of 1_B in direction
d. So "which direction maximises sojourn" is a Kakeya-type question: how much of a set can a
line see? codex's "four-residue multiplier Kakeya cone" (via THM-1148) is the same lens.

B is six boxes centred at the six permutations of (1/4, 1/2, 3/4). A box centred at
σ(1,2,3)/4 is threaded **lengthwise** by the direction σ(1,2,3) — so there are exactly **six
maximising rays**, which is precisely death-star-S58b's six-ray sojourn-max.

## (II) All six tie, exactly

| direction | sojourn | ×\|B\| | increasing? |
|---|---|---|---|
| **(1,2,3)** | 0.095230 | 28.28 | **YES** |
| (1,3,2) | 0.095230 | 28.28 | no |
| (2,1,3) | 0.095230 | 28.28 | no |
| (2,3,1) | 0.095230 | 28.28 | no |
| (3,1,2) | 0.095230 | 28.28 | no |
| (3,2,1) | 0.095230 | 28.28 | no |

Identical to six digits — 2/21 for every one.

## (III) Each ray hits its own centre, at the same parameter

| direction | its centre | hit at |
|---|---|---|
| (1,2,3) | (1/4, 1/2, 3/4) | **u = 3/4** |
| (1,3,2) | (1/4, 3/4, 1/2) | **u = 3/4** |
| (2,1,3) | (1/2, 1/4, 3/4) | **u = 3/4** |
| (2,3,1) | (1/2, 3/4, 1/4) | **u = 3/4** |
| (3,1,2) | (3/4, 1/4, 1/2) | **u = 3/4** |
| (3,2,1) | (3/4, 1/2, 1/4) | **u = 3/4** |

Exact rational identities, all at the same u. That head-on entry is what makes the chord —
and hence the sojourn — maximal.

## (IV) The ordering constraint selects one

The killers satisfy k₂ < k₃ < k₄, so d₂ < d₃ < d₄. **Exactly one of the six permutations is
increasing: (1,2,3).**

So (1,2,3) is not arithmetically privileged at all. It is the unique admissible
representative of a six-fold symmetric orbit, and the symmetry is broken by the problem's own
ordering convention rather than by number theory. That is why my beat, tube, standoff and
centre-hitting framings all failed: each looked for something special *about the vector*,
and there is nothing special about it.

## (V) Scale-invariance of the whole orbit

(2,4,6), (4,2,6), (6,4,2) all give 28.29×, so the full ray orbit ℝ₊·{σ(1,2,3)} ties at every
scale and the increasing-representative argument applies throughout — matching the
scale-invariance measured in THM-1151.

## Honest status

This resolves *why* (1,2,3), structurally, which four earlier framings of mine did not. What
it does **not** do is prove the Kakeya sup bound — that no direction outside the six-ray
orbit exceeds 2/21. THM-1152 found 96 (death-star: 117, counting permutations) other
directions that hit a centre, and they reach only 4–10×, so the evidence is strong; but a
supremum over all d is still not established. **Uniform r=5 remains open.**

## Credit

The six-ray structure is death-star-S58b's (MISTAKE-173 salvage). This confirms it, supplies
the exact tie and the common hit parameter u = 3/4, and adds the observation that the
ordering constraint is what selects (1,2,3) from the orbit.

## Named next — SUPERSEDED IN-SESSION BY THM-1211

The piece named here was the Kakeya sup bound: **sup over directions of the X-ray transform
of 1_B is attained on the six rays.** THM-1211 takes that up and reduces it, via the slack
identity f₁+f₂+f₃+f₄ ≡ 1/7, to the single combinatorial inequality **W'' ≤ 2/3** with
W'' = Σ_runs 1/max(d_exit(r), d_exit(r*)) — verified for all primitive d ∈ [1,24]³ and
tight at (1,2,3). That inequality does not mention mirror pairs, so it subsumes both
death-star's proved 1-pair class and their enumerated ≥2-pair residue. It is still
**verified, not proved**; see THM-1211's named next.

The chord-geometry guess recorded above is also worth retiring: B is six simplices, so
"longest chord of a box along its diagonal" is not the right picture.
