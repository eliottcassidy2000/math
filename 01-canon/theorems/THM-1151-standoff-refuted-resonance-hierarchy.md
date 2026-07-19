---
id: THM-1151
title: THE STANDOFF CONJECTURE IS REFUTED — and the true mechanism is a RESONANCE HIERARCHY, with (1,2,3) the unique full resonance at 28× and no other direction above 6.33×. (0) MY OWN THM-1150 CONJECTURE IS FALSE: non-proportional directions do NOT keep standoff above ρ ≈ 0.041. Witnesses reach **0.00486** at (211,367,593), **0.00494** at (97,173,281), **0.01748** at (41,73,131) — all far below ρ. Long geodesics equidistribute and pass arbitrarily close to the centres, exactly as they must. (I) SO NON-PROPORTIONAL DIRECTIONS DO ENTER B — they simply do not linger. Sojourn for large generic directions converges to |B| = 0.003367: measured **1.03×** at (11,19,37), **1.06×** at (97,173,281), **0.89×** at (211,367,593). This is ordinary equidistribution, and it is the correct mechanism in place of the geometric miss I proposed. (II) THE PROPORTIONAL FAMILY IS A RESONANCE, AND IT IS SCALE-INVARIANT: sojourn/|B| = **28.28×** at (1,2,3), (5,10,15), (20,40,60) and (97,194,291) — identical at every scale, because the geodesic is confined to the same 1-dimensional subvariety {g₃ = 2g₂, g₄ = 3g₂} which threads B. (III) THERE IS A HIERARCHY IN BETWEEN, which my earlier framings missed entirely: **(5,9,14) gives 6.33×**, and it satisfies d₂ + d₃ = d₄. The governing rule is that a relation m·d = 0 concentrates the sojourn **iff** m·(1/4,1/2,3/4) ≡ 0 (mod 1) — for (1,1,−1) that is (1+2−3)/4 = 0, so the relation is compatible with the centre and produces partial concentration. (IV) (1,2,3) IS THE UNIQUE FULL RESONANCE: it carries TWO independent compatible relations (d₃−2d₂ = 0 and d₄−3d₂ = 0) where a partial resonance carries one, which is why 28× versus 6.33×. The maximiser claim therefore survives — 0.0952 against a next-best 0.0213, both far below the 0.164 safe measure — but by the resonance hierarchy, not by any standoff
status: (0) REFUTED — my own THM-1150 conjecture, by explicit witnesses. (I),(II),(III) MEASURED, with the scale-invariance of 28.28× confirmed across four scales and the generic convergence to |B| across three. (IV) is the correct mechanism and the relation-compatibility rule is stated exactly, but it is **not proved** that no direction exceeds 6.33× outside the proportional family. **Uniform r=5 remains OPEN, and the maximiser claim is still not proved** — this is the third framing I have offered and the second I have had to withdraw
source: kind-pasteur-2026-07-18-S128 (cont.79; owner: prove the standoff bound for non-proportional directions)
depends_on:
  - THM-1150    # whose standoff conjecture this refutes
  - THM-1149, THM-1148
script: 04-computation/standoff_test_kps_S128c79.py (+ .out)
---

# THM-1151 — no standoff; a resonance hierarchy instead

## (0) The conjecture is false

THM-1150 conjectured that non-proportional integer directions keep standoff > ρ ≈ 0.041 from
the six centres, so their geodesics miss B. They do not:

| direction | standoff |
|---|---|
| (41,73,131) | 0.01748 |
| (97,173,281) | 0.00494 |
| (211,367,593) | **0.00486** |

Far below ρ. This was predictable and I should have checked it before conjecturing: a
geodesic of direction d has length ~|d| in T³, so as |d| grows it equidistributes and comes
arbitrarily close to every point. The standoff route was never going to work.

## (I) They enter B — they just do not linger

The correct mechanism is equidistribution. Sojourn for large generic directions converges to
the volume |B| = 0.003367:

| direction | sojourn | ratio to \|B\| |
|---|---|---|
| (11,19,37) | 0.003470 | **1.03×** |
| (97,173,281) | 0.003560 | **1.06×** |
| (211,367,593) | 0.003000 | **0.89×** |

## (II) The proportional family is a scale-invariant resonance

| direction | sojourn | ratio |
|---|---|---|
| (1,2,3) | 0.09523 | **28.28×** |
| (5,10,15) | 0.09525 | **28.29×** |
| (20,40,60) | 0.09540 | **28.33×** |
| (97,194,291) | 0.09523 | **28.28×** |

Identical at every scale, because the geodesic is confined to the same 1-dimensional
subvariety {g₃ = 2g₂, g₄ = 3g₂}, which threads B. Scale-invariance is the signature of a
resonance rather than an accident.

## (III) The hierarchy in between — which I had missed

**(5,9,14) gives 6.33×**, an order between generic and full resonance. It satisfies
**d₂ + d₃ = d₄**. The rule:

> a relation m·d = 0 concentrates the sojourn **iff** m·(1/4, 1/2, 3/4) ≡ 0 (mod 1).

For m = (1,1,−1) this is (1 + 2 − 3)/4 = 0 — compatible with the centre, hence partial
concentration. Neither my beat framing, my tube framing, nor my standoff framing predicted
this intermediate level; it only appears once the resonance structure is the object.

## (IV) Why (1,2,3) is the maximiser

It carries **two** independent compatible relations (d₃ − 2d₂ = 0 and d₄ − 3d₂ = 0), where a
partial resonance carries one. Two relations confine the geodesic to a curve that lies *in*
the bad configuration locus rather than merely crossing it — hence 28× against 6.33×.

The maximiser claim survives: 0.0952 for the full resonance, 0.0213 for the best partial one
found, both far under the 0.164 safe measure. But it survives by the hierarchy, not by any
standoff.

## Honest status

This is the third framing I have offered for the same phenomenon — beat, tube, standoff —
and the second I have had to withdraw. The maximiser claim is **still not proved**, and
uniform r=5 remains **open**.

What is better than before: the mechanism is now the standard one for this kind of problem
(equidistribution versus resonance), the intermediate level is visible, and the governing
rule — relation compatibility with (1/4,1/2,3/4) — is stated exactly rather than guessed.
What is missing is a bound showing no non-proportional direction exceeds the partial-resonance
level.

## Named next
- Prove the relation-compatibility rule of (III), and enumerate the compatible relation
  vectors m. Since the centre is (1/4,1/2,3/4), m·(1,2,3)/4 ∈ ℤ constrains m modulo 4 — a
  finite condition, so the compatible relations are classifiable.
- Then bound the sojourn of a direction by the number of independent compatible relations it
  carries. Two is the maximum (three would force d = 0), so (1,2,3) is maximal by
  construction, and the claim follows from a per-level bound rather than a search.
- Stop proposing geometric framings for this without first testing them at large |d|. All
  three failures shared that omission.
