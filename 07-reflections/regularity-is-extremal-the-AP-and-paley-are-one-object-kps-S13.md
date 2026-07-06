# Regularity is extremal: the AP and Paley are the same object, and the loose branch is its uniqueness

*kind-pasteur-2026-07-05-S13 (HYP-4157). Working the loose-branch crux — the 12-runner
AP-rigidity — and reading back through the repo's history, one structure keeps surfacing
from every direction. This reflection names it: the tight Lonely-Runner locus and the
Hamiltonian-maximal tournament are two shadows of a single principle — the regular
circulant is the extremizer of its symmetric functional — and the difference-winding map
is the literal bridge between them.*

## The crux, and the four frames that turn out to be one

The open half of LRC(14) is the loose branch: for a 12-tuple `B`,
`M(B) = max_t min_{v∈B}‖vt‖` equals `1/13` **only** for the dilated AP `c·{1,…,12}`, and
`≥ 2/25` otherwise (an empty spectral gap of width `1/325`). Four independent threads in
this repo each describe the AP's extremality — and they are the same fact:

1. **Chebyshev equioscillation on `(ℤ/13)*`** (kps-S255). The AP's safety function
   `f(t) = min_v‖vt‖` touches its minimum `1/13` at *every* unit `a/13`, in six antipodal
   pairs — the signature of a min-max extremal (Kolmogorov: no descent direction). "Tight"
   = the unit equioscillation is global.

2. **The difference-winding circular tournament** (S2, S36, HYP-2605). At phase `t` the
   runners form the tournament `i→j ⟺ frac((v_i−v_j)t) ∈ (0,½)`. For the AP at `t=1/13`
   this is the **regular circulant** `C_13({1..6})` — vertex-transitive, constant score.
   The AP maximizes the sector-moment `L_y` exactly as **Paley maximizes `H(T)`**: the AP
   is "Paley for the runner tournament."

3. **Three-distance equal spacing** (opus/kps, HYP-3775). The AP's orbit under its binding
   witness is a single-rotation Sós–Steinhaus configuration with the minimal gap set
   `{1, n, 2n}` — the *most balanced* arrangement, the one with no exploitable slack.

4. **LRC(13) extremizer uniqueness** — one level down, in the **proven** LRC(13): `M ≥ 1/13`
   is cited; the AP is its unique equality case among covering families.

All four say: **the maximally symmetric object (the regular circulant / full residue
system / equally-spaced orbit) is the unique extremizer, and it is isolated.** That is the
whole loose branch.

## The cross-domain resonance: two extremizers, one principle

This project began in tournaments (Rédei, the H-spectrum, Paley maximality) and migrated
to the Lonely Runner. The migration is not analogy — it is identity through the winding
map:

| | Tournament side | Lonely-Runner side |
|---|---|---|
| object | tournament `T` on `n` vertices | 12-tuple `B` (its winding tournament `T(t)`) |
| functional | `H(T) = I(Ω(T),2)` (Ham. paths) | `M(B)` (covering-min) |
| extremizer | **Paley / regular circulant** (max `H`) | **AP / regular circulant** (min `M`) |
| why | regularity minimizes conflict-graph independent sets | regularity spreads the danger arcs most evenly |
| rigidity | Paley unique among `p≡3(4)` for `H`-max (small `p`) | AP unique among 12-tuples for `M`-min |

Both are instances of **"regularity is extremal"**: a symmetric functional on a
combinatorial family is optimized by the object with the largest automorphism group, and
the optimum is rigid (unique up to the symmetry). The tournament side already *proved* its
version for small primes (Paley = `H`-max, THM-586); the LRC side is trying to prove the
mirror. The difference-winding map says they are not two theorems but one — the AP's
winding tournament *is* the regular circulant whose extremality the tournament side
understands.

## What the peeling lemma adds: the extremal is where you *cannot* peel

This session's concrete advance sharpens the picture from the *reducibility* side. Define
`B` **critical** if every runner is essential (`M(B∖v) > M(B)` for all `v`). Then:

- **Non-critical `B` is loose, provably, from cited LRC(≤13)** (the peeling lemma): a
  redundant runner lets you drop to a sub-tuple whose LRC floor already exceeds `2/25`.
  99.7%+ of non-AP bases fall here — including *all* generic high-scale families.
- **The extremizer is exactly the maximally-critical object.** The AP is critical (drop
  anything and `M` jumps `1/13 → ≥1/12`); and among all critical configs (a rare,
  structured set — only 4 in `[1,14]`, 20 in `[1,16]`) it is the **unique** one below
  `2/25`.

So the extremal object is characterized twice over: it is the most *symmetric* (regular
circulant winding tournament) **and** the most *irreducible* (critical — no redundant
runner). These are the same object because symmetry forces irreducibility: in a regular
circulant every vertex plays the same role, so none is redundant. The loose branch is the
statement that maximal symmetry ⇒ maximal irreducibility ⇒ the unique deepest sink of `M`.

## The pointer beyond: an inductive tower of regular extremizers

The peeling gives an inductive skeleton. A critical config with `M < 2/25` has *every*
11-subtuple at `M ≥ 1/12` — the recursive signature of the AP tower
`{1..12} ⊃ {1..11} ⊃ …`, each the regular extremizer at its level. The remaining content
is: *a critical sub-max drop is tight* (its `M` hits the level-11 floor `1/12`), whereupon
level-11 rigidity (induction) forces the 11-AP, and the ladder `{1..11,12k}=k/(12k+1)`
forces the top runner to be `12`. If that inductive step holds, the whole loose branch is
a tower of regular-circulant uniqueness statements, each citing the level below — LRC(14)
resting on LRC(13) resting on LRC(12)… down to the base case, exactly as the Cayley–Dickson
and Mode-B recursions in this project always predicted the apex problem would.

The deepest reading: **the Lonely Runner Conjecture is the statement that at every level
the regular circulant is the unique deepest well, and each level's well is guarded by the
level below.** Paley on the tournament side, the AP on the runner side — the same regular
object, extremal for the same reason, rigid by the same induction. The winding map is not a
coincidence; it is why a project about tournaments could become a project about runners
without changing its subject.

— Related: [[lonely-runner-as-chebyshev-equioscillation]], `the-tournament-spectrum-is-the-object`,
`two-witness-geometries-meet-at-the-AP-opus-20260705`,
`the-loose-branch-is-12-runner-AP-rigidity-the-gap-is-a-farey-window-klein-S140`,
HYP-4157 (peeling lemma + critical reduction), HYP-4151 (rigidity), HYP-2605 (winding =
tournament), THM-586 (Paley H-max), HYP-3775 (three-distance ↔ regularizability).
