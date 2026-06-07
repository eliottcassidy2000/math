---
id: THM-435
name: the-Hamming-cube-K3box3-is-angle-rigid-at-81-unit-distances
status: PROVED (elementary trig + collision lemma; verified numerically LEM-A/B/C)
date: 2026-06-07
session: monad-explorer-2026-06-07-S6
depends_on:
  - THM-432   # K3^[]3 = H(3,3) is the n=27 product tie (81 = 3*27), 6-regular
  - THM-433   # avgdeg additive under []; products beat 3N only at N=32
  - THM-431   # N* in [25,28]; OPEN-Q-057 frontier
relates_to:
  - OPEN-Q-057 # is u(27)=81 (=> N*=28) or u(27)>81 (=> N*<=27)?
---

# THM-435: the Hamming cube K₃^□3 = H(3,3) is ANGLE-RIGID at 81 unit distances

## Context

OPEN-Q-057: `N* =` smallest `N` with `u(N) > 3N`. Proven `N* ∈ [25,28]`. The single
construction that *ties* the threshold at `n = 27` is the Hamming cube
`H(3,3) = K₃^□3` — the Minkowski sum of three unit equilateral triangles at generic
angles: **27 points, 6-regular, exactly `81 = 3·27` unit distances** (THM-432). The
obvious next move toward `u(27) > 81` (which would give `N* ≤ 27`) is to **tune the
three rotation angles** so the cube acquires *accidental* (non-product) unit
distances. THM-435 proves this is **impossible**: any angle choice that would create
an extra unit distance simultaneously collapses two of the 27 points.

## Statement

Let `C(t₂,t₃)` be the realization of `K₃^□3` as
`{ v₁ + R(t₂)v₂ + R(t₃)v₃ : v_i ∈ T }`, where `T = {0, e^{i·0}, e^{i·60°}}` is the
unit triangle and `R(t)` is rotation by `t` (factor 1 unrotated, WLOG). Then

```
   for every (t₂,t₃):   (#distinct points = 27)  ⟹  (#unit distances = 81).
```

Equivalently: **no choice of angles yields 27 distinct points with more than 81 unit
distances.** The 3N tie at `n = 27` is *angle-rigid* within the cube family, so this
family cannot witness `u(27) > 81`.

## Proof

The `81` **product edges** — pairs of points differing in exactly one triangle
coordinate by a triangle-edge vector — have length `1` for *all* angles, since a
triangle edge is a unit vector and rotation preserves length. (Each point has one
neighbour for each of `3` coordinates × `2` other vertices = degree `6`;
`27·6/2 = 81`.) So `U ≥ 81` always, and any *extra* unit distance is a pair differing
in **≥ 2** coordinates whose displacement — a sum of one unit vector per differing
factor — has length exactly `1`.

**(A) Two differing factors `i < j`.** Displacement `= e^{iα} + e^{iβ}` with
`α ∈ {t_i + 60°k}`, `β ∈ {t_j + 60°k}` (the factors' hex-direction sets;
`t₁ := 0`). `|e^{iα}+e^{iβ}| = 1 ⟺ ∠(α,β) = 120° ⟺ t_i − t_j ≡ 0 (mod 60°)`.

**(B) Three differing factors.** Dividing by factor 1's unit vector, the length-1
condition is
```
        cos u + cos w + cos(u − w) = −1,
        u ≡ t₂,  w ≡ t₃,  u − w ≡ t₂ − t₃   (mod 60°).
```
Writing it as `A cos w + B sin w = C` with `A = 1+cos u`, `B = sin u`,
`C = −1−cos u`: if `u ≢ 180°` then `A,B` are not both zero, `√(A²+B²) = √(2(1+cos u))`,
and `C/√(A²+B²) = −cos(u/2)`, with phase `atan2(B,A) = u/2`. Hence
`cos(w − u/2) = −cos(u/2) = cos(180° − u/2)`, giving **`w = 180°` or `w = u − 180°`**.
If `u ≡ 180°` then `A = B = C = 0` and the equation holds for *all* `w` (degenerate
line). So the **complete solution set** is
```
        { u ≡ 180° } ∪ { w ≡ 180° } ∪ { w ≡ u − 180° }   (mod 360°),
   i.e.  t₂ ≡ 0   ∨   t₃ ≡ 0   ∨   t₂ − t₃ ≡ 0           (mod 60°).
```

**Collision lemma.** Each of `t_i − t_j ≡ 0 (mod 60°)` forces `< 27` distinct points:
both `T` and any `60°`-rotation of `T` lie in the same Eisenstein (triangular)
lattice, so two such factors have a Minkowski sumset of size `< 9`, and holding the
third coordinate fixed yields two coincident cube points. (Verified: the six
collision loci give `18` or `21` distinct points, never `27`.)

Cases (A) and (B) show every extra unit distance lives on a locus where some
`t_i − t_j ≡ 0 (mod 60°)`, which collapses the point set. Therefore, whenever the
`27` points are distinct, no extra edge exists and `U = 81`. ∎

## Verification

`04-computation/unit_distance_cube_angle_rigidity_monad_s6.py`
(`05-knowledge/results/…_s6.out`):
- **LEM-B** — dense `(u,w)` scan (4200 root crossings): every solution of
  `cos u+cos w+cos(u−w)=−1` lies on `{u≡180}∪{w≡180}∪{w≡u−180}`; 0 rogues. Plus
  100 000-trial algebraic check that `w=180°` and `w=u−180°` solve it identically.
- **LEM-A** — all six collision loci give 18 or 21 distinct points (< 27).
- **LEM-C** — global `(t₂,t₃)` grid (~90 000 cubes): max `U` over 27-distinct cubes
  `= 81`; **0** cubes with 27 distinct points and `U > 81`. Independent of the
  earlier 56 882-sample scan (`…_s6c.py`) which also peaked at 81.

## Scope and significance

- **Negative result, sharp.** THM-435 rules out the *most symmetric / most obvious*
  route to `u(27) > 81`. It does **not** prove `u(27) = 81` — a general (non-cube)
  Moser-lattice blob could still beat 81 (AMP's upper bound is only `u(27) ≤ 90`).
- **Complements THM-432/433.** Those showed *products* tie at 81 and first beat `3N`
  at `N = 32` (avgdeg additivity). THM-435 shows even *non-product* (angle-tuned, so
  irreducible) perturbations of the cube are stuck at 81: the cube's tie is not a
  generic-angle artifact but a genuine rigidity. This closes the "just tune the
  cube" idea and reinforces the structural story that `N*`'s extremal graph, if
  `≤ 27`, must be a genuinely irregular blob — *neither* product *nor* a tuned cube.
- **Why no accidental edge survives.** A 3-factor accidental edge needs three
  triangle-edge unit vectors to sum to a unit vector; the only way (in the
  hex-direction sets) is for two of the three factors to align mod 60°, which is
  exactly the collision condition. Accidental density and point-distinctness are in
  direct conflict for the equilateral cube.

## Files
- `04-computation/unit_distance_cube_angle_rigidity_monad_s6.py` (+ `.out`)
- `04-computation/unit_distance_augment_cube_monad_s6.py`,
  `04-computation/unit_distance_augment_cube_monad_s6c.py` (landscape: augmentation
  `H(3,3)+1` adds only `+2` ⟹ `u(28) ≥ 83` via this realization; greedy growth
  `27→30` gives `81,83,86,90`).
