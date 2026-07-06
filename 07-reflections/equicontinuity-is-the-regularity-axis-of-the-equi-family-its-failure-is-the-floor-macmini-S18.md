# Equicontinuity is the regularity axis of the equi-family — its failure IS the density floor, and it is dual to equidistribution

*mac-mini-2026-07-06-S18 (HYP-4462). Owner: see how equicontinuity relates to the
prior equidecomposability / equinumerosity / equioscillation / equidistribution
work; keep integrating and reframing creatively. This note places equicontinuity —
the least-explored of the "equi-" family (4 repo mentions vs 374 for
equidistribution) — as the *regularity axis*: the meta-property whose failure is
exactly the density floor, and whose failure is the dual of equidistribution.
Verified: the equioscillation count (`…S18`); builds on kps-S26 (compactness
bypass), kps-S255 (Chebyshev equioscillation), codex-S257 (forgetting triad).*

## The AP is *equi-everything* — one point, five symmetries

Every "equi-" notion the project developed is a face of the AP's maximal symmetry:

| notion | the AP's face | verified |
|---|---|---|
| **equioscillation** | `f_AP` touches its max `M=1/n` at **exactly φ(n) points = the units of ℤ/n** (12=φ(13); 6=φ(7)) — the unique Chebyshev minimax extremal | this note |
| **equidistribution** | at `t=1/n` the runners are the roots of unity, perfectly uniform (min star-discrepancy) | opus HYP-4074 |
| **equidecomposability** | the danger arcs *tile* the circle at the covering threshold (scissors class) | S599 |
| **equinumerosity** | maximal relation lattice / additive energy (the theta-sum) | opus HYP-4396 |
| **complement symmetry** | mirror `v ↔ n−v` about `t=1/2` (the `−1` of `(ℤ/n)*`) | kps-S255 |

Measured (equioscillation count `= #{t : f_S(t)=M(S)}`): **AP = φ(n) (maximal);
every other config = 2–4.** The doubled-apex and the n=7 gap member touch at only
2 points. So the equi-family is the structural fingerprint of the tight locus — the
AP is where all five coincide, and they *drop sharply* the moment you leave it.

**These are all structural and all universal** — they mark the AP special at every
`n`. By the codex-S257 "forgetting triad," they form a hierarchy of what a quotient
remembers: `equidistribution ⊂ equinumerosity ⊂ equidecomposability` (same limiting
measure ⊂ same count ⊂ same scissors class; each gap a resonance / Dehn invariant).

## Where equicontinuity sits: the regularity axis, and its failure is the floor

Equicontinuity is **not** another equivalence relation — it is the *regularity*
condition that decides whether the equivalences are **uniform**. And here is the
whole point (kps-S26): **`M` is NOT equicontinuous.** On the compact projective
space of directions, `M(v)=M(cv)` is continuous, but its modulus is
`L(d) ~ height(d)/13` — it *grows with height*, oscillating at frequency `~height`
near the tight locus. That is why the compactness/Dirichlet bypass fails, and why
the "height threshold" is non-uniform.

This is the analytic form of my S17 finding that **the floor is quantitative, not
structural**: the equi-symmetries (equioscillation, equidistribution, …) are
*equicontinuous invariants* — they see only the coarse structure — so they are
**necessary but not sufficient**. The density floor lives exactly in the part `M`'s
non-equicontinuity resolves and the equi-invariants forget: the fine
height-oscillation. The second gap being **n-specific** (nonempty n=7,8, empty
n=13) is the signature that no equicontinuous (structural) invariant can close it —
because such invariants are `n`-uniform, but the floor is not.

## The creative reframe: non-equicontinuity of `M` IS equidistribution of the far runner

The deepest link is a **duality**. `M` is non-equicontinuous *because* a far
(large) runner `v_far` oscillates at frequency `~v_far` — pointwise, `‖v_far t‖`
swings rapidly, so `M` is height-sensitive. But that *same* fast oscillation
**equidistributes**: `∫(1−g_β(v_far t))dt → 1−2β`. So

> **non-equicontinuity of `M` (pointwise, jagged) ⇔ equidistribution of the far
> runner (averaged, smooth)** — two faces of one fast oscillation.

The average is the good object: `safe(S,β)=∫∏(1−g_β(v_i t))dt` is *smooth*, and my
two-scale decorrelation (HYP-4402, `safe(A∪NB) → safe(A)safe(B)`) is precisely the
equidistribution/averaging side — it is why the multi-scale case *closed*
rigorously. The **max** `M` is the bad object: jagged, non-equicontinuous, and a
gap member would be a *narrow jagged spike* of `M` into `(1/13,2/25)` that the
smooth average cannot see. The floor is the statement that **no such spike survives
the averaging except as the AP's exact tiling** — the tension between `M`'s
non-equicontinuity (which permits spikes) and `safe`'s equicontinuity (which is
AP-rigid at 0).

This reframes the whole programme:
- **Structural / equi-invariant side** (equioscillation, equidistribution,
  equidecomposition, roots of unity): characterizes the AP. *Done, necessary.*
- **Non-equicontinuous / averaging side** (the far-runner oscillation, the
  height-modulus, the fine arithmetic): *is* the floor. The decorrelation lemma
  won the part of it that averages cleanly (multi-scale); the residue is the part
  where the spike and the average meet (single-scale, n-specific).

## What it buys — a sharper target

The floor is the **non-equicontinuity budget**: the maximum `M`-spike a covering
`12`-family can sustain into the gap before the smooth `safe` measure (the
equidistributed average) forces it back to `≥ 2/25` — and at `n=13` that budget is
`0` (only the AP's exact tiling reaches it), while at `n=7` it is positive (the
`5/33` spike survives). Concretely, the leave-one-out alignment (S17) is the
*equidecomposition* face — each 11-subfamily hole must tile into the dropped
runner's arcs — and its `1/(n(2n-1))` hole-width is exactly the non-equicontinuity
scale `1/height` at which `M` resolves. **Equidecomposition (tiling) and
equicontinuity (resolution) meet at the hole-width: the floor is that the tiling
tolerance `1/325` at n=13 is finer than any non-AP arc-lattice can align.**

## Net

- The five equi-notions are five faces of the AP's maximal symmetry
  (equioscillation at φ(n) units is the new, cleanest one). They characterize the
  tight locus — structural, universal, necessary.
- **Equicontinuity is the regularity axis, not an equivalence**: `M`'s
  non-equicontinuity (modulus `~height`, kps-S26) is *the density floor*, and is
  why the structural equi-invariants are necessary-not-sufficient (my S17) and the
  gap is n-specific.
- **Duality:** non-equicontinuity of `M` ⇔ equidistribution of the far runner. The
  averaged object `safe` is equicontinuous and AP-rigid (my decorrelation closed
  the multi-scale case); the max `M` is jagged and permits the gap spike. The floor
  is the meeting of the two, at the hole-width `1/(n(2n-1)) = 1/height` scale.

## Pointers

- equioscillation count script (`…S18`); kps-S26 (compactness bypass /
  non-equicontinuity), kps-S255 (Chebyshev equioscillation at units), codex-S257
  (forgetting triad), S599/S617 (equidecomposability vs equinumerosity), opus
  HYP-4074 (discrepancy inversion); mac-mini HYP-4402 (decorrelation), HYP-4452
  (leave-one-out alignment / quantitative floor), HYP-4412 (three-gap, n-corrected).
