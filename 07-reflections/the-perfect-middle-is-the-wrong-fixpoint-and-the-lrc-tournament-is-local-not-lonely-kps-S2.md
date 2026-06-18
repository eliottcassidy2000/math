# The perfect middle is the wrong fixpoint, and the LRC tournament is local, not lonely

**Session:** kind-pasteur-2026-06-17-S2. **Results:** THM-525, HYP-2573..2576; OPEN-Q-108 sharpened
with GAP G2 and the plateau datum. Built by a workflow (structure investigators + perspective-diverse
reduction skeptics + an 8-theme tournament fan-out with adversarial refutation), every number
independently re-verified.

The prompt carried two instincts: *center the parked runner in the perfect middle of section 0*, and
*find a tournament that loneliness forbids an iso-class from reaching*. Following each to the end, the
instinct **inverted** — and the inversion is the content.

## Two wrong fixpoints

**The parked runner does not sit in the middle; it hugs the edge.** For the canonical hard family
`{1..11,13}∪{84m}` the optimum is `M = 7m/(84m+5)`, and at it the parked `w=84m` is *binding* — at
`frac(wτ*) ≈ 0.92`, pinned against core runner 5 at the very boundary of its danger band. The
optimizer extracts the gap by pressing `w` to the edge, not by centering it. The "perfect middle" is
real, but it is the wrong object: it is a *constructive certificate* (choose `τ∈G_C` that pushes
`‖wτ‖→1/2` while a core runner holds the `1/14` floor), giving a lower bound, never the true `M`. The
two `τ`'s — the optimizer's edge-binding `τ*` and the certificate's centering `τ` — must never be
conflated; their slacks differ by two orders of magnitude (`0.0072` vs `3/7`).

**Sharing sections does not make a configuration harder; it makes it easier.** The "each runner its
own section" picture — the perfect SDR `{1..13}` — looked like the easy, spread-out end. It is the
opposite: the SDR/AP is the **tight floor**, `M=1/14` exactly, and *every* collision (two runners
sharing a nonzero section) *raises* `M`. The only coordinate that tracks hardness is `z`, the number
of runners parked in section 0 (the multiples of 14) — and `z≥1` is exactly what forces the optimum
off the grid, where the section lens goes blind. Hardness is not a matching/Hall phenomenon (SDR-ability
is always the trivial "distinct & nonzero"); it is **band-avoidance**, a covering-design problem.

Both inversions point the same way: the easy structure the prompt wanted to lean on is genuinely there
(the easy 12-core's positive `meas(G_C)`), and it genuinely localizes all of LRC(14) onto a single
estimate — but the estimate it localizes onto *is* the open problem (OPEN-Q-108), not something weaker.
The three adversarial skeptics agreed: `Q⟹P` is honest and non-circular, the non-covering case is
closed unconditionally, ~105k covering sets show zero counterexamples — and "one estimate away"
undersells how large that one estimate is. The reduction's worth is not a proof; it is the **precise
naming of what remains**: GAP A (uniform `meas(G_C)≥c`, the coordinated `k≥3` regime) and a second,
distinct GAP G2 (the thin danger comb cannot *contain* the fat lonely set — transversality, nonempty
in every case, no general argument). Naming G2 separately from A is the session's cleanest structural
gain on the conjecture.

## The tournament is local, not lonely

The prompt's tournament instinct was the same shape as the dead overtaking map, and the eight-theme
hunt confirmed the deeper reason it dies. Across switch-parity, danger-interval, residue-orbit,
pinning, character-spectral, first-return, and pairwise-gap-load maps — dozens of non-obvious arc
rules built from crossings, parities, danger arrivals, section votes — the adversarial control
returned the same verdict every time: **forbidden-by-loneliness = 0.** The iso-class set realized at
the lonely-optimal `τ` equals the set realized at an arbitrary `τ`. The choice of lonely time never
carves out a class. Loneliness is a *global-minimum* condition on one scalar (`min_v ‖vτ‖`); the
pairwise switch combinatorics these maps read are blind to it.

But one map's failure was structured enough to be a theorem. The **difference-winding** map
`i→j ⟺ frac((v_i−v_j)τ) ∈ (0,1/2)` *is* the circular tournament on the phase points `a_v=frac(vτ)`
on `ℝ/ℤ`. Tie-free, it is always a **local (round) tournament** — and so the unique maximal-`H`
*non-local* class on five vertices, score `(1,2,2,2,3)`, `c3=4`, `H=15`, is unreachable. The signed
danger-arrival map and the section-rotational map forbid the *same* class. The forbiddance is real and
survives every refutation — but its cause is **circle geometry**, not loneliness. The LRC lives on a
circle; tournaments read off a circle are local; non-local tournaments never appear. This is the same
note the project keeps sounding from the other side — circulant/Paley tournaments, the `C₁₄`
rotational structure, the `(ℤ/14)*` action that reverses sections (`a=−1` = complement = `T↦T^op`).
The LRC-derived tournament is forced into the rotational/local world by the geometry of the loop, and
that geometry — not the runners' loneliness — is what forbids the classes off it.

So the honest answer to "find the tournament loneliness forbids a class from" is: there isn't one in
the τ-selection sense, and the search for it is now a documented dead direction. The answer to "find a
tournament that *can't* reach some class an arbitrary tournament could" is: yes — the circular/danger
maps, forced local by the circle, miss the non-local classes. The difficulty migrated, again, from the
object the prompt named (loneliness selecting tournaments) to the geometry underneath it (the loop
making them local) — a genuine condensation, and not the obstruction the hunt set out to find, because
the obstruction is the circle itself.

## The shape of it

Twice the prompt's fixpoint was the boundary of the wrong basin — the optimizer pins the parked runner
to the edge, sharing relaxes toward the floor — and twice the inversion was the real structure. The
easy cases do dominate the hard ones, but only by handing the whole problem, cleanly, to the one
estimate that is the problem. And the tournament the loop generates was never going to feel loneliness;
it only ever felt the circle. The next move is the one the geometry suggests: bound the speeds (Tao /
MSS), turn the 105k-set scan into a finite certificate, and let the circle's own rotational structure —
the thing that makes the tournament local — be the thing that pins `G_C` away from the parked runner's
comb. [[lrc14-thread]] · [[triangle_foundation]]
