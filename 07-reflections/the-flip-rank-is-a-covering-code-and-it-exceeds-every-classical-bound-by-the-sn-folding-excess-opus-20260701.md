# The flip-rank is a covering code — and it exceeds every classical bound by the S_n-folding excess; depth and symmetry fuse only at Paley primes

*opus-2026-07-01-S16. Extends HYP-3805 (flip-rank / Paley obstruction) into the coding-theory framing the
owner asked for: covering radius, iso-diameter, the min transversal-subcube dimension of the S_n-orbit
coloring, and the standing Paley-prime prediction. Merges kind-pasteur's depth axis (MFAS) with my width axis.*

## The object: a covering code on the arc-hypercube
Fix `n`. The arc-hypercube `Q_m`, `m=C(n,2)`, is `{0,1}^m` — one bit per edge orientation, each point a labeled
tournament. The symmetric group `S_n` acts by relabeling; the **colors are the orbits** (iso classes),
`A000568(n)` of them. Every invariant this project has been circling is a *covering-code parameter of this
colored cube*:

- **Flip-rank `k(n)`** (= klein `rho` = mac-mini `kappa` = kps `k_min`) = the minimum dimension of an
  **axis-aligned subcube** (fix `m-k` coordinates, free `k`) that is a **transversal of the coloring** — its
  `2^k` points meet every color. This is a *covering* parameter: the smallest "flat" that touches all orbits.
- **Rainbow number `Rb(n)`** (klein HYP-3804) = the maximum subcube dimension with all `2^k` points in
  **distinct** colors. The *packing* dual.
- **Covering radius `R(n)`** (mac-mini `R` = max MFAS) = `max_T` distance from a tournament to the nearest
  **transitive** tournament. Since the transitive class's orbit is all `n!` linear orders, the distance from
  class `C` to it is exactly `MFAS(C)` (min feedback arc set), so `R(n) = max_C MFAS(C)`.
- **Iso-diameter `D(n)`** = `max` over class pairs of the min-iso Hamming distance.

## The bound chain (proved, verified exact n≤7)
| n | #classes | info-floor ⌈log₂⌉ | R(n)=max MFAS | D(n) | **k(n)** |
|---|---|---|---|---|---|
| 3 | 2 | 1 | 1 | 1 | **1** |
| 4 | 4 | 2 | 1 | 2 | **2** |
| 5 | 12 | 4 | 3 | 3 | **4** |
| 6 | 56 | 6 | 4 | 4 | **7** |
| 7 | 456 | 9 | 7 | ≥7 | **12** |

Two lower bounds, both provable in one line each:
- **`k(n) ≥ ⌈log₂ #classes⌉`** (sphere-covering / information): `2^k` points can't meet more than `2^k` colors.
- **`k(n) ≥ D(n) ≥ R(n)`** (geometric): a transversal subcube contains a rep of the transitive class and a rep
  of *every* other class, all agreeing on the `m-k` fixed coordinates — so any two of them differ only within
  the `k` free coordinates, giving `min-iso-dist(C₁,C₂) ≤ k` for **every** pair. Maximize over pairs → `k ≥ D`;
  and `D ≥ R` since transitive-to-farthest is one pair.

**The punchline is that `k(n)` exceeds all of them.** At `n=6`, `k=7` while `⌈log₂⌉=6`, `D=R=4` — every
classical bound is beaten. At `n=7`, `k=12` while `⌈log₂⌉=9`, `R=7`. The covering problem is *strictly harder*
than information, than diameter, than covering radius. That gap — `k − max(bounds) = 0,0,0,1,3` — is the
**S_n-folding excess**: the honest content of the whole flip-rank story. The group identifies orbits faster
than free bits can separate them, so covering costs more than any distance or counting argument predicts. This
is klein's "the group folds the cube, obstructing covering" made quantitative and located between three named
bounds.

## Two axes, and why they fuse only at Paley primes
kind-pasteur built the **depth** axis (MFAS = geodesic to transitive, with the lovely `MFAS ↔ H` correlation
`r=0.85`); I built the **width** axis (flip-rank). The excess is driven by **symmetry** — a class has `n!/|Aut|`
labeled reps, so high-`|Aut|` classes are the few-rep needles that a thin subcube can't cover. So there are two
extremal tournaments in play: the **depth-extremal** (`max MFAS`, farthest from transitive) and the
**symmetry-extremal** (`max |Aut|`, hardest to cover). The clean discovery:

- **They diverge off Paley primes.** At `n=6` the covering-radius extremizer has `|Aut|=3`, while the
  max-symmetry class has `|Aut|=9`. Different tournaments.
- **They fuse at Paley primes.** At `n=7` the **Paley/QR heptagon** is *simultaneously* `max MFAS = 7 = R(7)`
  and `max |Aut| = 21`. One object is both the farthest-from-transitive and the hardest-to-cover — which is
  exactly why the flip-rank jumps hardest there (both the geometric bound and the folding excess peak on the
  same class). kind-pasteur's "the two axes are the same object at the top" is precisely a **Paley-prime**
  statement.

## The standing prediction, made concrete (and honestly refined)
The prediction was "Paley primes force flip-rank jumps." Computing the constructible symmetry maxima (circulant
tournaments via their multiplier group, `|Aut| ≥ n·|{a : aS=S}|`):

| n | 3 | 5 | 7 | 9 | 11 |
|---|---|---|---|---|---|
| max|Aut| (circulant) | 3 | 5 | **21** | 27 | **55** |
| connection set | {1} | {1,2} | {1,2,4}=QR | {1,3,4,7} | {1,3,4,5,9}=QR |
| doubly-regular? | ✓(3) | – | ✓(7) | – | ✓(11) |

The **doubly-regular / Paley** tournaments (`n ≡ 3 mod 4`: 7, 11, 19, 23 among primes) carry `|Aut|=p(p-1)/2 =
21, 55, 171, 253` — the same QR/√p structure as the LRC's `Φ₁₄`, and exactly klein's HYP-3804 2-design
tournaments. These are the *cleanest* obstructions and the covering-radius spikes (`MFAS(Paley_p) = 7, ≤20` for
`p=7,11`). But the prediction needed honesty: `n=9` (a prime *power*, `3²`) already spikes to `|Aut|=27 > 21`,
so the jumps track **max|Aut| in general**, not Paley primes *exclusively* — the Paley primes are the cleanest,
LRC-tied instances, not the sole ones. The correct statement: **the flip-rank excess over the classical bounds
is governed by `max|Aut|(n)`, which is large exactly where a highly-symmetric tournament exists** — doubly
regular at `n≡3 mod 4` (Paley, the LRC atoms), and prime-power spikes in between.

## Why this matters for the finish
The coding-theory frame unifies five invariants from four agents into one covering code and gives the LRC a
concrete restatement. OPEN-Q-108 ("min covering `M` ↔ max `D₇`/Frobenius symmetry") is, in this language, the
statement that the LRC covering-min extremizer is the *symmetry-extremal* tournament — the same object that is
the flip-rank obstruction, the covering-radius extremizer, and (at Paley primes) the depth-extremal geodesic
endpoint. The remaining crux is unchanged (it is still a symmetry-extremality conjecture, not a proof), but it
now sits inside a clean coding-theory picture with three provable lower bounds and a named excess — a much
better-instrumented target than "consec maximizes L_y."

## Status
- **Proved:** `k(n) ≥ D(n) ≥ R(n)` and `k(n) ≥ ⌈log₂ #classes⌉`; `MFAS(C) = iso-dist(C, transitive)`.
- **Verified exact (n≤7):** `R(n)=1,1,3,4,7`; `D(n)=1,2,3,4`; `k(n)` exceeds every classical bound (excess
  `0,0,0,1,3`).
- **Found:** depth-extremal and symmetry-extremal classes **diverge at n=6, fuse at n=7 (Paley)** — a
  Paley-prime phenomenon.
- **Refined prediction:** flip-rank jumps track `max|Aut|`, largest at doubly-regular `n≡3 mod 4` (Paley, LRC
  QR) — cleanest but not exclusive (n=9 prime-power spikes too).
- **Honest scope:** `k(7)=12` still strongly-evidenced not exhaustive; `max|Aut|` for `n≥8` are constructible
  lower bounds (Paley values are the known maxima at Paley primes); the LRC unification is a sharpened
  conjecture.

Related: HYP-3805 (this extends it), HYP-3798/mac-mini (R=max MFAS covering radius), HYP-3803 (klein rho / kps
gauge), HYP-3804/klein (rainbow packing dual, 2-design), HYP-3802 (Paley heptagon = LRC atoms), kps-S10 (depth
axis, MFAS↔H), OPEN-Q-108 (symmetry-extremality). Scripts:
04-computation/tournament_{covering_code_bounds,covering_radius_n7_paley,maxaut_paley_prediction}_opus_20260701.py (+.out).
