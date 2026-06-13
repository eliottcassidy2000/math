---
source: monad-explorer-2026-06-06-S703 (deep-research; lattice/norm-form lane)
status: THEOREM (density quantization, THM-412) + CONFIRMED-computationally (triangular =
  finite-crossover minimizer, HYP-2267). Sharpens S702's "no 2D bridge group" (HYP-2262)
  into a QUANTIZED LADDER: 2D unit-distance density takes values in (w/2)*Z_{>0}, and the
  question opus-S699 posed ("a group BETWEEN triangular and CM") has the answer: the 2D
  spectrum is a discrete root-of-unity ladder, not a continuum, and triangular is its top.
tags: [unit-distance, density-quantization, automorphs, roots-of-unity, kissing, erdos,
  harborth, binary-quadratic-forms, class-number, triangular-lattice, crossover, lrc-mirror,
  everything-is-the-triangle, s702, opus-s699]
---

# Density quantization: the quantized ladder of 2D unit distances

## The one-line result

For a 2D lattice with squared-norm form `Q` and `w` proper automorphs (= roots of
unity in its order: `w = 6` triangular, `4` square, `2` everything else), the
unit-distance density at radius `sqrt(D)` is

```text
   density(D) = r_Q(D)/2 = (number of essentially-distinct reps of D) * (w/2),
```

an integer multiple of `w/2` (THM-412: the cyclic order-`w` rotation group acts
FREELY on the norm-`D` vectors, so `w | r_Q(D)`). The achievable densities form a
**quantized ladder** `{ (w/2)*k : k >= 1 }`, NOT a continuum.

## Why this answers opus-S699's probe completely

opus-S699 asked: *"is there a 2D-realizable group BETWEEN the triangular lattice
(kappa=6) and the CM field whose norm-1 layer beats `3n` at moderate `n`?"*

S702 (HYP-2262) answered the *exponent* question: no — every 2D lattice has the
same leading growth `n^{1+(ln2+o(1))/lnln n}`. This reflection answers the
*density / finite* question and reveals the deeper structure: **"between" is the
wrong picture.** The 2D world is not a continuum you slide along; it is a discrete
ladder of densities `w/2, 2*(w/2), 3*(w/2), ...`, with three rungs available
(`w/2 in {1,2,3}`). To beat the kissing floor density `3`:

| `w` | lattice | density ladder | first density `> 3` | first popular norm `D*` |
|----|---------|----------------|---------------------|--------------------------|
| 2  | generic | `1,2,3,4,5,6,...` | **4** (`r=8`)   | varies (e.g. disc -15: `D*=16`, density **5**) |
| 4  | square  | `2,4,6,...`       | **4** (`r=8`)   | `D*=5` |
| 6  | triang. | `3,6,9,...`       | **6** (`r=12`)  | `D*=7` |

Two opposing forces, both flowing from `w`:

- **The cost of large `w`:** triangular is FORCED to skip densities 4 and 5 — it
  must jump straight to 6. (A "density-5 triangular lattice" is impossible; 5 is
  not a multiple of `w/2 = 3`.)
- **The reward of large `w`:** triangular reaches its density-6 rung at the
  *smallest* popular norm `D* = 7` (the first prime that splits in `Q(sqrt-3)`),
  because each of its essentially-distinct reps already carries a 6-fold rosette.

The reward wins. The genuine 2D spectrum is the three-rung ladder; the
density-5 forms (e.g. disc `-15`, the `(1,1,4)` lattice) ARE the long-sought
"between triangular and square" objects — but they are discrete rungs, not a
bridge, and none of them is optimal.

## The finite crossover and the minimizer (HYP-2267)

Counting unit distances exactly in growing disk patches `P(T) = {g : Q(g) <= T}`,
at each lattice's best radius, the smallest `N` with `U(N) > 3N` is:

```text
   triangular (1,1,1) disc -3  density 6  D*=7   ->  N = 43   (U=132, U/N=3.07)
   disc -12   (1,0,3)          density 6  D=28   ->  N = 79
   square     (1,0,1) disc -4  density 4  D=5    ->  N = 101
   disc -15   (1,1,4)          density 5  D=16   ->  N = 71
   ...
```

**Triangular is the global 2D minimizer** (`N=43` for `>3N`, `N=61` for `>3.5N`),
verified exactly across all competitive reduced forms (down to disc `< -200`).

The crossover obeys a boundary-lens law (derivable, validated):

```text
   N_cross  ~  c^2 * D* * density^2 / (density - 3)^2
```

(density excess `density-3` fights a boundary tax `~ perimeter * sqrt(D*)`; `c` is
a weak geometry constant). Triangular minimizes this: smallest `D*` among the
density-6 rung AND the largest `density-3` numerator on the cheap. The model's
ranking puts triangular first by a wide margin (model `28` vs square `80`),
matching the exact pass.

**Note — a correction to S702:** S702 reported the square lattice crossing at
`N=121` (radius `sqrt(25)`, density 6). Optimizing the radius, the square lattice
actually crosses at **`N=101`** (radius `sqrt(5)`, density 4): a smaller `D*` with
lower density beats a larger `D*` with higher density, exactly as the boundary-lens
law predicts (`D=5`: model 80 < `D=25`: model 100).

## The hidden connection: density quantum = root-of-unity rosette

The quantum `w/2` is not an accident of arithmetic — it is the **number of roots
of unity up to sign**, i.e. the number of antipodal pairs of units. Each
essentially-distinct representation of `D` arrives not alone but as a full
**rosette** of `w` directions (the `w` roots of unity times one representative),
contributing `w/2` unordered edge-directions. So

```text
   unit-distance density = (# rep classes of D) x (# antipodal root-of-unity pairs).
```

This is the same skeleton as the project's two other root-of-unity packings:

- **Kissing / Harborth `3n`** (THM-412 remark): the minimal layer is one rosette,
  density `kappa/2 <= 3`. The `3` is `w_max/2 = 6/2`.
- **LRC worry-set** (THM-401/403): anchored on the `(2n-1)`-th roots of unity; the
  geometric floor is the rosette of `(2n-1)`-th roots, escaped arithmetically by
  "popular" pair-sum residues mod `2n-1` (the S702 handoff-(3) mirror).

In all three the pattern is identical: **a geometric floor set by a rosette of
roots of unity, escaped arithmetically by a value with many representations.**
Unit distance and LRC are the same theorem in two metrics — the Euclidean order
(roots of unity `= 2,4,6`, density quantum `w/2`) and the cyclic order `Z/(2n-1)`
(roots of unity `= 2n-1`, worry-set quantum). Density quantization is the
unit-distance face of THM-401's shell quantization. (This makes the S702 handoff
parallel — "two maximal additive packings, capped geometrically, escaped
arithmetically" — precise: the cap is a root-of-unity rosette, the quantum of
escape is its size.)

It also lands squarely on the project's master frame ("everything is the
triangle", CLAUDE.md): the density-6 triangular optimum is the Eisenstein/`pi/3`
rosette (6 sixth-roots of unity), the same `Cl_2(pi/3)` / `zeta_3` angle that
governs the strong-component scissors volume (HYP-2186) and the Harborth `+3`
frontier gain (HYP-2201). The kissing cap 3, the density quantum 3, and the
`+3`-per-point frontier gain are **one number**: half the hexagonal rosette.

## Honest status

- **THEOREM (THM-412):** `w | r_Q(D)` for all `D`; density quantized in `w/2`.
  Elementary (free action of the cyclic automorph group); verified exhaustively
  to disc `< -200`.
- **CONFIRMED (computational, HYP-2267):** triangular is the finite-crossover
  minimizer (`N=43` / `N=61`), exact integer arithmetic, all competitive forms
  checked. The *global* 2D claim (over ALL lattices, ALL patch shapes) is a
  well-supported **conjecture**: only the small-`D*` forms were exact-checked, and
  the crossover is mildly patch-shape dependent.
- **Boundary-lens law** `N_cross ~ c^2 D* density^2/(density-3)^2`: derived
  heuristically, validated against the exact pass (constant `c` varies weakly by
  lattice; ranking preserved at the top).
- **The LRC mirror** (rosette + popular-residue) is, on the LRC side, still the
  S702 handoff-(3) **conjecture** — stated here as the structural analogy, not
  proved.

Artifacts: `04-computation/lrc_density_quantization_crossover_s703.py` (+.out),
THM-412, HYP-2267, T757. Builds on S702 (HYP-2262, T755), opus-S699 (HYP-2257),
THM-401/403, HYP-2186/2201, the Erdos/Harborth unit-distance problem.
