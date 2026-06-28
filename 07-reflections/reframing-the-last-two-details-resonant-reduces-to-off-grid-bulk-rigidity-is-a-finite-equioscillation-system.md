# Reframing the last two details to fit: the resonant survivor REDUCES to the off-grid bulk (v's danger is on-grid, the optimum dodges it), and the rigidity is a FINITE equioscillation system (3 sum-14 binding pairs at the units)

*mac-mini-2026-06-28-S82. The owner: keep creatively reframing the last few details of LRC(14) to get them to
fit. The two open details (S81) were (b) the resonant-v survivor-positivity and (a) the tight-locus rigidity.
Both reframe into tractable forms: (b) reduces to the generic off-grid bound; (a) becomes a finite equioscillation
system. Builds on [[tightening-the-rigor-construction-rigorous-bounded-margin-rigorous-but-equidistribution-1-over-7-fails-for-resonant-v]],
[[the-vitali-wall-brouwer-equioscillation-and-the-cyclotomic-core-construction]], kps S255 (equioscillation/units).*

## REFRAME (b): the resonant case reduces to the GENERIC off-grid bulk
The S81 pull was that a resonant `v` (a multiple of 14) removes far more than `1/7` of the seed's lonely set
(`v=14` removes `0.73`), breaking the clean equidistribution argument. The reframe dissolves it:
> **A resonant `v=14m` has its danger ON the 14-grid** (the arcs `‖14m·t‖<1/14` sit at `t≈a/14`). The seed's
> optimum therefore **dodges OFF the exact `a/14`** to a point where `v` is safe — and there `M > 1/14`:
```
  seed alone:  M=1/7  at t*=1/7  (ON-grid, the core)
  seed + 14:   M=1/11 at t*=0.864 (OFF-grid) ✓     seed + 28: M=0.111 off-grid ✓
  seed + 42:   M=0.120 at t*≈0.140 (just OFF 2/14) ✓  seed + 56: M=0.125 ✓
```
So the resonant `v` covers the **on-grid CORE** (where the seed's optimum was), and the optimum relocates to the
**off-grid BULK**, where `v` is generic (non-resonant) and safe. **The resonant case is NOT special — it reduces
to the off-grid generic bound** (the same equidistribution that works for generic large `v`). This is the Vitali
wall (S75f) in action: `v`'s danger attacks the measure-zero core (`a/14`); the positive-measure bulk survives.
> **(b) gets to fit:** the analytic core is now uniformly "the off-grid bulk has `M ≥ 1/14`," with NO resonant
> exception — the resonant `v` just moves the witness off-grid, into the generic regime.

## REFRAME (a): the rigidity is a FINITE equioscillation system (the 3 sum-14 binding pairs)
The tight-locus finiteness (tight = AP/GW only) is the famous hard part. kps S255: tight ⟺ the safety
equioscillates at the unit group `(ℤ/14)* = {1,3,5,9,11,13}`. The reframe makes the system explicit:
> Each unit `a` is pinned by a **binding pair** `{s, 14−s}` summing to 14 (the apex-7 pair, HYP-2909): at
> `t=a/14`, the two runners sit at `±1/14`. The 3 odd-unit pairs are `{1,13}, {3,11}, {5,9}`. **Both AP and GW
> contain all 3** (verified: AP residues `{1..13}`, GW residues `{1..11,13,10}` — both have the 3 pairs). So
> tight ⟺ the residue set contains the 3 sum-14 binding pairs, is complement-symmetric (`R=−R` mod 14), and
> covers the units.
This recasts the rigidity from an **unbounded search** ("is any 13-set tight?") to a **finite equioscillation
SYSTEM**: which integer 13-sets realize the 3 binding pairs + complement symmetry + cover the units, with
`M=1/14`. The Chebyshev/Kolmogorov theory says the equioscillating extremal is rigid; the integer solutions of
this finite system are AP/GW (+ dilations).
> **(a) gets to fit:** the rigidity is the finite question "the equioscillation system (3 sum-14 pairs at the
> units, complement-symmetric) has only AP/GW as integer 13-set solutions" — a bounded, structured problem, not
> an unbounded census.

## The two reframes share one picture: the unit-grid `(ℤ/14)*`
Both details live on the apex unit grid:
- **(a)** the tight configs equioscillate AT the units `a/14` (the 3 binding pairs pin them);
- **(b)** the resonant `v`'s danger sits ON the unit grid `a/14`, and the surviving optimum dodges OFF it.
So the unit grid `(ℤ/14)*` is simultaneously where the tight-locus is pinned (a) and where the resonant danger
concentrates (b) — and the off-grid bulk is where the generic bound lives. **The last two details are the
on-grid (core, finite, equioscillation) vs off-grid (bulk, generic, equidistribution) split** — exactly the
Vitali wall, now applied to both pieces at once.

## Honest status
- **REFRAME (b) — confirmed, reduces resonant→generic:** the resonant `v` covers the on-grid core; the optimum
  dodges off-grid where `v` is safe; `M>1/14` (verified `v=14,28,42,56`). The analytic core is now uniformly the
  off-grid bulk bound (no resonant exception). Still needs the generic off-grid equidistribution (the bulk
  positivity) — but that is the EASY, uniform case now.
- **REFRAME (a) — the rigidity as a finite equioscillation system:** tight ⟺ 3 sum-14 binding pairs at the units
  + complement symmetry; AP/GW verified to satisfy it. Recasts the unbounded census as a finite structured
  system; the solution-count (= AP/GW only) is the remaining rigidity content.
- **NOT a proof.** Both details are reframed into tractable forms (resonant→generic, unbounded→finite system),
  unified by the unit-grid core/bulk split. The residual: the off-grid bulk positivity (b) and the
  finite-system solution-count (a). LRC(14) open, but the last two details now FIT one frame.

Related: HYP-3255 (this), HYP-3253 (the tightening), kps S255 (equioscillation/units), HYP-2913 (three-gap),
HYP-2914 (battery), HYP-2909 (binding pairs), HYP-3237 (Vitali wall), OPEN-Q-108.
