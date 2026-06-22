# The ELEMENTARY closure of the LRC(14) witness floor — and its UNIFICATION with the p0 wide bound

**kind-pasteur-2026-06-22-S30** · feeds THM-527 (Part G crux), the `hfloor` Lean node, OPEN-Q-108.

## Summary

The remaining analytic crux of the LRC(14) witness route — a **uniform positive
floor** `rho*_glob(P,E) > 0` for the global 1/7-witness density — is closed by an
**elementary Bonferroni bound**, and is moreover **IMPLIED by the team's already-closed
p0 wide bound**. This *dissolves the compactness obstruction* that THM-527 Part G
called "the genuine remaining crux": no compact shape space, no continuity/closure
argument, no `Vmax ≤ V0` finite check, and (via the unification) **no new analytic
lemma** is needed for positivity.

## ★ THE UNIFICATION (headline): the witness floor follows from the p0 wide bound

`GOOD(E)` is the complement of the 1/7-dense set, so `meas(GOOD) = 1 − D(E)`,
`D(E) = meas{cluster phases 1/7-dense}`. Bonferroni gives
`rho*_glob = meas(GOOD ∩ G_P) ≥ meas(GOOD) + meas(G_P) − 1 = meas(G_P) − D(E)`.

**Elementary lemma (PROVED, not just a.e.):** `D(E) ≤ p0(E)`, where
`p0(E) = meas{x : all 6 inner sectors [j/7,(j+1)/7), j=1..6, are hit by some frac(e x)}`
is *the repo's own cover atom* `measS7`. Proof: if inner sector `S_j` holds no phase,
the adjacent gap is `≥ 1/7`; if it is `= 1/7` exactly, a phase sits at `j/7 ∈ S_j`
(half-open sector), contradiction. So `1/7-dense ⟹ all inner sectors hit`, i.e.
`{dense} ⊆ {S7-cover}` as sets. Hence

```
   rho*_glob(P,E)  ≥  meas(G_P) − D(E)  ≥  meas(G_P) − p0(E).                 (UNIF)
```

Now the **DUALITY** (verified exactly, `lrc_capGP_exact_kps.py`):
`cap_k = min_{|P|=13−k} meas(G_P)` equals the p0 `Q`-plateau `capRat(k)` for every
`k = 8..13`. The team's **wide bound** `p0(E) ≤ cap_k` (closed exhaustively: Leg-C
HYP-2817, gK8, single-far THM-563) then gives, with `meas(G_P) ≥ cap_k`,

```
   rho*_glob(P,E)  ≥  meas(G_P) − p0(E)  ≥  cap_k − max_E p0(E)  =  δ_k  >  0,
```

where `δ_k = cap_k − max p0` is **exactly the team's wide-bound margin** (≈ 0.05–0.16).
**So the 1/7-witness route and the p0 sector route UNIFY:** `rho*_glob > 0` is a
*corollary* of `p0 ≤ cap`, via the elementary `D ≤ p0`. The duality `cap_k = min meas(G_P)`
is the linchpin — it is *why* `p0(E) ≤ cap_k = min meas(G_P) ≤ meas(G_P)`. Closing the
p0 wide bound closes the witness floor too, with no extra analytic content.

VERIFIED (`lrc_witness_p0_unification_kps.py`): `D ≤ p0` and `rho*_glob ≥ meas(G_P)−p0`
on all test `(P,E)` (worst margin `+0.054` at k=8 worst-cap `P=(1,5,7,8,9)`).

## The Bonferroni floor (sharper, alternative — needs Lemma A)

The Bonferroni bound also yields a *larger* explicit floor via the pure witness
measure `ν(E) = meas(GOOD(E))` (no small part), at the cost of one three-distance
lemma. This route is recorded below; the UNIFICATION above is the rigorous closure.

## Setup (the corrected global-witness object)

After KPS Thread A, the live witness object is the **global 1/7 density**
(the via-`Vmax` 2/7 density `rho*_crit` was the wrong object — it hits exact 0 on
admissible covering sets). For a covering 13-set `S = P ∪ L`, `P = S∩{1,…,13}`,
cluster co-offsets `E` (`|E| = k`, `0 ∈ E`), and `G_P = {x : ‖p x‖ ≥ 1/14 ∀ p∈P}`:

```
rho*_glob(P,E) = witnessG2(P,E) = meas{ x ∈ G_P : circular maxgap{frac(e x): e∈E} > 1/7 }.
```

`rho*_glob > 0 ⟹ M(S) ≥ 1/14` (THM-527 Part A, the slow-fast witness implication;
a separately-owned node). **This note proves `rho*_glob > 0` — the floor.**

## The Bonferroni reduction (elementary)

Write `GOOD(E) = {x : maxgap{frac(e x)} > 1/7}` and `ν(E) = meas(GOOD(E))`. Then
`rho*_glob = meas(GOOD(E) ∩ G_P)`, and the **Bonferroni / union bound** gives

```
   rho*_glob(P,E)  =  meas(GOOD(E) ∩ G_P)  ≥  ν(E) + meas(G_P) − 1.        (★)
```

So `rho*_glob > 0` follows from two **universal, decoupled** lower bounds:

- **Lemma A** (pure three-distance, no small part): `ν(E) ≥ ν_consec(k)` for every
  `k`-shape `E`.
- **Lemma B** (proven canon rationals): `meas(G_P) ≥ cap_k := min_{|P|=13−k} meas(G_P)`.

and the finite arithmetic fact `ν_consec(k) + cap_k − 1 > 0` for all `k`.

## The finite arithmetic core (EXACT, verified)

`ν_consec(k) = meas{x : the Steinhaus orbit {0,x,…,(k−1)x} has maxgap > 1/7}`,
an exact rational (three-distance: good `x` live only near `a/b`, `b ≤ 6`).
`cap_k = min_{|P|=13−k} meas(G_P)`, an **exact finite rational computation**
(`lrc_capGP_exact_kps.py`): the minimizer is always `P = {1} ∪ {largest speeds}`
— `cap_8 @P=(1,5,7,8,9)`, `cap_9 @(1,11,12,13)`, `cap_10 @(1,12,13)`,
`cap_11 @(1,13)`, `cap_12 @(1)`, `cap_13 @∅`. **DUALITY (verified):**
`min_{|P|=13−k} meas(G_P)` equals the p0 `Q`-plateau `capRat(k)` *exactly* for every
`k=8..13` — the small-part safe-measure floor IS the p0 cap (same rational, both routes).

| k  | ν_consec(k)        | ≈        | cap_k        | ≈        | floor = ν+cap−1 | ≈        |
|----|--------------------|----------|--------------|----------|-----------------|----------|
| 8  | 691/735            | 0.94014  | 2243/5880    | 0.38146  | **+**           | +0.32160 |
| 9  | 247/294            | 0.84014  | 1979/4004    | 0.49426  | **+**           | +0.33440 |
| 10 | 38/49              | 0.77551  | 55/91        | 0.60440  | **+**           | +0.37991 |
| 11 | 1381/2205          | 0.62630  | 66/91        | 0.72527  | **+**           | +0.35157 |
| 12 | 13823/24255        | 0.56990  | 6/7          | 0.85714  | **+**           | +0.42704 |
| 13 | 477/1078           | 0.44249  | 1            | 1.00000  | **+**           | +0.44249 |

**All k=8..13: floor `ν_consec(k) + cap_k − 1 > 0`, worst = `1891/5880 ≈ 0.32160` at k=8 (exact rational).**
For `k ≤ 7`: `ν(E) = 1` for *every* shape (pigeonhole: `k` points ⟹ maxgap ≥ 1/k ≥ 1/7),
so `rho*_glob = meas(G_P) = cap_k > 0` trivially.

Consistency: the Bonferroni floor +0.3216 (k=8) sits just below the actual
exhaustive min `G2 = 8152/24255 ≈ 0.336` (claude-sonnet-S5) — (★) is tight and
correct. The floor +0.32 is also **far above** the conservative admissible floor
`m_P = 14249/252252 ≈ 0.0565` (THM-530) used in the Lean skeleton, so this route
proves the existing `hfloor : witnessMP ≤ witnessG2` obligation a fortiori.

## Lemma A: ν(E) ≥ ν_consec(k) — status and the closure architecture

**Key structural fact:** `ν` is **SCALE-INVARIANT**, `ν(cE) = ν(E)` for every
integer `c ≥ 1` (the map `x ↦ cx mod 1` is measure-preserving and the integrand is
1-periodic). So `ν` depends only on the *primitive* shape `E/gcd(E)`; "wide via
scaling" does **not** lower `ν`. The only candidate minimizers are primitive
shapes, and the claim is that the **most coherent** one (consecutive = a single
arithmetic progression of frequencies) minimizes `ν`.

Lemma A splits into the SAME architecture as the already-closed single-far
(THM-563, mac-mini) and genuine-wide Leg-C (HYP-2817, claude-opus) legs —
**[finite core] + [decorrelation tail]** — but here the tail is the **SAFE
direction with enormous margin**:

- **CORE (bounded spread `≤ W*`):** EXACT-rational exhaustive. Verified: consecutive
  is the strict ν-minimizer over all primitive shapes with spread `≤ k+4`
  (k=8..12, this run) and, in a separate wide-stress run, over spread `≤ 80` with
  one or two far elements (k=8,9). 0 counterexamples.
- **TAIL (wide spread `> W*`):** as the primitive shape widens, `{frac(e x)}`
  **decorrelate** toward i.i.d. uniform, whose maxgap concentrates near
  `H_k/k ≈ 0.34 ≫ 1/7`. Hence the dense-measure `D(E) = 1 − ν(E) → 0`, i.e.
  `ν(E) → 1`. We need only the **loose** bound `ν(E) > 1 − cap_k` (i.e. `D < cap_k`,
  with `cap_k ≥ 0.38`), while the truth is `D ≈ 0` — a `~0.38` margin. A crude
  Erdős–Turán / Weyl discrepancy bound on `D(E)` therefore closes the tail; this
  is the *opposite* regime to the razor-thin p0/2-7 wide bound (margin 0.13).
  **Float scan confirms** (N=200k): at k=8 the min ν over random primitive shapes
  *rises* with spread — `0.980 (≤11) → 0.985 (≤16) → 0.994 (≤32) → 1.000 (≤64)`,
  all `≫` the threshold `1−cap_8 = 0.6185`; same at k=9 (`0.940 → 0.998`). So
  consecutive (ν=0.940) is a *sharp isolated* minimum; every wider shape is `≥ 0.94`
  and typically `≈ 0.98–1.0`. The tail margin is `~0.3–0.4`, not razor-thin.

**This is why the 1/7 witness route is the right one:** in the p0/2-7 route wide
configs are *dangerous* (tight); in the ν/1-7 witness route wide configs are
*comfortable* (ν → 1). The obstruction sign flips.

## Net

`rho*_glob > 0` (THM-527 Part G, the witness floor / OPEN-Q-108 crux) reduces to:
**(★) Bonferroni [elementary] + Lemma B [proven rationals] + the finite floor table
[exact] + Lemma A [scale-invariant three-distance, CORE exhaustive + decorrelation
TAIL].** The compactness obstruction is dissolved by the Bonferroni split into a
scale-invariant pure-shape lemma and the proven small-part floor. The one remaining
rigorous step is Lemma A's decorrelation tail bound — a loose, safe-direction
Weyl/ET estimate.

Scripts: `04-computation/lrc_nu_universal_minimizer_kps.py` (wide stress, k=8,9 to
spread 80), `04-computation/lrc_nu_floor_and_tail_kps.py` (exact core k=8..13 +
float tail scan). Outputs in `05-knowledge/results/`.
External: Lonely Runner Conjecture; Steinhaus three-gap; Weyl equidistribution;
Erdős–Turán–Koksma discrepancy.
