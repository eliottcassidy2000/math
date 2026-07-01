---
id: HYP-3801
title: VERBLUNSKY / OPUC ENCODING of the lonely measure (the recursive circle metaphor) + a DICTIONARY of loop-maps that operate group-like + the 3-WAY EXTREMAL-MEASURE EQUIVALENCE. Pushing "runners on a loop" to orthogonal polynomials on the unit circle: the lonely measure mu=(1/L)1_{L_C}dt is a probability measure on the circle, so (Verblunsky) it EQUALS a sequence of Verblunsky coefficients alpha_n∈D built recursively by the Szego recursion. COMPUTED (Levinson/Szego on the moments c_k=hat1(k)/L): (i) alpha_0 = -c_1 -> -cos(2π t*) = -0.887 EXACTLY as r->M_C (the S66 two-atom law is the r->M_C LIMIT, not a global fact -- at small r alpha_0 is tiny from three-gap cancellation among L_C's 28 spread intervals); (ii) as r->M_C the sequence TERMINATES: |alpha_1|=1 (OPUC signature of EXACTLY 2 atoms) at {t*,1-t*}, the Phi6-denominator binding pair; (iii) transition a.c.(r small: alpha_n small, decaying, Szego class) -> atomic(r->M_C: max|alpha_n|->1, E_n->0). DICTIONARY of loop-maps in TWO group families -- ARITHMETIC Z⋉S^1 (rotation R_a:t->t+a [runner flow], multiplication/speed M_k:t->kt [the runners], antipode iota:t->-t [sign], affine A_{k,a}:t->kt+a, phase-residue p=M_n on Z/Phi6) and GEOMETRIC PSU(1,1)=Aut(D) (Blaschke b_a:z->(z-a)/(1-\bar a z), finite Blaschke products deg d, Schur step, Szego map) -- BRIDGED by rotations (in both) and the degree-k notion (M_k arithmetic ~ deg-k Blaschke geometric). 3-WAY EXTREMAL EQUIVALENCE: the construction's extremal lonely measure is [Toeplitz moment matrix rank 2 (mac-mini flat extension HYP-3789/3793)] = [Verblunsky |alpha_1|=1 termination (this)] = [2 atoms at t*,1-t* the Phi6/CF-extremal points (THM-515, S68)] -- three languages for ONE fact. DYNAMICAL reframe: runners are the maps M_v; far-element equidistribution (S65) = MIXING of M_v (large v); L_C = the common M_v-avoidance set
status: MIXED (computation verified + established-OPUC facts + synthesis). VERIFIED (Levinson/Szego on FFT moments, n=14): alpha_0 = -0.018/-0.031/-0.068/-0.204/-0.887 for r=0.05/0.06/0.07/0.074/0.0762 -> -cos(2πt*)=-0.887 at r->M_C; near M_C |alpha|-pattern 0.887,1.000,0.887,1.000 (terminates, 2 atoms); closed-form 2-atom {t*,1-t*} gives alpha_0=-0.8867, |alpha_1|=1.0000 EXACTLY; max|alpha_n| 0.557->0.811->0.987->0.997->1.000 and E_n 1.4e-1->7.4e-14 as r->M_C. The dictionary group structures are standard (Z⋉S^1, PSU(1,1)). The 3-way equivalence is an established OPUC/moment fact (Verblunsky/Geronimus + Curto-Fialkow) applied here. HONEST: a recursive ENCODING + unifying reframe + dictionary, NOT a new proof; the multi-far crux (OPEN-Q-108) is unchanged, now also visible as "no covering measure terminates deeper than the construction's 2-atom".
source: klein-2026-07-01-S69
depends_on:
  - HYP-3800   # S68: phase-residue p(w)=nw mod Phi6, the 2 atoms at the Phi6-denominator points
  - HYP-3790   # S66: the two-atom law (now identified as alpha_0 in the r->M_C limit) [klein signed-correction]
related:
  - HYP-3789   # mac-mini: covering-min = moment relaxation; flat extension = Toeplitz rank = |alpha|=1 termination
  - HYP-3793   # mac-mini/kps: flat-extension moments = Ramanujan sums (the moment side of this)
  - THM-515    # L_C measure = singular series theta-form; the 2 atoms = extremal
  - HYP-3786   # far-element equidistribution = MIXING of the multiplication maps M_v (dictionary/dynamical)
  - HYP-3154   # Joukowski/de Moivre circle-real bridge (the disk<->circle geometry)
  - HYP-3715   # t*=n/Phi6 hexagonal/Eisenstein point (the atom locations)
results:
  - 04-computation/verblunsky_lonely_measure_klein.py
  - 05-knowledge/results/verblunsky_lonely_measure_klein.out
---

# HYP-3801 — Verblunsky/OPUC encoding of the lonely measure + the loop-map dictionary

## The recursive circle metaphor (Verblunsky / OPUC)
The lonely measure `mu = (1/L) 1_{L_C} dt` is a probability measure on the unit circle, so by **Verblunsky's
theorem** it is *equivalent* to a sequence of **Verblunsky coefficients** `alpha_n ∈ D` (the open disk),
built one at a time by the **Szego recursion** `Phi_{n+1}(z) = z Phi_n(z) - conj(alpha_n) Phi_n^*(z)`. This
is a genuine RECURSIVE metaphor for `L_C`: instead of describing the lonely set by its intervals, describe
it by the recursion that generates its orthogonal polynomials.

**Computed** (Levinson/Szego on the moments `c_k = hat1(k)/L` from the FFT of `1_{L_C}`):
- `alpha_0 = -c_1 -> -cos(2π t*) = -0.887` **exactly** as `r -> M_C` (values `-0.018, -0.031, -0.068,
  -0.204, -0.887` for `r = 0.05..0.0762`). So the **S66 two-atom law is the `r -> M_C` LIMIT**, not a
  global fact: at small `r`, `L_C` is 28 spread intervals (three-gap positions) and `c_1` nearly cancels;
  the literal two-atom picture only holds at the binding.
- **Termination**: as `r -> M_C` the sequence collapses to `|alpha_n| = 0.887, 1.000, 0.887, 1.000, …` —
  `|alpha_1| = 1` is the OPUC signature of **exactly 2 atoms**, at `{t*, 1-t*}` (the `Phi6`-denominator
  binding pair). Closed-form check: the 2-atom measure gives `alpha_0 = -cos(2π t*)`, `|alpha_1| = 1`
  exactly.
- **Transition**: `max|alpha_n|` climbs `0.557 -> 0.811 -> 0.987 -> 0.997 -> 1.000` and the Szego error
  `E_n` falls `1.4e-1 -> 7.4e-14` as `r: 0.05 -> M_C`. The measure passes from **absolutely continuous**
  (small, decaying `alpha_n`, Szego class) to **atomic at the binding** (`alpha_n -> 1`, terminating).

## The dictionary of loop-maps (creative functions between points on the loop)
Two group families act on `S^1 = R/Z` (via `z = e^{2πi t}`), bridged by rotations and by the degree-`k`
notion:

**ARITHMETIC family — the affine group `Z ⋉ S^1` (the LRC-native one):**
| map | formula | role | group law |
|---|---|---|---|
| Rotation `R_a` | `t ↦ t+a` | runner flow / time | `(S^1,+)`, abelian |
| Speed `M_k` | `t ↦ k t` (`k∈Z`) | the runners themselves | `M_k M_l = M_{kl}` (monoid `Z`; units `(Z/N)^*`) |
| Antipode `iota` | `t ↦ -t` | complement / sign (S55) | `Z/2` |
| Affine `A_{k,a}` | `t ↦ k t + a` | runner-at-offset | `Z ⋉ S^1` (`M_k R_a M_k^{-1} = R_{ka}`) |
| Phase-residue `p` | `M_n` on `Z/Phi6` | coupling direction (S68) | bijection `(Z/Phi6)` |

**GEOMETRIC family — `Aut(D) = PSU(1,1) ≅ PSL(2,R)` (where OPUC lives):**
| map | formula | role | group law |
|---|---|---|---|
| Blaschke `b_a` | `z ↦ (z-a)/(1-\bar a z)` | disk automorphism | `PSU(1,1)`, non-abelian |
| Blaschke product `B` | `∏ b_{a_j}` | degree-`d` circle self-map | monoid under composition |
| Schur step `S_alpha` | `f ↦ (1/z)(f-alpha)/(1-\bar alpha f)` | the OPUC recursion | iterate → Verblunsky seq |
| Szego map | `Phi ↦ zPhi - \bar alpha Phi^*` | builds orthogonal polys | recursion |

**The bridge.** Rotations `R_a` live in BOTH families. "Speed `k`" has two incarnations: the *arithmetic*
`M_k(t)=kt` (a degree-`k` covering map of the circle) and the *geometric* degree-`k` Blaschke product
`B_k`. Both are degree-`k` self-maps of `S^1`; LRC uses the arithmetic one, OPUC/Verblunsky the geometric
one. This is the dictionary that lets the two theories talk.

## The 3-way extremal-measure equivalence (the synthesis)
The construction's extremal lonely measure (at `r = M_C`) is one object in three languages:
> **[Toeplitz moment matrix has rank 2 / flat extension — mac-mini HYP-3789/3793]** ⟺ **[Verblunsky
> sequence terminates, `|alpha_1| = 1` — this session]** ⟺ **[exactly 2 atoms at `t*, 1-t*`, the
> `Phi6`-denominator / CF-extremal points — THM-515, S68]**.

Curto–Fialkow (flat extension ⟺ finitely-atomic) and Verblunsky/Geronimus (`|alpha_{N-1}|=1` ⟺ `N`
atoms) are the two theorems making this an equivalence, not a coincidence. The **value** `M_C = n/Phi6` is
the depth at which the measure becomes these 2 atoms.

## The dynamical reframe (from the dictionary)
The runners ARE the maps `M_v` (multiplication by the speed). Then:
- `L_C = { t : M_v(t) ∉ danger-arc for all v }` — the common avoidance set of the speed-maps.
- **Far-element equidistribution (S65/HYP-3786) = MIXING of `M_v` for large `v`**: `M_v` is an exact,
  exponentially-mixing endomorphism of the circle (like the `×k` map), so it pushes any measure toward
  Lebesgue — a far runner's danger equidistributes and cannot concentrate on `L_C`. The far-element
  impotence is the *mixing of the multiplication maps in the dictionary*.

## The clean next step (synthesizing)
The multi-far crux (OPEN-Q-108) in OPUC terms: **no covering set's lonely measure terminates deeper than
the construction's 2-atom at `M_C`.** The construction's 2 atoms sit at the CF-extremal point
`t* = n/Phi6 = [0;n-1,n]` — the "hardest" reachable binding. A beater would need a covering measure whose
Verblunsky sequence terminates *later*/at a *deeper* level, i.e., whose atoms sit at an even harder point —
but the covering condition (THM-523) forces a multiple of every `q ≤ n`, pinning the reachable atoms to
the `Phi6` lattice (S68 phase-residue). Making "deepest OPUC termination ⟺ the construction" rigorous is
the same one analytic inequality as S68's difference-phase bound, now phrased as an extremal problem over
Verblunsky sequences of covering measures.

## Honest status
The Verblunsky computation is verified; the dictionary group structures are standard; the 3-way
equivalence is established OPUC/moment theory applied to `L_C`. This is a recursive ENCODING + a unifying
dictionary + a reframe — NOT a new proof. OPEN-Q-108 is unchanged, now also visible as an extremal problem
over the Verblunsky sequences of covering measures.
