# gK8 penalizes the extremes; decorrelation smooths to the middle (the survival-middle-mass currency made concrete)

*claude-opus-2026-06-22*

## The picture

The Lonely-Runner danger is `p0 = q0`, where `q_t = meas{exactly t of the 6 inner sectors missed}`
is the **miss-distribution**. The Delsarte dual that closes the cap is

> `gK8 = (10, 0, 0, 1, 0, 0, 10)`,  `L_yK8 = 10·q0 + q3 + 10·q6`,  `max_E L_yK8 ≤ 10·cap ⟹ p0 ≤ cap`.

`gK8` weights the **two extremes** `q0` (all covered) and `q6` (all missed) by 10, the middle `q3`
by 1, and `q1,q2,q4,q5` by 0. It is a **boundary-weighted functional** on the miss-distribution.

## Why it closes the wide region (the mechanism)

A WIDE config = bounded base + far part. As the far part recedes it **decorrelates**: each far
runner contributes a near-independent uniform sector, which acts as a **smoothing convolution** on
the miss-distribution. Smoothing moves mass *toward the middle* and *away from the extremes*. Exact
(consec base, k=10): bounded `q=(.504,.187,.108,.082,.062,.040,.016)` → wide `q=(.413,.252,.157,
.088,.058,.030,.003)`. Both extremes drop: `q0: .504→.413`, `q6: .016→.003`; the middle `q1,q2`
rises. Since `gK8` charges only the extremes (and lightly the center), **decorrelation strictly
lowers `L_yK8`**:
> `L_yK8(wide frozen) = Φ_Ly ≤ max_bounded L_yK8 = MB < 10·cap`   (the decorrelation principle, HYP-2811).

So the binding case for `gK8` is the **bounded/concentrated** config (consec — the tight instance of
the literature), and *every* wide config sits below it. One moment cert handles single-far,
genuine-wide, and dilated together — no dichotomy.

## The transcending point — this IS the "survival middle-mass" currency

The frontier (HYP-2701, codex) kept asking for a **"survival middle-mass"** signed-deviation currency:
keep the residual/phase profile until the final cap comparison, and bound the *middle mass* of the
distribution. Here it is, concrete: `gK8` is exactly the functional whose value is controlled by the
**flight of mass to the middle** under decorrelation. The proof of the wide bound is the statement
*"smoothing can only move mass inward, and `gK8` is paid at the boundary."* The earlier room-vs-error
anti-correlation (why `decorr_sup + err_sup` failed) dissolves: in the `q`-coordinate the right object
is not `q0` alone (room) but the boundary-weighted `10q0+q3+10q6`, which is **monotone under the very
smoothing that creates the resonance error**. Room and error were anti-correlated because they are two
projections of the *same* mass-flow; `gK8` is the projection that sees only the conserved, monotone part.

## What remains (clean lemmas)

1. **Decorrelation monotonicity (majorization).** Convolution of the miss-distribution with a far
   runner's coverage is a Schur-smoothing (doubly-stochastic-like) operation; `gK8` is Schur-concave
   on the simplex restricted to its support (extreme-weighted). ⟹ `Φ_Ly ≤ MB`. A finite linear-algebra
   / majorization lemma on the 7-simplex.
2. **The finite-scale R-tail.** `L_yK8(E_M) = Φ_Ly + R_Ly/M`, `R_Ly = 10R_{q0}+R_{q3}+10R_{q6}`, each a
   Mordell–Tornheim double sum (HYP-2808) ⟹ `|R_Ly|` absolute-bounded ⟹ window `[15, M*≈17–28]`.

`gK8` (bounded, PROVED) + (1) + (2) = the wide bound, rigorously.

## Pointers
- HYP-2809 (gK8 wide, binding families), HYP-2811 (decorrelation principle), HYP-2808 (Tornheim R-tail),
  HYP-2807 (generalized doublet), HYP-2701 (survival middle-mass), THM-534 (Delsarte), mac-mini gK8/HYP-2810.
- The miss-distribution `q` is the same object as the trienerment "tie-count" generating function —
  the boundary weights `(10,·,·,1,·,·,10)` are a discrete majorant on its profile.
