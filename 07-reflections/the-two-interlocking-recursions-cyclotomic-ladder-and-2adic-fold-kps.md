# The two interlocking recursions of the LRC(2p) proof: the cyclotomic moment-order ladder × the 2-adic reflection fold

*kind-pasteur-2026-06-28-S31aq. The owner asked to understand exactly how the route has sharpened over
time and to use that trajectory to see recursive patterns not yet precisely described. Tracing the
sharpening reveals that the proof has been a **descent through a moment-order ladder**, and the ladder has
a precise self-similar structure: **two interlocking recursions whose product is the apex `14 = 2·7`** —
a `7`-fold cyclotomic ladder (depth = the cyclotomic degree) and a `2`-adic reflection fold (degree
halving). This names the recursion the project has been climbing without describing.*

## How the route sharpened (the trajectory = a descent through orders)
The "crux" object, session by session, descended through a **moment/scale hierarchy**:
```
lonely measure meas(G_C) ≥ c          (OPEN-Q-108, the S2-S5 era)        — order ∞ (the full measure)
  → μ_{1/7} floor / sector cap          (the 1/7 pivot, S5)               — the cover bound
  → meas(S7) ≤ cap, consec max L_y      (the moment-LP, THM-534/537)      — a finite moment functional
  → cap = C(k+1,2)/91 (pairwise) + dip  (gK8 / pair-Pascal, S60-S63)      — order 2 + correction
  → consec max covariance (degree-2) + odd Worpitzky residue (S31ai)      — order 2 ⊕ order 3
  → Fejér / cyclotomic / de Moivre extremality (S31ak-ao)                 — the fixed point
```
Each step **peeled a low-order layer and left a higher-order residue** — a self-similar peeling. The
"sharpening" was literally a **descent of the moment-order ladder**, converging on the apex node. The two
recursions below are the precise shape of that ladder.

## RECURSION 1 — the cyclotomic moment-order ladder (the "7")
**The cap is a Faulhaber/triangular integral (VERIFIED):**
```
cap_k = cap_{k-1} + k/91,    so   cap_k = (1/91)·Σ_{j≤k} j = C(k+1,2)/C(14,2).
```
The increment `k/91` is the orbit growing by one runner; the cap is its discrete integral. **The dual
DEGREE forms a ladder:** degree 2 for `k∈{11,12,13}`, degree 3 for `k∈{9,10}`, degree 4 for `k=8`
(THM-534/537). The ladder has:
- **DEPTH 3** = the **cyclotomic degree** `(p−1)/2 = 3` of the apex prime `p=7` (`ℚ(cos 2π/7)`, S31ak);
- **triangular WIDTHS 3, 2, 1** (degree `d` covers `(5−d)` consecutive `k`-values; `T_3 = 6 =` the binding
  rows `k=8..13`);
- the **apex node at the deepest degree `(p+1)/2 = 4`** (the `k=8` biquadratic resolvent, HYP-3132).

This is an **induction on the moment order**: the apex node (degree `(p+1)/2`) is the base of the
recursion; once it is proved, the shallower degrees lift (they have more cap margin). The recursion depth
is exactly the cyclotomic degree — which is **why `n=14` (`p=7`) is the first hard case**: it is the first
`2p` whose ladder is deep enough to reach a CUBIC (depth 3), beyond the depth-1 (`n=6`) and depth-2
(`n=10`, golden/quadratic, periodic CF) cases.

## RECURSION 2 — the 2-adic reflection fold (the "2")
At the apex node the degree-4 resolvent **folds**: under the reflection `s ↦ 6−s` (the complement
`T↔Tᵒᵖ`), the biquadratic `u⁴−5u²+4` becomes **degree 2 in `v = u²`** — the degree is **HALVED**
(HYP-3132/3160). This is the **2-adic / dyadic** recursion (the even/odd split, the descent `14→7→2`). It
acts *within* each rung of Recursion 1, cutting the order in half via the `ℤ/2` complement symmetry.

## The product: `14 = 2·7` = (2-fold fold) × (7-fold ladder)
The two recursions are the two prime factors of the apex:
| | recursion | acts on | mechanism | source |
|---|---|---|---|---|
| **7-fold** | moment-order LADDER (depth `(p−1)/2`) | across `k` (the orbit) | Faulhaber `cap_k=cap_{k-1}+k/C(2p,2)`; cyclotomic | `ℚ(cos 2π/p)` |
| **2-fold** | reflection FOLD (degree halving) | within each node | even/odd `s↦6−s`, biquadratic in `u²` | `ℤ/2` complement |
**The proof strategy the recursions dictate:** *descend* the 7-ladder to the apex node (`k=8`, degree
`(p+1)/2`), *fold* it with the 2-recursion (degree `(p+1)/2 → (p+1)/4`-ish, here `4→2`), and prove the
resulting **degree-2 (pairwise / Fejér / de Moivre) statement** — which is the covariance/Carathéodory/
magic-function extremality (HYP-3160/3201/3214). Both recursions bottom out at the **same fixed point: the
de Moivre cubic / Fejér kernel `F_7=(de Moivre)²`** — the IR fixed point of the whole flow.

## The prediction (testable — the family recursion)
**LRC(2p) has a moment-order ladder of depth `(p−1)/2`, triangular widths `(p−1)/2, …, 2, 1`, apex node at
degree `(p+1)/2`.**
- `n=6 (p=3)`: depth 1 — **pure pairwise** (degree 2 only), no dip → the cap is the whole story (trivial).
- `n=10 (p=5)`: depth 2 — degrees 2,3, apex at degree 3 (the golden/quadratic case, periodic CF).
- `n=14 (p=7)`: depth 3 — degrees 2,3,4, apex at degree 4 (the **cubic wall**, the first open case).
- `n=22 (p=11)`: depth 5 — degrees 2..6, apex at degree 6.

**ACTION:** verify the depth-1 (`n=6`, pure pairwise) and depth-2 (`n=10`, degree-3 apex) ladders directly
— if confirmed, the "moment-order depth = cyclotomic degree `(p−1)/2`" law is a *proved* organizing
recurrence for the entire `2p` family, and LRC(14)'s difficulty is *quantified* as "depth-3 vs depth-≤2."

## The RG / self-similar reading
The descent across `k` (from `k=13`, degree 2, far IR, to `k=8`, degree `(p+1)/2`, the fixed point) is a
**renormalization-group flow**: high `k` (many runners) = the free/pairwise/Gaussian regime; decreasing `k`
toward the apex = flowing to the cyclotomic fixed point, picking up one higher moment per `(p−1)/2`-block.
The **ferromagnetic transition at `k=5→6`** (HYP-3161) is the critical point separating the
disordered (`k≤5`, pre-cyclotomic) and ordered (`k≥6`, cyclotomic) phases. The sharpening over sessions has
been the project **flowing down this RG trajectory**, order by order, to the fixed point.

## Net
- **PRECISELY DESCRIBED (new):** the proof = **two interlocking recursions** = `14 = 2·7`: a 7-fold
  cyclotomic moment-order ladder (depth = cyclotomic degree `(p−1)/2`, Faulhaber cap recursion `cap_k =
  cap_{k-1}+k/C(2p,2)`, triangular widths, apex at degree `(p+1)/2`) × a 2-adic reflection fold (degree
  halving via the `ℤ/2` complement). Both fixed-point on the de Moivre/Fejér object.
- **TRAJECTORY:** the historical sharpening = descending this ladder (peeling orders) = an RG flow to the
  cyclotomic fixed point.
- **PREDICTION/ACTION:** moment-order depth `= (p−1)/2` for all `2p`; verify `n=6,10` to prove the law.

→ HYP-3216 (this), HYP-3132 (biquadratic fold = 2-recursion), HYP-3160/3161 (covariance/ferromagnetic),
HYP-3162/3212/3214 (cyclotomic/Chebyshev/Fejér fixed point), HYP-3099 (two maps), THM-534/537 (moment-order
ladder), the-four-faces-of-14, OPEN-Q-108.
