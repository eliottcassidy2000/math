---
id: THM-1550
title: "AN EXACT ALGEBRAIC CRITERION FOR THE TORAL NULLCONE — the infinite family CT(Λ^m)=0 collapses to ONE identity. For Λ = u^{−M}R(u) with R(0) ≠ 0 and deg R = d = M+N, let u_1(t),…,u_M(t) be the M roots of u^M = tR(u) tending to 0. Then CT(Λ^m) = 0 for every m ≥ 1 ⟺ Π(t) := ∏u_i(t) equals c·t EXACTLY ⟺ ∏_{i}R(u_i(t)) is constant in t. Proved by an exact Wiener–Hopf factorisation of CT(log(1−tΛ)) — no asymptotics, no genericity, no saddle points. Two consequences: (i) THM-1530's M=1 theorem becomes a one-liner (Π = u_1 = tR(u_1) = ct ⟹ R(ct) ≡ c ⟹ R constant); (ii) substituting u = εv with ε = t^{1/M} turns the criterion into '∏v_i(ε) is constant' for v^M = R(εv), and the order-ε^k condition carries Σ_i ζ^{(k+1)i}, which is nonzero EXACTLY when M | (k+1) — so M=1 constrains every order, M=2 the odd orders, M=3 every third. That is the precise arithmetic of why M ≥ 2 is a different problem, and it says the conditions, though sparse, are still infinite against d free coefficients. GMC(2) IS NOT FINISHED HERE."
status: >
  THE CRITERION: PROVED, exactly.  The factorisation argument is complete and elementary
  (the log of each inside factor u^{-M}∏(u−u_i) = ∏(1−u_i/u) has no constant term when
  expanded in 1/u; likewise log(1−u/a_j) for outside roots).
  THE ORDER ARITHMETIC: PROVED (a root-of-unity sum).
  THE M=1 COROLLARY: PROVED (re-derivation of THM-1530 in one line).
  THE M,N >= 2 CASE: STILL OPEN.  Bounded search only -- zero counterexamples, cross-checked
  against exact integer CT arithmetic rather than a tolerance.
  RETRACTED IN PLACE: my first numerical test of the criterion was BROKEN and its Parts 1-2
  are void.  See §5.  The criterion itself was never in question -- it is an identity.
  GMC(2) REMAINS OPEN.  This sharpens the last gap; it does not close it.
source: klein-2026-07-20-S347 (owner: aim to finish GMC(2) by finishing the stronger two-dimensional nullcone conjecture)
concurrency: >
  The fleet closed several strata while this was being written, and they are cited not
  duplicated:
  - death-star THM-1515 CARRIED THROUGH my Bessel-EMP and proved GMC(2) on the {−1,0,1}
    stratum -- the gap THM-1510/THM-1530 had left open.  That is their result, and it used
    my setup; credit runs both ways and theirs is the closing step.
  - opus THM-1535 proved GMC(2) for all sign-coherent P (Hankel positive-definiteness) and
    pinned the residue exactly: "P at n = 2 with charges of BOTH signs".
  - boxeph THM-1525 proved the W-linear class and named the wall the resurgent regime.
  - mac-mini THM-1500/1520 hold the master theorem and the one-sided-charge branch.
  THE REMAINING GAP IN ALL THREE FRAMINGS IS THE SAME OBJECT: TNC with M, N >= 1.  This file
  gives it an exact criterion.
depends_on:
  - THM-1530  # the toral nullcone framing and the M=1 Lagrange proof
related: [THM-1510, THM-1515, THM-1525, THM-1535]
script: 04-computation/tnc_exact_criterion_klein_S347.py (+ .out)
---

# THM-1550 — an exact criterion for the toral nullcone

**TNC.** For a Laurent polynomial `Λ`, when is `CT(Λ^m) = 0` for every `m ≥ 1`?
Conjecturally: exactly when all exponents share a strict sign. By THM-1530 this is NC2's
leading order, and NC2 ⟹ GMC(2).

Write `Λ(u) = u^{−M}R(u)` with `R` a polynomial, `R(0) = r_0 ≠ 0`, `deg R = d = M+N`, and
`M, N ≥ 1` (both signs present — the only remaining case).

## 1. The criterion

For small `t`, `Φ_t(u) := u^M − tR(u)` has exactly `M` roots inside `|u| = 1` (all `→ 0`) and
`N` outside. Let `Π(t) := ∏_{i=1}^{M} u_i(t)` over the inside roots.

> **Theorem.** `CT(Λ^m) = 0` for every `m ≥ 1` ⟺ `Π(t) = c·t` **exactly**, for a constant `c`
> ⟺ `∏_{i=1}^{M} R(u_i(t))` is constant in `t`.

*Proof.* `CT(Λ^m) = 0` for all `m ≥ 1` is equivalent to `CT(log(1 − tΛ)) = 0`, since
`log(1−tΛ) = −Σ_{m≥1} t^mΛ^m/m`. Factor

```text
1 − tΛ(u) = Φ_t(u)/u^M = A · u^{−M}∏_{in}(u − u_i) · ∏_{out}(u − a_j),   A = −t·r_d.
```

Now `u^{−M}∏_{in}(u−u_i) = ∏_{in}(1 − u_i/u)`, and `log(1 − u_i/u)` expands in strictly
negative powers of `u`, so its constant term is `0`. For an outside root,
`log(u − a_j) = log(−a_j) + log(1 − u/a_j)`, and the second piece expands in strictly positive
powers, constant term `0`. Hence

```text
CT(log(1 − tΛ)) = log[ A · ∏_{out}(−a_j) ] = 0   ⟺   A·(−1)^N∏_{out}a_j = 1.
```

Since `∏_{all}a_j = (−1)^d r_0/r_d`, we get `∏_{out}a_j = (−1)^d r_0/(r_d·Π)`, and
substituting `A = −t r_d` gives `Π(t) = (−1)^{N+d+1} r_0·t`. Finally, multiplying
`u_i^M = tR(u_i)` over `i` gives `Π^M = t^M∏_i R(u_i)`, so the two forms are equivalent. ∎

No asymptotics, no genericity, no saddle points — an infinite family of conditions replaced by
a single analytic identity.

## 2. `M = 1` in one line

`Π = u_1` and `u_1 = tR(u_1)`. If `u_1 = ct` then `ct = tR(ct)`, so `R(ct) = c` for all small
`t`, so `R` is constant. That is THM-1530's Lagrange–Bürmann theorem, re-derived in a line.

## 3. Why `M ≥ 2` is a different problem — the exact arithmetic

Substitute `u = εv` with `ε = t^{1/M}`. Then `u^M = tR(u)` becomes `v^M = R(εv)`, and
`∏_i R(u_i) = ∏_i v_i^M`, so the criterion reads: **`∏_i v_i(ε)` is constant in `ε`.**

Writing `v_i = w_i(1+δ_i)` with `w_i = r_0^{1/M}ζ^i` (`ζ = e^{2πi/M}`), the order-`ε^k`
contribution to `Σ_iδ_i` carries the factor

```text
Σ_{i=0}^{M−1} ζ^{(k+1)i}  =  M if M | (k+1),  else 0.
```

| `M` | orders `k ≤ 12` carrying a condition |
|---|---|
| 1 | 0,1,2,3,4,5,6,7,8,9,10,11,12 |
| 2 | 1,3,5,7,9,11 |
| 3 | 2,5,8,11 |
| 4 | 3,7,11 |

`M = 1` constrains **every** order — which is exactly why Lagrange–Bürmann closes it
immediately. `M ≥ 2` constrains only a sparse arithmetic progression. This is the precise
mechanism behind THM-1530's observation that the `M ≥ 2` case is *structurally* different, and
it also says the conditions remain **infinite** in number against `d = M+N` free coefficients —
which is why TNC should still hold; a proof needs the sparse subsequence to be independent,
and that is what remains.

## 4. Verification

The criterion agrees with ground truth (`CT(Λ^m)` by exact integer arithmetic, `m ≤ 10`) on
every control: single weights `2u^{−1}`, `2u^{−3}` give `Π/t` constant; `u^{−1}+1`,
`u^{−1}+u`, `u^{−2}+1`, `u^{−2}+u^{−1}+1+u` all give `Π/t` non-constant and all have some
`CT(Λ^m) ≠ 0`. Search over `M ≤ 3`, `N ≤ 2`, coefficients in `{0,±1,±2,±i}`: **zero
counterexamples** (46,000+ polynomials), each candidate cross-checked with exact `CT`.

## 5. Retraction — a broken test, caught by its own control

My first run of this test reported `Λ = u^{−1}+u` as being in the nullcone (it is not:
`CT(Λ²) = 2`), and reported "hits" at `M = 1` where §2 **proves** there are none. Both were
false positives from one design error:

> `Π(t)/t → c` as `t → 0` for **every** `R` — that is just the leading asymptotic. Evaluating
> at `t = 10^{−3}…10^{−5}` and asking whether the values agree to `10^{−6}` measures
> *convergence*, not *constancy*. The reported "min spread" figures `1e−6, 1e−9, 1e−12` were
> simply `O(t^k)`.

Redone at moderate `t` (`0.05…0.01`), where a non-constant `Π/t` separates by `~10^{−2}`, the
controls all pass. **The criterion itself was never in question — it is an identity; only my
test of it was wrong.** Recorded because the failure mode is specific and reusable: *when a
quantity has a limit, testing "is it constant" near the limit point tests nothing.*

## 6. Honest status

GMC(2) is **not** finished. What the fleet now holds: one-signed/charge-definite (opus,
mac-mini), two-weight in every degree (THM-1510), `{−1,0,1}` (death-star, carrying my
Bessel-EMP), `W`-linear (boxeph), leading symbol with extreme weight `±1` (THM-1530). The
residue in every framing is the same object — **TNC with `M, N ≥ 2`** — and this file reduces
it to a single exact identity plus a stated independence question. That is a sharper handle,
not a closure.

*Files: `04-computation/tnc_exact_criterion_klein_S347.py` (+ `.out`).*
