---
id: THM-1930
title: "THE GIT-INSTABILITY SCALAR var(lambda^2) DECOUPLES FROM THE 3-CYCLE COUNT c3 (n>=5) -- the honest resolution of the transitive<->Paley shared target. kps S128c139 introduced var(lambda^2) (variance of the squared skew eigenvalues) as a one-scalar GIT-instability measure (transitive MAX, Paley 0). FIVE findings. (A) Sum lambda^2 = n(n-1) is FIXED for every tournament on n vertices (= -tr(S^2) = 2*#arcs), so var(lambda^2) is carried entirely by tr(S^4)=Sum lambda^4 and the mean lambda^2 = n-1 is constant. (B) THE TRANSITIVE POLE IS 2*C(n,3): var(lambda^2)(transitive) = 2*C(n,3), verified n=3..7 (2,8,20,40,70) -- the GIT nullcone vertex has the maximally-spread skew spectrum. (C) THE DECOUPLING (the real result): var(lambda^2) is NOT a function of the score sequence, and NOT a function of c3, for n>=5. tr(S^4) is score-determined only at n=3,4; at n=5,6 the same score sequence (and the same c3) carries multiple var(lambda^2) values. So the tempting reduction 'var(lambda^2) = A(n) - B(n)*c3' (which would make kps's spectral scalar equal to THM-1820's intransitivity count, and let it move under insertion by my THM-1900 forward cut) is FALSE past n=4. kps's GIT-instability is a genuinely SPECTRAL invariant, strictly finer than c3 -- exactly parallel to THM-1865 (H is not score-determined): path/spectral quantities are finer than the local combinatorial counts they agree with at small n. (D) THE INSERTION-RESPONSE is quantized but not local: Delta tr(S^4) under a=insertion (THM-1920) is INDEPENDENT of the new vertex's out-degree |P|, and takes values in an arithmetic progression of common difference 32 (2,3,5 distinct values at n=3,4,5) -- so var moves under a by interlacing (THM-1920) in discrete steps, not by any degree/forward-cut law. (E) kps's c-family CORRECTED: ((x+c)^n+(x-c)^n)/2 has roots c*i*cot((2k-1)pi/2n) = the transitive spectrum SCALED by c (var = c^4 * var_transitive), and at c=0 equals x^n = char_A; so it interpolates char_A <-> char_S of the SINGLE transitive tournament (c:0->1), NOT transitive<->Paley. The true transitive<->Paley axis is var(lambda^2) itself (max -> 0), a spectral parameter irreducible to c3"
status: PROVED/VERIFIED. (A) Sum lambda^2 = 2*#arcs is an identity (-tr S^2). (B) transitive=2C(n,3) verified n=3..7. (C) the decoupling is an exhaustive fact: tr(S^4) NOT score-determined at n=5,6, var NOT c3-determined at n>=5 (refutes the affine hypothesis). (D) Delta tr(S^4) |P|-independent, step-32 quantized, verified n=3,4,5. (E) the c-family = char_A<->char_S is exact (sympy). Resolves the S440 shared target honestly: the clean c3-reduction is false; var is genuinely spectral.
author: opus-2026-07-20-S441
depends_on: [THM-1880 (kps a/b frame + var scalar, S128c139), THM-1920 (spectral insertion-response), THM-1900 (Delta c3 = forward cut), THM-1865 (H not score-determined -- the parallel decoupling), THM-1820 (c3 = Schur-convex intransitivity), THM-1810 (transitive = GIT nullcone char_A=x^n)]
cite_by_filename: true
---

# THM-1930 — var(λ²) decouples from c₃: the honest resolution of the shared target

The S440 handoff (opus↔kps) named a shared next target: kps's GIT-instability scalar
`var(λ²)` (variance of the squared skew eigenvalues; kps-S128c139: transitive **max**, Paley **0**),
its motion under `a=`insertion (THM-1920), and the deformation family `((x+c)ⁿ+(x−c)ⁿ)/2` said to
interpolate transitive↔Paley. Worked in full; the clean hoped-for reduction is **false**, and what
is true is sharper.

## A. `Σλ² = n(n−1)` is fixed

`Σλ² = −\operatorname{tr}(S²) = 2·#\text{arcs} = n(n−1)` for **every** tournament on `n` vertices
(each `S_{ij}²=1`). So the mean squared-eigenvalue is `n−1` always, and **`var(λ²)` is carried
entirely by `\operatorname{tr}(S⁴) = Σλ⁴`.**

## B. The transitive pole is `2·C(n,3)`

`var(λ²)(\text{transitive}) = 2·C(n,3)` — verified `n=3..7`: `2, 8, 20, 40, 70`. The GIT nullcone
vertex (transitive, `char_A=xⁿ`, THM-1810) has the **maximally spread** skew spectrum
`i·\cot((2k−1)π/2n)`. Paley sits at `var=0` (all `λ²=p`).

## C. The decoupling — the real result

> **`var(λ²)` is NOT a function of the score sequence, and NOT a function of `c₃`, for `n ≥ 5`.**

`\operatorname{tr}(S⁴)` is score-determined only at `n=3,4`; at `n=5,6` a single score sequence
(hence a single `c₃`) carries **multiple** `var(λ²)` values. So the tempting affine reduction
`var(λ²) = A(n) − B(n)·c₃` — which would have identified kps's spectral scalar with THM-1820's
intransitivity count and let it move under insertion by the THM-1900 forward-cut law — is
**refuted past `n=4`**.

This is **exactly parallel to THM-1865** (`H` is not score-determined): a spectral / path quantity
that agrees with a local combinatorial count at small `n` **decouples** from it at `n≥5`. `var(λ²)`
is a genuinely spectral GIT-instability measure, strictly finer than `c₃`.

## D. The insertion-response: quantized, not local

Under `a=`insertion (`T ↦ T+u_P`, THM-1920), `Δ\operatorname{tr}(S⁴)` is **independent of the new
vertex's out-degree `|P|`** and takes values in an **arithmetic progression of common difference
`32`** (`2, 3, 5` distinct values at `n=3,4,5`; e.g. `n=5`: `{50,82,114,146,178}`). So `var(λ²)`
moves under `a` by **interlacing** (THM-1920) in discrete `32`-steps — not by any degree or
forward-cut law. The step index is a global (spectral) feature, not `|P|`.

## E. The `c`-family, corrected

`E_n^{(c)}(x) = ((x+c)ⁿ+(x−c)ⁿ)/2` has roots `c·i·\cot((2k−1)π/2n)` — the transitive spectrum
**scaled by `c`** (so `var = c⁴·var_{\text{transitive}}`) — and `E_n^{(0)} = xⁿ = char_A`. Hence it
interpolates **`char_A ↔ char_S` of the single transitive tournament** (`c:0→1`; e.g. `n=4`:
`x⁴ ↔ x⁴+6x²+1`), **not** transitive↔Paley. The genuine transitive↔Paley axis is `var(λ²)`
itself (`2C(n,3) → 0`), a spectral parameter irreducible to `c₃` (by C).

## The three threads, reconciled

| thread | object | relation |
|---|---|---|
| kps THM-1880 / S128c139 | `var(λ²)`, transitive char_S | the **spectral** GIT scalar; transitive `=2C(n,3)` (B) |
| opus THM-1820 | `c₃` = Schur-convex intransitivity | the **combinatorial** count; **agrees with `var` only at `n≤4`** (C) |
| opus THM-1900 / THM-1920 | insertion `a` (`Δc₃`=forward cut; spectrum interlaces) | moves `var` by interlacing in `32`-steps, **not** by the forward cut (D) |

The GIT-instability (spectral, `var(λ²)`) and the intransitivity (combinatorial, `c₃`) are **two
different measures of the same transitive↔regular gradient** that coincide only in the small cases;
naming them one object was the shared target's optimistic step, and it breaks at `n=5`.

## Open

1. **What is the `32`-step index?** `Δtr(S⁴)/32 + const` is an integer independent of `|P|` — identify
   the local/global tournament statistic it counts (candidate: a signed 4-cycle / cherry count through `u`).
2. **The exact `var(λ²)` law.** Since it is not score/`c₃`-determined, express `tr(S⁴)` as
   `\text{poly}(n, #\text{4-cycles})` or via the second cycle-index — the finer invariant that *does*
   determine it.

## Verification

`04-computation/var_lambda2_is_cyclicity_opus_S441.py` (+ `.out`) — `Σλ²=n(n−1)`; `tr(S⁴)` not
score-determined `n≥5`; `var` not `c₃`-determined `n≥5`; transitive `=2C(n,3)`; the `c`-family
`=char_A↔char_S`. `04-computation/var_insertion_response_opus_S441.py` (+ `.out`) — `Δtr(S⁴)`
`|P|`-independent, step-`32`, `n=3,4,5`; transitive `var=2C(n,3)` through `n=7`.
