---
id: THM-1625
title: "TNC SPLIT INTO THREE, TWO OF THEM CLOSED. (1) VANDERMONDE LEMMA, PROVED: if the DOMINANT saddle values w_j are DISTINCT then TNC holds — vanishing for all m forces sum_j c_j w_j^m = 0 for m = 1..k, a system with det = prod_{i<j}(w_j-w_i) * prod_j w_j != 0, hence every c_j = 0; but c_j = 1/(u_j sqrt(2 pi g''(u_j))) is NEVER zero at a nondegenerate saddle. So the ONLY escape is COLLIDING dominant values. (2) COLLISIONS ARE REAL BUT NOT FATAL: a 500-quartic sweep at N=2 found 30 (6%) with colliding dominant values, and in every one the prefactor sum is EXACTLY 0 — yet CT(Lambda^m) is still nonzero (R = u^4-2u^2-2 gives -2, 0, 16, -56, 48, 384, -1920). Leading-order cancellation is NECESSARY, NOT SUFFICIENT; the subleading term takes over. (3) SYMMETRIC COLLISIONS DESCEND, verified exactly: if R(u) = S(u^k) then [u^{Nm}]R^m = 0 unless k | Nm and otherwise EQUALS [v^{Nm/k}]S^m, a smaller instance — at N=2, k=2 that is exactly N'=1, which THM-1530(B) PROVED by Lagrange-Burmann (matched m=1..6 on the nose). (4) HONEST RESIDUAL: NOT every collision is symmetric — the sweep's k-fold symmetries are {1,2,4}, so k=1 ASYMMETRIC collisions exist. Those, and only those, are what remains of the TNC"
status: PROVED (1); VERIFIED (2),(3); (4) is the honest open residual
author: opus-2026-07-20-S417
depends_on: [THM-1615 (a genuine saddle always exists), THM-1530 (klein: N'=1 by Lagrange-Burmann), THM-1595 (boxeph: ladder closures), THM-1550]
---

# THM-1625 — TNC split into three, two of them closed

## Setup

From THM-1615, `CT(Λ^m) = [u^{Nm}]R^m = (1/2πi)∮e^{m·g(u)}du/u` with
`g = log R − N log u`, and a **genuine saddle always exists**. The standard expansion gives

```
CT(Λ^m)  ~  Σ_j  c_j · w_j^m / √m ,     w_j = R(u_j)/u_j^N ,
c_j = 1 / ( u_j · √(2π g″(u_j)) )
```

summed over the **dominant** saddles (those of maximal `|w_j| = ρ`).

## 1. The Vandermonde lemma (PROVED)

> **Lemma.** If the dominant saddle values `w_j` are **distinct**, then TNC holds for `Λ`.

*Proof.* Vanishing for all `m` forces `Σ_j c_j w_j^m = 0` for `m = 1,…,k`. That system has
matrix `M_{mj} = w_j^m = w_j·w_j^{m−1}`, i.e. `M = V·diag(w_j)` with `V` Vandermonde, so

```
det M  =  ∏_{i<j}(w_j − w_i) · ∏_j w_j  ≠  0
```

for distinct nonzero `w_j`. Hence every `c_j = 0`. But at a nondegenerate saddle
`u_j ≠ 0` and `g″(u_j)` is finite and nonzero, so `c_j ≠ 0`. Contradiction. ∎

**So the only possible escape is colliding dominant values.**

## 2. Collisions are real — and prefactors do cancel — but that is not enough

A sweep of 500 quartics `R` at `N = 2`:

| | count |
|---|---|
| with saddle data | 500 |
| **with colliding dominant values** | **30 (6.0 %)** |

and in **every** collision case the prefactor sum is **exactly `0`** — the leading order
cancels. Examples: `u⁴−2u²−2`, `u⁴−u²−2`, `u⁴−2`, `u⁴+u²−2`.

**But the constant terms are not zero.** For `R = u⁴−2u²−2`:

```
CT(Λ^m),  m = 1..7  =  −2, 0, 16, −56, 48, 384, −1920
```

> **Leading-order cancellation is necessary but not sufficient.** When the dominant term
> dies, the subleading term takes over and `CT` survives. Any attempt to close TNC by
> exhibiting prefactor cancellation alone is therefore incomplete — this is the sharpest
> practical warning in this note.

## 3. Symmetric collisions descend (verified exactly)

If `R(u) = S(u^k)` then `R^m` is a polynomial in `u^k`, so

```
[u^{Nm}]R^m  =  0                     unless k | Nm,
             =  [v^{Nm/k}] S(v)^m     when k | Nm.
```

The right side is a **strictly smaller TNC instance**. At `N = 2, k = 2` it is exactly
`N′ = 1` — the min-exponent `−1` case **proved by THM-1530(B)** via Lagrange–Bürmann.

**Verified on `R = u⁴−2u²−2`, `S(v) = v²−2v−2`:**

| `m` | `[u^{2m}]R^m` | `[v^m]S^m` |
|---|---|---|
| 1 | −2 | −2 |
| 2 | 0 | 0 |
| 3 | 16 | 16 |
| 4 | −56 | −56 |
| 5 | 48 | 48 |
| 6 | 384 | 384 |

Exact match. This also explains the `m = 2` zero: `k ∤ Nm` fails there. **Symmetry kills a
congruence class of `m`, never all of them** — consistent with THM-1530(C).

## 4. The honest residual: asymmetric collisions

I hoped every collision came from a symmetry. **It does not.** Computing the largest `k` with
`R(u) = S(u^k)` across the 30 collision cases gives

```
k-fold symmetries observed:  {1, 2, 4}
```

`k = 1` means **no symmetry at all**, yet the dominant values still collide. So the clean
two-case proof (distinct ⟹ Vandermonde; colliding ⟹ symmetric ⟹ descent) **does not close**.

> **What remains of the TNC is exactly: `Λ` whose dominant saddle values collide for a reason
> other than a `u ↦ ζu` symmetry.** In the sweep these still had `CT ≠ 0` — the subleading
> term saved them — but that is an observation, not a proof.

## 5. Status

| case | status |
|---|---|
| dominant values distinct | **PROVED** (§1, Vandermonde) |
| collision from symmetry `R = S(u^k)` | **descends** to a smaller instance (§3); at `N=2,k=2` lands on THM-1530(B) |
| collision without symmetry (`k = 1`) | **OPEN** — the whole residual |

This is a genuine narrowing of HYP-8450 (which asked the vanishing-sums question in full): the
distinct-value case is now closed outright, and the symmetric case reduces. It is **not** a
closure of TNC, and should not be cited as one.

## 6. Next

1. **Characterise asymmetric dominant-value collisions.** `w(u) = R(u)/u^N` taking equal
   critical values at critical points not related by a symmetry — a question about the
   critical values of a rational function, i.e. about its **branch data**. The natural tool
   is the monodromy/Hurwitz description of `w`, not more elimination.
2. **Second-order Vandermonde.** When the leading term cancels, the subleading coefficients
   satisfy their own linear system. If *that* matrix is also nonsingular for distinct `w_j`,
   the argument iterates and §4 closes too. This looks the most promising and is concrete.

## Verification

`04-computation/tnc_vandermonde_opus_S417.py` (the lemma; the collision sweep; the
prefactor-sum and `CT` values), `04-computation/tnc_descent_opus_S417.py` (the `k`-fold
census showing `k ∈ {1,2,4}`; the exact descent match `m = 1..6`). Outputs in
`05-knowledge/results/`.
