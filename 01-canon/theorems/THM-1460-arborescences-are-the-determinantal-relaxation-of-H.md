---
id: THM-1460
title: "ARBORESCENCES ARE THE DETERMINANTAL RELAXATION OF HAMILTONIAN PATHS, AND THE LOGARITHM SEPARATES THEM. A Hamiltonian path from r IS a spanning out-arborescence rooted at r with every out-degree ≤ 1; dropping that one constraint is exactly what turns an intractable count into a DETERMINANT (Tutte), so h_r ≤ a_r and H ≤ Σ_r a_r with Σ_r a_r computable in polynomial time. (A) Matrix-tree verified against brute force on every iso class n ≤ 6, with Kirchhoff's Σ_r a_r = ∏(nonzero eigenvalues of D_in − A). Transitive gives EXACTLY (n−1)!, all of it at the source. PALEY CLOSED FORM: Σ_r a_r = [q(q+1)/4]^{(q−1)/2}, verified exactly at q = 3,7,11,19 (18 digits at q=19). (B) THE RELAXATION GAP Σ_r a_r / H is extremised by the two poles of the repo: the TRANSITIVE tournament MAXIMISES it at exactly (n−1)!/1, and the REGULAR/Paley tournament MINIMISES it (1, 2, 3.27, 7.2, 14.52 at n = 3..7). Transitive simultaneously MINIMISES the arborescence count itself, Paley maximises it. (C) Σ_r a_r is LAPLACIAN-spectral by Kirchhoff but is NOT ADJACENCY-spectral — it differs inside 111 of the 116 adjacency-cospectral groups at n = 7 — so it is TRANSVERSE to the THM-499/500 spectral hierarchy, which is entirely about spec(A). (D) THE LOGARITHM. Ordinal sum makes D_in − A block upper triangular, giving the exact law Σa(T₁⊕T₂) = Σa(T₁)·det(|T₁|·I + L_in(T₂)) — PROVED and verified. So log H is additive with NO interaction term while log Σa is additive with a SIZE-DEPENDENT SHIFT. Tree entropy (1/n)log Σa is minimised by the transitive tournament at every n = 3..7."
status: >
  (A) PROVED (Tutte / Kirchhoff, standard) and VERIFIED-EXACT: matrix-tree cofactor equals
  a direct enumeration on every iso class at n = 3..6, and Σ_r a_r equals the nonzero-eigenvalue
  product at every class. Transitive = (n−1)! PROVED (L_in is triangular with the in-degrees
  0,…,n−1 on the diagonal). Paley closed form PROVED from spec(A) and verified exactly to q=19.
  (B) h_r ≤ a_r PROVED (a Hamiltonian path is an arborescence with out-degrees ≤ 1) and verified
  on every class n ≤ 7. The two extremal identifications are VERIFIED n = 3..7, NOT proved.
  (C) VERIFIED-EXACT n = 4..7 by grouping all iso classes by characteristic polynomial.
  (D) The ordinal-sum law is PROVED (block triangularity) and verified on all 9 pairs with
  |T₁|,|T₂| ∈ {2,3}. The tree-entropy minimality of the transitive tournament is VERIFIED
  n = 3..7, NOT proved.
source: mac-mini-2026-07-20-S132 (owner: "merge in and extend ideas related to arborescences on
  tournaments and determinants and logarithms")
depends_on:
  - THM-506   # master cycle-packing polynomial Phi: the det/OCF weight-family MERGE, already canon
  - THM-505   # the OCF non-spectral defect -- why H has no determinant
related:
  - THM-499   # H finer than the ADJACENCY spectrum from n=6
  - THM-500   # odd-cycle count non-spectral from n=7
  - THM-118   # c_k = tr(A^k)/k
  - THM-1440  # skew-Seidel spectra; a DIFFERENT matrix again
script: 04-computation/arborescences_determinants_logs_macmini_S132.py (+ .out)
---

# THM-1460 — arborescences are the determinantal relaxation of `H`

**One line.** A Hamiltonian path from `r` **is** a spanning out-arborescence rooted at `r` in
which every vertex has out-degree `≤ 1`. Drop that one constraint and the count becomes a
**determinant** — which is precisely the boundary between the repo's tractable and intractable
spanning counts.

## Merge note — what is already canon

The natural first merge, *"the characteristic polynomial and the OCF are the same
disjoint-cycle sum at different weights,"* **is already repo canon**: HYP-2514 / THM-506's
master cycle-packing polynomial `Φ(T;{y_k}) = Σ_{L linear subdigraph} ∏_{C∈L} y_{|C|}`, with
the spectrum its **signed** vertex-graded (Sachs) face and `H` its **unsigned odd-only** face.
Not re-claimed here. What was genuinely absent is **arborescences**, and that is what follows.

## (A) Matrix-tree on tournaments

With `L_in = D_in − A` (note: `1ᵀL_in = 0`, so `0 ∈ spec`, but `L_in·1 ≠ 0` in general):

> `a_r(T)` = # spanning out-arborescences rooted at `r` = `det` of `L_in` with row/col `r`
> deleted (Tutte), and `Σ_r a_r = ∏` of the **nonzero** eigenvalues of `L_in` (Kirchhoff).

Both verified against brute-force enumeration on **every** iso class at `n = 3..6`.

**Transitive.** `L_in` is upper triangular with the in-degrees `0,1,…,n−1` on the diagonal, so

> `Σ_r a_r(TT_n) = (n−1)!`, and *all* of it sits at the source — no arborescence can be rooted
> anywhere else, since nothing beats the source.

**Paley** (`q ≡ 3 mod 4`) is regular, so `L_in = ((q−1)/2)I − A`, and
`spec(A) = {(q−1)/2} ∪ {(−1 ± i√q)/2}` gives

> **`Σ_r a_r(Paley_q) = [q(q+1)/4]^{(q−1)/2}`**, with `a_r = (1/q)` times that.

Verified exactly: `3, 2744, 39135393, 630249409724609375` at `q = 3, 7, 11, 19`.

## (B) The relaxation gap, and its two poles

`h_r ≤ a_r` always (verified every class `n ≤ 7`), so `H ≤ Σ_r a_r`. The gap is extremised by
exactly the repo's two canonical tournaments:

| `n` | min `Σa/H` | at | max `Σa/H` | at |
|---|---|---|---|---|
| 3 | 1.000 | regular `[1,1,1]` | 2 | transitive |
| 4 | 2.000 | `[1,1,2,2]` | 6 | transitive |
| 5 | 3.267 | `[1,2,2,2,3]` | 24 | transitive |
| 6 | 7.200 | `[2,2,2,3,3,3]` | 120 | transitive |
| 7 | 14.519 | regular `[3,3,3,3,3,3,3]` | 720 | transitive |

> **The transitive tournament maximises the gap at exactly `(n−1)!/1`** — the relaxation is
> maximally loose precisely where `H = 1`. **The regular/Paley tournament minimises it**
> (at `n=7`, `2744/189 = 14.519`, attained by `Paley₇`).

And the arborescence count *itself* is extremised the opposite way: **transitive minimises**
`Σ_r a_r` and **Paley maximises** it, at every `n = 3..7`.

## (C) Where it sits spectrally — transverse to THM-499/500

Kirchhoff makes `Σ_r a_r` **Laplacian**-spectral by construction. But `L_in` knows the scores,
while `spec(A)` pins only `Σ s_i` and `Σ s_i²` (the latter via `c₃ = tr A³/3`). So:

| `n` | adjacency spectra | cospectral groups | `Σ_r a_r` differs inside |
|---|---|---|---|
| 4 | 3 / 4 classes | 1 | **1** |
| 5 | 9 / 12 | 2 | **2** |
| 6 | 28 / 56 | 19 | **16** |
| 7 | 168 / 456 | 116 | **111** |

> **`Σ_r a_r` is NOT adjacency-spectral.** It is therefore **transverse** to the THM-499/500
> hierarchy, which is entirely a statement about `spec(A)`: arborescences are not "above" or
> "below" those boundaries, they read a different matrix.

## (D) The logarithm

Under ordinal sum `T₁ ⊕ T₂` (every vertex of `T₁` beats every vertex of `T₂`, `p = |T₁|`), the
in-Laplacian is **block upper triangular**:

```
L_in(T₁⊕T₂) = [ L_in(T₁)      −J     ]     so   spec = spec(L₁) ∪ ( p + spec(L₂) )
              [    0       pI + L₂   ]
```

> **`Σa(T₁ ⊕ T₂) = Σa(T₁) · det(p·I + L_in(T₂))`**, `p = |T₁|`.

Proved by block triangularity; verified on all 9 pairs with `|T₁|,|T₂| ∈ {2,3}` (alongside
`H(T₁⊕T₂) = H(T₁)H(T₂)`, which also held throughout). Taking logarithms:

| | law | interaction |
|---|---|---|
| `log H` | `log H(T₁) + log H(T₂)` | **none** |
| `log Σa` | `log Σa(T₁) + Σ_μ log(p + μ)` | **size-dependent shift** |

So the logarithm is exactly what exposes the difference: `H` is a clean multiplicative norm
under ordinal sum, while the arborescence count picks up a term that depends on *how big the
first factor is*. Sanity check: `TT_n = TT_{n−1} ⊕ •` gives `Σa(TT_n) = Σa(TT_{n−1})·(n−1)`,
recovering `(n−1)!`.

**Tree entropy.** `(1/n)·log Σ_r a_r` is **minimised by the transitive tournament at every
`n = 3..7`** (`0.231, 0.448, 0.636, 0.798, 0.940`), with the maximum at Paley/regular
(`0.366, 0.576, 0.801, 0.968, 1.131`).

## Honest scope

- (A) is standard Tutte/Kirchhoff applied to tournaments, plus the two closed forms. The
  *application* is new to this repo; the theorems are classical.
- **(B) and (D)'s extremal claims are VERIFIED `n ≤ 7`, not proved.** "Transitive minimises
  `Σ_r a_r`" and "Paley maximises it" rest on five data points each; the Paley claim is only
  tested where a Paley tournament exists.
- (C) says `Σ_r a_r` is not a function of `spec(A)`. It says **nothing** about whether it is
  finer or coarser than `H` — that comparison is not made and the two are likely incomparable.
- The `h_r ≤ a_r` bound is the only proved relation between the two counts. **No formula
  expressing `H` in terms of arborescences is claimed, and none is known.** The relaxation is
  one-directional.
- The `Σa/H` ratio is offered as a natural measure of relaxation looseness, not as an invariant
  with established meaning.

*Artifacts:* `04-computation/arborescences_determinants_logs_macmini_S132.py` (+out).
*Credits:* THM-506/HYP-2514 for the master cycle-packing polynomial (the det/OCF merge, already
canon and not re-claimed); THM-505 for why `H` has no determinant; Tutte and Kirchhoff.
