---
id: THM-1580
title: "THE ARBORESCENCE COUNT STRICTLY DOMINATES H AS A FINGERPRINT, AND IT IS POLYNOMIAL-TIME. (A) Σ_r a_r and H are INCOMPARABLE for n ≥ 5 — neither determines the other — and Σ_r a_r is far finer: 298 distinct values against H's 77 over the 456 iso classes at n = 7, with the PAIR reaching 417. Since Σ_r a_r is a determinant (Tutte/Kirchhoff) while H is not, this is a strictly better fingerprint that is also strictly cheaper — the same shape as THM-506's (char, perm) result but with the computational advantage running the right way. (B) MY OWN S138 FRAMING IS REFUTED: the size-dependent shift in the ordinal-sum law is NOT about the prime 2. The law Σa(T₁⊕T₂) = Σa(T₁)·det(|T₁|·I + L_in(T₂)) holds for both parities of |T₁| and never mentions 2. The real cause is CROSSING MULTIPLICITY — a Hamiltonian path crosses the cut exactly once with no freedom, while an arborescence lets each of the |T₂| vertices independently choose a parent among the |T₁| across the cut. Evenness of Σa is a consequence, not a cause. (C) HYP-8315 EXTENDED TO n=8 AND ITS WORDING CORRECTED: the minimum is the transitive tournament at exactly (n−1)!, now confirmed n = 4..8; the maximiser is REGULAR at odd n and NEAR-REGULAR at even n — 'Paley maximises' does not even parse at even n, where no regular tournament exists. Maximiser score sequences: [1,1,2,2], [2,2,2,2,2], [2,2,2,3,3,3], [3,3,3,3,3,3,3], [3,3,3,3,4,4,4,4]. (D) The 2-adic law is Legendre and is not deep: v₂ is additive along the ordinal-sum tower, giving v₂(Σa(TT_n)) = (n−1) − s₂(n−1)."
status: >
  (A) VERIFIED-EXACT over all iso classes at n = 4,5,6,7 with exact integer determinants.
  Incomparability is exhibited at n = 5,6,7 (at n = 4 the class count is too small and
  Σa happens to refine H).
  (B) The ordinal-sum law is PROVED (THM-1460, block triangularity) and re-verified on 15
  factor pairs spanning both parities of |T₁|. The refutation of the prime-2 framing is a
  logical observation about that law, not a computation. The crossing-multiplicity account
  is an EXPLANATION, stated as such, not a separate theorem.
  (C) The minimum is EXHAUSTIVE at n ≤ 6 and hill-climbed (60 restarts) at n = 7,8 —
  so at n = 7,8 the MIN is confirmed against the transitive value but the MAX is a LOWER
  BOUND, not certified. Six data points now, not five.
  (D) PROVED (v₂ additive along the tower) and verified n = 3..11 against Legendre.
source: mac-mini-2026-07-20-S139 (owner: "work the arborescence leads")
depends_on:
  - THM-1460  # arborescences as the determinantal relaxation of H; the ordinal-sum law
related:
  - THM-506   # (char, perm) strictly dominates H -- the same shape, one tier up
  - THM-499   # H is finer than the adjacency spectrum
supersedes_in_part: HYP-8390   # part (b) refuted; part (c) answered
script: 04-computation/arborescence_leads_macmini_S139.py (+ .out)
---

# THM-1580 — the arborescence count is a better fingerprint than `H`, and cheaper

## (A) Incomparable, finer, and polynomial-time

| `n` | classes | distinct `H` | distinct `Σ_r a_r` | distinct pair | `Σa` refines `H`? | `H` refines `Σa`? |
|---|---|---|---|---|---|---|
| 4 | 4 | 3 | 4 | 4 | yes | no |
| 5 | 12 | 7 | 11 | 12 | **no** | **no** |
| 6 | 56 | 19 | 43 | 51 | **no** | **no** |
| 7 | 456 | **77** | **298** | **417** | **no** | **no** |

> For `n ≥ 5` the two invariants are **incomparable** — neither determines the other — and
> `Σ_r a_r` is far finer, separating **298** of the 456 classes at `n = 7` where `H` separates
> only 77. The pair reaches 417.

This is the shape of THM-506's `(char, perm)` result, **but with the computational advantage
running the right way**: `Σ_r a_r` is a determinant (Tutte/Kirchhoff, polynomial time) while
`H` is not. So the cheaper invariant is also the stronger one, and the pair is stronger still.
That is the practically useful sentence in this file.

## (B) The size-shift is crossing multiplicity, not the prime 2 — my own framing refuted

S138 (HYP-8390b) guessed that the prime-2 split explains why `log Σa` gains a size-dependent
shift while `log H` does not. **That is wrong, and the law itself says so:**

`Σa(T₁ ⊕ T₂) = Σa(T₁)·det(|T₁|·I + L_in(T₂))`

holds for **both parities** of `|T₁|` (verified on 15 factor pairs) and **never mentions 2**.
The actual cause is combinatorial:

- **`H`:** in `T₁ ⊕ T₂` every Hamiltonian path runs through all of `T₁` and then all of `T₂`,
  since nothing in `T₂` beats anything in `T₁`. It **crosses the cut exactly once, with no
  choice about where.** Hence `H(T₁⊕T₂) = H(T₁)H(T₂)` — no interaction term.
- **`Σa`:** an arborescence rooted in `T₁` lets **each of the `|T₂|` vertices independently
  choose a parent**, from the `|T₁|` vertices across the cut or from inside `T₂`. `|T₂|`
  independent crossings, each with `|T₁|` options — and that is exactly where `p = |T₁|`
  enters the shift.

So **evenness of `Σa` is a consequence (the option count `p` can be even), not a cause.**
HYP-8390b is refuted as stated.

## (C) The extremals at `n = 8`, and a correction to the wording

| `n` | transitive `= (n−1)!` | min found | min = transitive? | max found | maximiser scores |
|---|---|---|---|---|---|
| 4 | 6 | 6 | yes | 10 | `[1,1,2,2]` |
| 5 | 24 | 24 | yes | 55 | `[2,2,2,2,2]` |
| 6 | 120 | 120 | yes | 333 | `[2,2,2,3,3,3]` |
| 7 | 720 | 720 | yes | 2744 | `[3,3,3,3,3,3,3]` |
| 8 | 5040 | 5040 | yes | 23363 | `[3,3,3,3,4,4,4,4]` |

**The minimum is the transitive tournament at exactly `(n−1)!`, now on six data points.**

**The wording of HYP-8315 needs correcting.** It said "regular/Paley maximises," but at even
`n` **no regular tournament exists**, so that does not parse. The maximiser is **regular at
odd `n` and near-regular at even `n`** — the score sequence is as balanced as the parity
allows. (At `n = 7` the maximiser value 2744 is the Paley value, consistent with THM-1460's
closed form.)

## (D) The 2-adic law — clean, and not deep

`v₂` is additive along the ordinal-sum factorisation, and for the transitive tower
`TT_n = TT_{n−1} ⊕ •` the shift is exactly `p = n−1`, so

`v₂(Σa(TT_n)) = v₂((n−1)!) = (n−1) − s₂(n−1)` — Legendre. Verified `n = 3..11`.

There is no special 2-adic structure beyond Legendre acting on the size factors. HYP-8390c has
an answer and the answer is that the question was not deep.

## Honest scope

- **The `n = 7, 8` maxima are hill-climbed (60 restarts), hence LOWER BOUNDS, not certified.**
  The minima at those `n` are confirmed only against the transitive value — a smaller value
  could in principle exist outside the basins explored. `n ≤ 6` is exhaustive.
- (A)'s incomparability is a statement about **values**, not about which invariant is more
  *useful*: `H` carries the Rédei/OCF structure that `Σa` does not, and nothing here suggests
  replacing it. The claim is narrowly that as a *separator* `Σa` is finer and cheaper.
- (B)'s crossing-multiplicity account is an **explanation of a proved law**, not itself a new
  theorem. What is established is the negative: the prime-2 framing does not explain the shift.
- (C) still does not **prove** either extremal. HYP-8315 remains open; it now rests on six
  points rather than five, with corrected wording.
- No claim that `(H, Σa)` is a complete invariant — it separates 417 of 456 classes at `n = 7`,
  so **39 classes remain unseparated**.

*Artifacts:* `04-computation/arborescence_leads_macmini_S139.py` (+out).
