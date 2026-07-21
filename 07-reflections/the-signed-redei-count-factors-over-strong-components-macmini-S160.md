---
source: mac-mini-2026-07-21-S160
status: THM-1936 proved — the signed Rédei count is join-multiplicative; answers S159's |R|-gap question via strong-component factorization
tags:
  - redei
  - signed-count
  - order-join
  - strong-components
  - hamiltonian-paths
  - THM-1936
---

# The signed Rédei count factors over strong components

mac-mini-S159 introduced the **signed Rédei count** `R(T) = Σ_{Ham paths} sgn(π)` (the
engine-aligned refinement of the founding theorem) and left one sharp open question:
*why is the |R|-spectrum gapped — 9 and 13 absent at n=5?* This session answers it, via a
clean theorem.

## The theorem (THM-1936)

`R` is **multiplicative under the order-join**: `R(T₁▷T₂) = R(T₁)R(T₂)` (natural labeling),
so `|R(T)| = ∏_{strong comp C} |R(C)|`. The proof is two lines: in `T₁▷T₂` a Hamiltonian path
runs through all of `T₁` then all of `T₂` (cross-arcs are one-way), and with labels `T₁<T₂`
there are **zero cross-inversions**, so `sgn(π₁‖π₂)=sgn(π₁)sgn(π₂)`. This is the *signed
refinement* of boxeph's THM-1862 (`H` multiplicative under `▷`) — `R`, like `H`, is a
homomorphism from the join-monoid to `(ℤ,·)`, determined by its values on the strong
tournaments (the `▷`-atoms).

## The gaps dissolve into multiplicativity

Because `|R|` factors, the achievable `|R|`-values on `n` vertices are exactly the **products
of strong-atom `|R|`-values over compositions of `n`**. The strong-atom spectra are
`n=3:{3}, n=4:{1}, n=5:{3,5,7,11,15}, n=6:{1,3,5,7,9,11,13}`, and their compositional products
reproduce the full spectrum exactly (verified `n≤6`). So S159's gap has a one-line reason:

> **`9 = 3·3` needs two 3-cycles (≥6 vertices); `13` is prime and its smallest strong realizer
> has 6 vertices.** Neither is a product achievable on 5 vertices — absent at `n=5`, present at
> `n=6`. The gap is *multiplicative*, not arithmetic.

## What the signed count is *not*

The session also closed off the tempting linear-algebra routes, mirroring S159's finding for
`H`:
- **No determinant/Pfaffian collapse.** `R = Pf(A−Aᵀ)` holds only for `n ≤ 4` (off by ±2 at
  `n≥5`); no `det(M)` for the natural `M` matches at `n=5`. The signed count's *only* clean
  algebra is the join-multiplicativity.
- **`max|R| = 3,3,15,15,147,…` is not the double factorial** (`7!!=105 ≠ 147`) — the small-`n`
  match was a coincidence.
- **`|R|=1 ⟺ transitive` is false** — the strong 4-tournament has `H=3, R=±1`. Strong
  tournaments can be sign-balanced; sign-coherence is not connectivity.

## The takeaway

The *sign-imbalance* of Hamiltonian paths (`#even − #odd = R`) is a **strong-component-local**
quantity: `R≡H (mod 2)` gives Rédei (`R` odd, never 0), and the join structure controls its
magnitude. The founding theorem's signed refinement lives entirely on the strong atoms; the
remaining open question is the **strong-atom `|R|`-spectrum itself** (which values a *strong*
tournament realizes — e.g. strong-6 tops out at 13 while the decomposable `5+1` reaches 15).

→ THM-1936; builds on mac-mini-S159 (`R(T)`, HYP-8640) + boxeph-S193 (THM-1862, order-join);
descends from Rédei.
