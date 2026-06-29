---
id: HYP-3544
title: The antipodal Euler number SC(n)=P_n(-1) refines into Z_2-equivariant Betti numbers of the level-graded S_n-invariant chain complex of the arc-hypercube; the R-odd Betti numbers are the Borsuk-Ulam obstruction — computable by Kaczynski-Mischaikow-Mrozek computational homology (engineering deliverable)
status: OPEN (klein-2026-06-29-S2) — Euler characteristic (=SC) known in closed form (THM-587); the boundary operator and Betti numbers not yet computed
source: klein-2026-06-29-S2
depends_on:
  - THM-587   # P_n(-1)=SC(n) = the antipodal Euler characteristic
  - THM-584
related:
  - THM-224   # simplicial up-Laplacian (the boundary/Laplacian machinery)
  - HYP-3538
  - HYP-2983  # the repo's prior "Kaczynski" (analytic-sieve / iterated-projection) thread
  - HYP-3543
---

# HYP-3544 — Equivariant homology refines the antipodal Euler number

## Claim

`SC(n) = P_n(-1) = sum_k (-1)^k mult(k)` (THM-587) is an Euler characteristic: the alternating sum of the
dimensions of the level-graded `S_n`-invariant spaces `C_k = W_k^{S_n}` (`dim C_k = mult(k)`, `W_k =`
level-`k` Boolean-Fourier subspace of `Q_d`). Therefore there is a chain complex `(C_*, partial)` with

```
chi(C_*) = sum_k (-1)^k dim C_k = SC(n),
```

and `SC(n) = sum_k (-1)^k b_k` where `b_k` are its Betti numbers. The natural boundary is the cube's
signed "remove one arc" down-operator restricted to `S_n`-invariants (the simplicial boundary of the
`d`-simplex / the THM-224 up-Laplacian's adjoint), twisted by the antipodal `Z_2`. Conjecture: the
`Z_2`-equivariant homology splits `H_* = H_*^{R-even} (+) H_*^{R-odd}`, and the **`R`-odd Betti numbers are
the Borsuk-Ulam obstruction** — the homological content the Euler number only sees with sign.

## Why it matters

- The Euler number `SC(n)` is the *shadow*; the Betti numbers are the *object*. Two complexes with the
  same `chi=SC` differ in their homology, and the `R`-odd homology is where the project's signed
  obstruction (HYP-3538 `M_odd`, the witness `phi`) becomes a topological invariant rather than a number.
- **Engineering deliverable.** This is a concrete `Z_2`-equivariant homology computation on a level-graded
  complex — exactly the Kaczynski-Mischaikow-Mrozek computational-homology toolkit, and the project's
  `circulant_homology` / `tournament_tda` mandate. The complex is small per level (`dim = mult(k)`, all
  computable from `n!` permutations via THM-587), so n=7,8,9 are reachable.

## Next steps

1. Build the boundary `partial_k: C_k -> C_{k-1}` (signed down-operator on invariants) and verify
   `partial^2 = 0`; confirm `sum (-1)^k dim C_k = SC(n)` (already true by THM-587).
2. Compute `b_k` for n=4..8; split by `R`-parity (even/odd levels); identify the `R`-odd homology.
3. Relate to THM-224 (up-Laplacian spectrum) — the metagraph adjacency `A = U+D` is `up+down`; the
   homology is `ker/im` of `D` on invariants.
4. Connect the `R`-odd Betti numbers to Rédei's odd index / Ky Fan's alternating count (the
   `the-per-level-signed-cycle-index-borsuk-ulam-ky-fan` reflection) — is `b_k^{odd}` the obstruction the
   LRC/witness threads have been seeking?

## Note on the two "Kaczynski" readings

The repo's prior HYP-2983 uses "Kaczynski" for an analytic-sieve / exponential-sum (and iterated
projection) template. That meets this thread at the **Reynolds projection** `(I+R)/2`: projecting onto the
`R`-even merged metagraph is the Kaczmarz-style iterative projection onto the invariant subspace, the
linear-algebra twin of the 2-adic descent's iterated Reynolds step. The computational-homology reading
(Kaczynski-Mischaikow-Mrozek) is the topological refinement of the same `R`-even/`R`-odd split.
