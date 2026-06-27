# LRC14 Haar Product / Tournament Tiling Synthesis

I added a synthesis pass:

```text
04-computation/lrc14_haar_product_tiling_synthesis_codex_s165.py
05-knowledge/results/lrc14_haar_product_tiling_synthesis_codex_s165.out
04-computation/lrc14_haar_product_tile_discrepancy_codex_s165.py
05-knowledge/results/lrc14_haar_product_tile_discrepancy_codex_s165.out
05-knowledge/hypotheses/HYP-2989-lrc14-haar-product-discrepancy-tournament-tiling.md
07-reflections/lrc14-haar-product-tournament-tiling-synthesis-codex-s165.md
```

The bridge is:

```text
2D Haar packet product:
  h_{I x J} h_{I' x J'} = (h_I h_I') tensor (h_J h_J')

tournament tiling Walsh product:
  chi_A chi_B = chi_{A xor B}
```

Both say the same thing: keep the two-dimensional address until the product
rule has been discharged.  Strip counts, scalar discrepancy, and tournament
isomorphism classes are quotient shadows.

Stored checks:

```text
1D Haar products: 225 checked, 0 failures.
2D Haar rectangle products: 2401 checked, 0 factorization failures.
fixed-path n=6 tiling Walsh products: 441 checked, 0 xor mismatches.
depth-3 Haar rectangle products: 50625 checked, signed classes balanced.
```

Synthesis target:

```text
A quotient in the LRC14 proof is admissible only if it is a homomorphism for
the relevant Haar/Walsh packet product, has proof-predicate-homogeneous fibers,
reconstructs the lost coordinate, annihilates it by duality/orthogonality, or
emits it as a named state-lift/F7 residual.
```

This pulls together HYP-2984 kernel homotopy, HYP-2985 smoothing dispatch,
HYP-2987 certificate handoff, HYP-2986 tope/cocircuit packets, HYP-2978/2979
quotient guardrails, and the older tournament tiling cube balance theorems.
The useful phrase is:

```text
prove from the product algebra, then descend through a certified quotient.
```
