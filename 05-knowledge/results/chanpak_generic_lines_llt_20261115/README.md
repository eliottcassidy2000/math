# Chan–Pak: generic collision lines, a local limit law for αH + βP₂, and host operators

Session klein-2026-11-15-llt. **Chan–Pak remains OPEN. Neither the generic-line equidistribution nor a local limit
theorem is proved here.** What is proved, measured and reduced:

1. **[Generic single-point lines](GENERIC_LINES.md)** (the open core of HYP-9752 after THM-5232).
   - *PROVED.* For a fixed signed difference σ between the two middle blocks of a collision, p(x) − p(y) is
     **affine in the shared bits**, with explicit coefficients: partial forward and reversed differences. So the
     exact collision count is a sum of **Riesz products** Π(1 + e(θc_k)) against two Dirichlet kernels.
   - *NUMERICAL.* Generic cores hit their lattice lines at 1.00–1.15 times the random-divisibility rate (spans
     10–13, suffixes 0–4); divisibility itself is unbiased.
   - *Why it resists.* Per signature the expected count is about 2^(−a/2), so equidistribution can only come from
     cancellation **across** signatures.
   - *REDUCED.* The missing estimate is stated precisely as a large-sieve bound for these digital Riesz products,
     [HYP-9782](../../hypotheses/HYP-9782-generic-line-large-sieve-the-open-core-of-the-row-collision-lemma-is-a-large-sieve-estimate-for-digital-riesz-products-across-signatures.md).
     It is plausible because the Mahler measure of 1 + z is 1.
2. **[A local limit law for αH + βP₂](LLT.md).**
   - *Exact data (orders 8–10).* Per-line class counts of 3H + P₂ fluctuate by √(D/μ), with dispersion D ≈ 2–4.
     The relative fluctuation falls from 1.54 to 0.23 to 0.17, with the worst bulk line at 0.25–0.37 of the smooth
     density.
   - *Flat singular series.* All residue classes modulo 8, 16, 32, 3, 9, 27, 5, 7 and 24 are within 1%.
   - *No local proof.* S lies on 6ℤ + 3 at odd order, and a unit change of S is a vanishingly rare local move. So
     **no switching proof through single moves** is possible.
   - *Conjecture.* A ratio local limit law, with the L² and bilinear routes laid out.
3. **Host operators (PROVED).** The double cone acts linearly on the path-cover polynomial:
   P_(D(B))(q) = q[P_B(q+1) + qP_B(q) − P_B(q−1)].
   - Nested double cones are explicit linear statistics of the core; for example H(D²(TT_k)) = 2·3^k + 5·2^k − 1.
   - Hosts form a semigroup of functionals, and HYP-9772 is its first coordinate.
4. **Exact repairs.** Local moves cannot repair at unit scale, and bilinear swaps are coarse. So exact placement
   must come from equidistribution, or from explicit additive digit families: carrier rows, near-joins, and wave 37's
   fibre swaps in steps of 24.

## Replay

| script | what it checks |
|---|---|
| `affine_lemma.py` | the affine lemma, exhaustive to a = 7, sampled to 10 |
| `equi_test2.py`, `residue_bias.py` | the random-divisibility law |
| `llt_exact.py` | the exact per-line fluctuation; needs S_n.txt from `paircensus.c` on `gentourng` |
| `llt_test.py 10 3 1` | the residue profile; needs uniform samples from `tiltsample 10 0.5 ...` |
| `fliptest2.c` (`./fliptest2 9 300`, `./fliptest2 10 120`) | local moves |
| `doublecone.py` | the operator identity |

Credits:
- the two-ended equation, signature rigidity and moment bounds come from the concurrent
  [TWISTED_LATTICES.md](../chanpak_cone_collision_placement_20261115/TWISTED_LATTICES.md);
- fibre swaps are wave 37;
- the tilted ensemble is THM-5127 (B);
- the cyclic join is THM-5172 (C) and THM-5222.
