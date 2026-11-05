# Independent review: pair-entry response barrier

**PASS.** No mathematical or source correction was required. The scope is an exact rooted-state obstruction for the specified ordinary one-seam move, not an obstruction to its scalar Hamilton-path counts.

## Frozen source reviewed

All paths are under `session197-hybrid-constructions/pair-entry-barrier/`:

- PROOF.md: d69c4597b5ed1b4de2b850293320efe3c2261a571f162d7a1dbe0a206fb33d17
- core.py: dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042
- verify.py: 8ec820b69c54bb445b453cc3271b4c1296d4743b852889f7f8187392875269ae
- control_graphs.json: d4603d925bf5d647d98fd6e828ea4d64189e09b762ba8a7aceaaecda60f227d9
- VALIDATION.json and VALIDATION_OPTIMIZED.json: 45d50632038949e643d4aed9f3277f4384816d3e023f59ddec64fdcc23012cde

The input files were copied before replay and left unchanged. Both ordinary and optimized Python runs reproduce the frozen receipt exactly.

## Mathematical checks

The universal cone follows for arbitrary nonempty physical inputs, not only the rooted C1 grammar. I checked the separate L=1, positive even L, and odd L>=3 cases. The gadget inequalities 2g-Q>=1 (and >=3 at odd L>=3), B<=2L-2, R>=h, and e<=h give the claimed nonpositive defect. No assumption that the input already satisfies the output cone is used; the singleton is therefore harmless.

For the adjacent outdegree-one vertex a with sole successor b, the two-cover split proving R_b=2(h-e_a)-e_b is exact. The forbidden reserved singleton b accounts for precisely H(T-{a,b})=e_b. The first clone is the newly born f with f->a->b; the old a is the last clone. Substituting the independently reviewed marked profiles gives exactly R_a-h+e_b, which is positive. Re-marking the same graph does not preserve the response tuple, and the proof correctly makes no scalar exclusion.

The inside-cone profile (13,2,18) is physical. For L>=2, its ending count would force old H=2; for L=1, the unique inverse is (7,5,11). The stated parity/SCC-product/strong-order-bound reduction excludes these at every order. The small strong-order census is only the final base check after that reduction, not a finite-to-infinite inference.

## Independent controls

I rebuilt Hamilton paths and ordered two-covers by literal vertex permutations and cut positions, without importing the source graph DP. All 5,405 rooted labelled tournament controls through order five satisfy the physical input inequalities. Every one of the 654 adjacent-clone cases in the labelled C1 inputs through order five satisfies the exact identity, positive first-clone excess, and last-clone cone. The observed excesses are 1,2,5. The literal marked H13 profile and the complete strong-order-three/four path sets {3},{5} also agree.

The independent source and receipts are in this review directory:

- independent.py: 1cd7dfeae8b71b3f697d6de0bf4df211898ce3a7020869d16f6da999bafeab71
- independent.json and independent-optimized.json: c0397a36fd79779ad27c940226c9ee9a11dd5cb6d6b9f44c64d958baee3f4e36

Normal and optimized independent runs agree. No unresolved issue remains in the stated theorem.
