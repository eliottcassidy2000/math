# Session Result

## Task Chosen

I chose a tiny independent sanity check of the Paley character tournament
values recorded in `05-knowledge/results/number_theory_tournament_atlas_s651.out`.
The stated rows are:

```text
p=7:  residues [1, 2, 4],      score_hist {3: 7},  c3=14, SCC [7],  H=189
p=11: residues [1, 3, 4, 5, 9], score_hist {5: 11}, c3=55, SCC [11], H=95095
```

## What I Did

I ran a transient Python check using the Paley edge rule on vertices modulo
`p`: orient `i -> j` exactly when `j-i` is a nonzero quadratic residue mod
`p`. I then computed:

- the score histogram directly from the adjacency matrix,
- directed 3-cycles by scanning all triples,
- SCC sizes by forward/backward reachability, and
- `H(T)` by subset dynamic programming over terminal vertices.

Tournament Analysis declaration:

- Pairwise observable: the quadratic character of `j-i mod p`.
- Switch/gauge: residue membership gives the binary tournament edge.
- Tie Hamiltonian path: not needed; Paley comparisons for `p=7,11` have no
  ties between distinct vertices.
- Fingerprints recorded: score histogram, directed 3-cycle count, SCC sizes,
  and Hamiltonian path count.

Assumption challenge:

- I kept tournament vertices as residues modulo `p` because the claim being
  checked is explicitly a finite-field character tournament claim.
- Alternate vertices such as arcs, residues classes of gaps, Fourier modes, or
  proof obligations could expose extra arithmetic structure, but would destroy
  the direct check of the stated tournament fingerprints.

## Concrete Result

The independent computation reproduced the recorded values exactly:

```text
p=7:  residues [1, 2, 4],       score_hist {3: 7},  c3=14, SCC [7],  H=189
p=11: residues [1, 3, 4, 5, 9], score_hist {5: 11}, c3=55, SCC [11], H=95095
```

Thus the Paley carrier rows in `number_theory_tournament_atlas_s651.out` pass
this independent spot-check for `p=7` and `p=11`.

## Confidence Note

Confidence is high for this bounded verification. The computation is small,
uses a direct adjacency construction from quadratic residues, and counts
Hamiltonian paths with a DP independent of the prior S651 artifact.
