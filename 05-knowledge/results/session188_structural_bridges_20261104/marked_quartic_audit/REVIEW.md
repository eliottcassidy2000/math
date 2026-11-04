# Independent marked-quartic review

Status: PASS. No mathematical correction required. Producer proof 14f9df3ba3f336f94974340e83b5f62e6d337e3b20388aba7f9c94bce761efa8 and checker d2cbf4cc743856598ae10ce8f150dc985d67384bc65c7e6a6798fa21d97082cd are pinned in both independent receipts.

## Analytical checks

The marked ordered covers use a specified path slot, not a sum over slots. Their endpoint normalization is consequently exactly q_r, without a factor r. In a quotient word for two separate path slots, the run containing a prescribed endpoint occupies one fixed position among that child's run occurrences. Permuting the ordered internal cover slots transports that position to the first slot bijectively. This validates the profile functional without either omitting or adding a multiplicity. Distinct marked children are essential: one same-child joint profile is not replaced by two independent marginals.

Quotient reversal transposes the run matrix. It sends c_i r_j to c_j r_i, and swapping the physical endpoint order swaps E_u,S_v to S_u,E_v with the opposite kernel sign. Hence the symmetric ordered-cover contrast is exactly the antisymmetric profile tensor paired with K_ij. Applying the inherited one-shared-label deletion identity gives the stated physical quartic, with positive deletion-card signs. Deleting a name preserves the quotient only when every child stays nonempty, as explicitly assumed.

The next-to-last marked coordinate has precisely one two-vertex internal path. Its marked-path contribution is (n-2)! times the incoming or outgoing degree; all singleton-mark contributions cancel. Thus endpoint balance implies regularity, and every surviving vertex on an even deletion card has a nonzero next-to-last imbalance. Vertex transitivity supplies endpoint balance by uniform endpoint counts. A homogeneous directed pair has unequal outdegrees, so a nontrivial vertex-transitive tournament is pair-free.

If the repeated name lies in neither different marked child, both marked profiles remain balanced and the tensor vanishes exactly, independent of the exterior kernel. No conclusion is inferred for same-child marks, disjoint quartics, cubic coordinates, or outer-kernel non-annihilation. In particular, a nonzero deletion-card input is not asserted to be a response separator by itself.

## Independent executable evidence

The audit imports no producer code. It constructs every marked profile by literal vertex permutations and cut sets rather than subset cover convolution. It checks all 1,099 labelled child tournaments through order five and the score formula and normalization throughout.

For 512 configurations of all three-name quotient orientations and singleton, both oriented-pair, and directed-triangle children, it compares every different-child physical endpoint factorization against an independent full-graph cover calculation. The run kernel is formed by literal bounded quotient-word enumeration, independently of the producer's count-vector recurrence.

For all four retained named examples it reconstructs physical substitutions and deleted cards, checks every marked profile and every kernel coefficient, and recovers quartics -24, 0, -170964, 0. A separate C++ predecessor dynamic program recomputes all 128 actual physical Hamilton counts. It receives independently reconstructed graph rows rather than embedding producer graph literals. The augmented orders are at most 17; the stated integer bounds are safe.

Both producer modes replay exactly; both independent modes are byte-identical. These finite checks validate the proof's conventions and examples. They do not prove general prime-quotient separation or any exact collision.
