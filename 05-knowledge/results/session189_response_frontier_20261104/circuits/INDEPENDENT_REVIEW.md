# Independent circuit proof review

**Status: approved for the stated shared-DAG model.** I reviewed the complete proof, checked the relevant verifier implementations, and independently replayed `verify.py` normally and with `python3 -O`. The two replay outputs were byte-identical and reported PASS.

## Mathematical checks

1. The exact loss is `beta(G)-beta(G-C)=E+a-kappa+1`. Connectivity follows from the live single-output assumption. Every outside component except the output component must have a directed edge into C; this proves the correction is nonnegative and gives exactly the stated equality condition.
2. The cycle-space quotient retains every physical incident edge, including parallel edges after contracting outside components. Its kernel is the direct sum of the outside component cycle spaces. The four-gate example tests precisely the multiplicity loss that would invalidate a simple-graph quotient.
3. Parent-pin aliases admit the stated branch-set minor certificate. After deleting C, later simplifications introduce no vertices and remain minors of G-C. This justifies telescoping the full Q charge, rather than only the outgoing fanout charge E. Arbitrary nonparent aliases, resynthesis, and newly introduced internal computations are correctly excluded.
4. The capacity identity follows from the C-parent indegree partition into zero, one and two and the inherited boundary identity. Both correction terms lambda and kappa-1 are nonnegative. Q therefore decomposes existing certified deletion rather than adding deletion credit. The displayed hostile has two deletions, Q=1, and no valid third deletion under the stated cleanup.
5. The shielding construction transfers p's fanout excess to the last ladder gate while leaving all other excesses unchanged. Reverse-order restriction deletes exactly the designated two gates and restores the shorter ladder. Any forbidden new constant or direct-pin literal at an old core gate would persist after setting remaining controls to zero and contradict the original core hypothesis. Fresh independent controls protect the remaining ladder gates.
6. In the XOR grid, the two parents of V(i,j) have distinguishable extremal boundary coefficients L_i and T_j. Thus each is nonzero and their XOR is neither constant nor either parent or its negation. The graph counts and grid minor are valid at every size. Minimum circuit size and better restriction schedules remain expressly outside scope.

## Computational and scope checks

The implementation counts distinct gate edges, not repeated pins; preserves isolated outside vertices when computing kappa; uses explicit failures under optimization; and checks branch-set connectivity and crossing edges. The truth-table, graph-topology and shielding checks have the declared finite scopes and do not substitute for the all-size proofs. The cited classical cycle-rank formula is used only for its standard mathematical content; this review makes no separate historical-priority claim.

No repair was needed. Source SHA-256 pins:

- PROOF.md: 5f0a750f12740cc519fe58c7101ae48bb81299711e7e1e795bdb73de316b70bb
- verify.py: 34b524bcf0c8a702ed124ce5b5bf38fb67a76873736de5150e04d43c513691b5
- RECEIPT.json: c2e272c2ae530ac568b9c9f0de0122ffd0d7f220ebae6a1dd0a24e122fcd362d
