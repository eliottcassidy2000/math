# Exact unsigned-curl realization phases: independent audit

Verdict: PASS for `session36-context/REALIZATION_PHASES.md`. The encoding is elementary Boolean constraint compilation; its stated scope does not claim new general SAT machinery, polynomial-time recognition, or NP-hardness.

## Signs and bijection

For i<j<k, the cyclic-direction literals are (x_ij,x_jk,1−x_ik), with x_ij=1 meaning i→j. The triple is cyclic exactly when those three literals are equal. Positive triples therefore impose x_ij XOR x_jk=0 and x_ij XOR x_ik=1, while negative triples impose their signed not-all-equal clause.

A consistent positive system partitions edge variables into affine parity components, one free bit per component. Substitution into the negative clauses gives an exact bijection between satisfying phase assignments and all labelled tournament realizations. The long-edge negation is correct and necessary.

Repeated-literal simplifications and global complement pairing are correctly stated, including the one-vertex exception. The imported C3-structure theorem classifies the full fibres but is not being misrepresented as a recognition algorithm.

## Two obstructions and the phase gap

The four-vertex three-positive-triangle example is already XOR inconsistent. In contrast, the five-vertex target {014,023} leaves six positive phase bits, yet is unrealizable. The normalized-cycle contradiction is valid: avoiding013,024,124 forces3→1,4→2,1→2 and then produces forbidden123. All eight negative clauses genuinely involve three distinct positive-phase variables.

This is the smallest order at which positive consistency does not suffice, as established by the exhaustive input census. Zero or one demanded positive triple is always realizable by an ordinal sum of singleton blocks and, when needed, the prescribed cyclic triple.

For the target {012,014,023}, the displayed edge expressions are correct. Negative clauses force d=1−a, then b=a, then c=1−a, giving exactly masks117 and906. Their indecomposability is independently checked. Thus four positive phase components need not mean four surviving orientation bits.

## Independent finite certificate

`independent_realization_phases.py` uses exhaustive orientation bitsets and GF(2) Gaussian elimination, not parity union-find and not producer imports. It verifies all candidate cyclic-triple hypergraphs through order5, including impossible inputs. The order-five split is exactly:

-771 positive-XOR-inconsistent;
-45 positive-consistent but residual-NAE-unsatisfiable;
-208 realizable.

All stated phase histograms, indecomposable-fibre histograms, and converse pairings are independently reproduced. The literal convention is checked on10,504 actual triple instances. Target mask12 has6 phases and no solution; target mask13 has4 phases and exactly117/906.

Outputs are `independent_realization_phases.out` and `realization_phase_receipt.json`. Normal and optimized outputs are byte-identical.
