# Session232 subgroup-spectrum verification receipt

Status: **FINITE-EXACT controls passed; all-order claims are proved from the inherited inputs explicitly named in PROOF.md.** The inherited THM-4812 all-order finite certificates were not rerun in this packet.

## Deliverables

- `PROOF.md`: all-subgroup omission of7/21, low-quotient prime cost, and the full equality classification.
- `SHARP_STRATA.md`: no intermediate cost between b and3b; exact low seven/five/twenty-five symmetry strata; optional stronger distinguished-module interaction bounds.
- `verify.py`, `receipt.json`, `normal.out`, `optimized.out`: exhaustive n≤7 controls, independent uncoloured graph-encoding controls, and explicit order10 sharp examples.
- `verify_strata.py`, `strata_receipt.json`, `strata_normal.out`, `strata_optimized.out`: threefold separation, sharp strata and nonempty ordered two-path-cover constants.

Copy those files together. Neither verifier relies on the repository layout, temporary absolute paths, network access or third-party packages. The second imports the first from its own directory. Generated `__pycache__` files are not part of the deliverable.

## Reproduction

    python3 verify.py > normal.out
    python3 -O verify.py > optimized.out
    cmp normal.out optimized.out
    python3 verify_strata.py > strata_normal.out
    python3 -O verify_strata.py > strata_optimized.out
    cmp strata_normal.out strata_optimized.out

Both comparisons passed. Python bytecode compilation passed. Every executable gate uses an explicit exception through `require`; optimized Python therefore retains it.

## Exact universes and totals

Orders1,…,7 contain respectively1,1,2,4,12,56,456 ordinary tournament classes:532 in total. The verifier checks every literal automorphism subgroup of each chosen representative:640 class/subgroup pairs. It independently counts the literal Hamilton words and their literal subgroup orbits, and verifies them against subset DP and the restriction-index formula.

There are168 pairs with q_A<23. Every prime-cost inequality, every equality iff, and every threefold-separation gate passes. The prime-cost equality occurs33 times. The sharp five-stratum has1,2,3 instances at orders5,6,7; the sharp seven-stratum has one instance at order7. The separate two-R5 control has H225, Aut25 and quotient9.

The ordered nonempty two-path-cover counts are90 for R5 and1638 for P7, independently reproduced by subset products and vertex permutations with a nontrivial cut. The C3 substitutions with those modules and two singleton modules have (H,Aut,q)=(135,5,27) and(2205,21,105).

The parallel graph worker's smaller gadget is independently checked on Q7 and P7 by a generic full uncoloured automorphism enumerator without gadget-role colors. Each graph has112 vertices,126 edges,42 leaves. Group orders are7 and21, all graph edge-reversal signs are positive, and the admissible traversal counts are175 and189. Ordinary graph Hamilton-path counts are0 because each graph has more than two leaves.

## Output hashes

`normal.out` and `optimized.out`:

    74cc51c8073b6897b64f1819f246c55ddf3b130ba464d8dbf22b15065d443aa0

`strata_normal.out` and `strata_optimized.out`:

    ce523df106b36cd407f49123ffea983dd7dc6149172886719fe311a475eb9923

Each JSON receipt also contains the SHA256 of the code that produced it. The supplement receipt records both verifier hashes. The normal and optimized files have identical bytes, including hashes.
