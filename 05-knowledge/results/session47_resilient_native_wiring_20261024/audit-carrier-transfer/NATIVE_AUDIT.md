# Independent native audit: complete order-four-to-six suspension

24 October 2026. **PASS:** complete finite native carrier and readout reconstruction; exact all-modulus harmonic lift. **Not claimed:** an all-size suspension theorem, any order-eight-to-ten transfer, target-570 native harmonic realization, or physical-row move preservation.

## What was independently rebuilt

The producer packet is `session47-carrier-transfer`, specifically `native_suspension_probe.py`, `native-normal.json`, and `SUSPENSION.md`. This audit does not import its code. It uses every permutation of all four or six vertices, not the producer's score-partition canonicalizer. For each actual word it computes the minimum labelled tournament code and the converse code, checks the entire class partition against the stored owner map, and authenticates every stored canonical representative's class.

Every tournament has a directed Hamiltonian path: inductively insert a new vertex before the first path vertex it beats, or append it if it beats none. Consequently normalization to a fixed Hamiltonian path reaches every tournament class. Enumerating all free nonconsecutive bits therefore covers the entire fixed-path native carrier, not an induced subset or a selected reservoir.

The independent reconstruction verifies:

- Order 4: 8 literal words, 4 ordinary classes, 2 self-converse classes, 3 converse-merged nodes, 4 complement lines, and 3 edge-support pairs.
- Order 6: 1,024 literal words, 56 ordinary classes, 12 self-converse classes, 34 converse-merged nodes, 512 complement lines, and 183 edge-support pairs.
- Both full carriers are connected. Rooted matrices retain every nonroot coordinate.

To independently count edges, this audit chooses the closed endpoint-bit representative of each complement pair, not the producer's numerical smaller-mask representative. It reconstructs every edge multiplicity and every entry of the full literal Laplacian. Self-loops are retained in the edge export and correctly contribute zero to this Laplacian.

The staircase reflection sends chord (i,j) to (n-1-j,n-1-i). A word fixed by reflection is blue. Other literal lines form black reflection pairs. The independently reconstructed color totals are:

- Order 4: 2 blue and 2 black lines; 1 blue loop, no black loops.
- Order 6: 32 blue and 480 black lines; 2 blue loops and 24 black loops.

Every merged black edge multiplicity is even. The full Laplacian uses blue plus literal black multiplicity, not blue plus half the black count. This agrees with the inherited convention G=L+2B, where B alone divides paired black lines by two. The audit checks the reflection/converse owner relation for all 1,032 words.

## Actual suspension records, not a guessed class map

For each of the eight old words, the audit independently forms S(T)=a->T->b with a->b and checks the actual Hamiltonian path (a,0,1,2,3,b). It complements all nonpath arcs and verifies that deletion of a,b gives the literal old complement tournament. Both new source and new target owners are found by the independently reconstructed full canonical class map.

The resulting source/target endpoint evaluations are exactly the nine distinct (new node, old node) readout pairs in the producer receipt. There is no collision assigning conflicting old nodes to a selected new node. This checks both sides separately and does not identify the complemented target with S(K_P(T)).

## Nonvacuous exact all-modulus certificate

The reconstructed old grounded matrix is

    A = [[2,-2],[-2,3]].

Its determinant is 2 and its entries have gcd 1, so its Smith factors are (1,2). At every even modulus m its grounded harmonic bank is exactly {(0,0),(m/2,0)}. At odd m the bank is zero. Thus the positive result includes a genuine nonzero old character.

Let N be the entire 33-by-33 new grounded Laplacian. Stack all nine endpoint constraints below it to obtain M, and let J contain the matching old readouts below zero harmonic rows. Construct Z independently: its first column is the indicator of new canonical nodes 3 and 21, and its second column is zero. Compute

    Q = (M Z - J) A^{-1},
    A^{-1} = [[3/2,1],[1,1]].

Every entry of Q is an integer, and independent exact multiplication verifies

    M Z - J = Q A.

The independently constructed Z and Q equal the producer's full saved matrices. Therefore, for every integer modulus m and every old harmonic x, Zx is harmonic on every new physical vertex and has every prescribed endpoint evaluation. This is an all-modulus theorem from an integer identity, not extrapolation from a finite sample.

The audit also rechecks all saved modulo-2/4/8 witnesses and tests the fixed Z on every old harmonic at each modulus 2 through 12. It does not claim that Z is an ambient right inverse of the readout on nonharmonic old vectors.

## Reproduction and receipt

    python3 session47-audit-carrier-transfer/audit_native_probe.py
    python3 -O session47-audit-carrier-transfer/audit_native_probe.py

The independent normal and optimized receipts are byte-identical and contain 3,147 explicit checks, the complete colored edge tables, all eight actual suspension records, the nine readouts, and the independently computed integer Z,Q. The producer native checker was also replayed separately and its receipt matched exactly. Input and output hashes are frozen in `REVIEWED_INPUTS.json` and `MANIFEST.json`.

The numerical inherited old order-eight group was checked separately against `ELECTRICAL_CODES.md` and the earlier independent filtration audit. Its Smith levels 1,2,8 yield the sufficient moduli 2,4,256 in the general transfer theorem. That dependency does not turn this order-four-to-six native certificate into an order-eight-to-ten one.

## Portable publication layout

The checker accepts `--producer-dir` pointing to the directory containing `native-normal.json`. Without the option it resolves either the published sibling `carrier-transfer/` or the original sibling `session47-carrier-transfer/`. A clean renamed publication layout was tested from an unrelated working directory, both with autodiscovery and with the explicit option; both receipts match the frozen audit exactly. `LAYOUT_CHECKS.json` records the check.
