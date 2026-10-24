# Square-sum G27 graph-birth response

Theorem and scope: THEOREM.md. All paths and deletion witnesses: results.json.

Key finite result: G27 has3 endpoint extensions of G26 Hamiltonian paths and32 internal joins of distinct compatible two-path covers, totaling35. Its neighbors are9,22. The32 compatible covers are among1116 total spanning two-path covers of G26.

General result: exact linear-forest generating-polynomial birth law, minimum-cover trichotomy, and sufficient endpoint-pair-state transfer. These are elementary identities with no priority claim; the degree-two deletion correspondence is inherited from Session32.

Verification:

    python3 check.py > results.json
    python3 -O check.py > results.optimized.json
    cmp results.json results.optimized.json
    python3 test_rejections.py > test.out
    python3 -O test_rejections.py > test.optimized.out
    cmp test.out test.optimized.out

No external packages. Normal and optimized output agree. The separate vertex-path and selected-edge-cover searches compare full compatible cover sets. General-state tests cover1099 labeled simple graphs with33646 linear forests. SOURCE.json pins the internal mathematical inheritance; no live external literature or all-order existence claim is needed for these proofs.
