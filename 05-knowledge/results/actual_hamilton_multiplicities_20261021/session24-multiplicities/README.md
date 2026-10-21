# Actual-H multiplicities and the two-degree gap defect

This packet keeps actual Hamiltonian-path count H separate from presentation count h=H/Aut and ordinary classes separate from converse-merged classes.

Main result: fixed-H ordinary multiplicities are eventually polynomials with exact degree equal to the maximum number of nontrivial strong SCCs. Their onset is at most floor(3 log_5 H); the symmetry/merged parity quasipolynomial onset is at most one more. The two forbidden values7,21 force an exact two-degree deficit along H=7*3^a. The leading constants and the contrasting typical metagraph degrees1 versus63 are explicit. See PROOF.md.

Finite-exact data: all classes through n8, including complete global strong-atom inventories for odd H<=73. n8 contains6880 ordinary classes,176 self-converse classes and3528 converse-merged classes. All320 weighted H histogram bins match the independent frozen full2^28 labelled census. multiplicities_n3_n8.csv contains the ordinary/SC/merged/labelled channels. low_H_formulas.json contains exact rational functions and all-order ordinary polynomials for H<=73.

Inheritance and novelty boundaries are explicit in the proof. Ordered SCC words and square-divisor centers are inherited; no literature-priority claim or proof that7,21 are the only missing odds is made. Independent review is required before promotion.

## Replay

Use Python3 with SymPy and a C++17 compiler. Set MATH_REPO to the repository containing the frozen source atlases if the packet is not inside or next to that repository.

    export MATH_REPO=/path/to/math
    python3 extract.py > /tmp/multiplicity-extraction.txt
    g++ -O2 -std=c++17 augment8.cpp -o /tmp/multiplicity-augment8
    /tmp/multiplicity-augment8 n7_representatives.txt > /tmp/multiplicity-n8.csv
    cmp /tmp/multiplicity-n8.csv n8_classes.csv
    python3 analyze.py > /tmp/multiplicity-analysis.txt
    python3 check.py > /tmp/multiplicity-check.json
    python3 -O check.py > /tmp/multiplicity-check-O.json
    cmp /tmp/multiplicity-check.json receipt.json
    cmp /tmp/multiplicity-check-O.json optimized-receipt.json

The active independent checker performs2063 checks even under python -O. The producer scripts use assertions in their ordinary replay. The stored n8 classes encode each8-vertex adjacency matrix in a64-bit integer with bit8*i+j indicating i->j. All arithmetic is exact; no floating-point fit is used.

Sources are fingerprinted in census.json. The n8 generator augments456 ordinary n7 representatives by128 new-vertex orientations. Every n8 tournament has an n7 deletion, so this is exhaustive; canonicalization enumerates every degree-preserving relabelling. The independent frozen labelled histogram tests every H bin after weighting by8!/Aut. Source files were read only, this packet was staged outside the repository, and no Git action was performed by its producer.
