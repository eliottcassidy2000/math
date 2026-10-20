# Validation receipt

2026-10-20. Producer-side exact tests PASS. Independent review is a separate artifact.

- `python3 verify_generator_free_naturality.py --output verification.json > verification.out`
- `python3 -O verify_generator_free_naturality.py > verification.optimized.out`
- `cmp verification.out verification.optimized.out`

Both runs completed and produced byte-identical output. There are **183,074 explicit checks**, all retained under Python optimization.

Universe: every skew word and every symmetric even-primitive-root word for m=1,...,13, quotient by all rotations. Every unit decimation modulo 2m is checked against independently computed fixed-type counts. For every odd m in this range, every input is tested for output membership, exact two-element labelled fibers, type bijection, rotation and Aut(C) equivariance. For even m≥4 the explicit moved-source/fixed-target obstruction is checked, including fixed-point-free strengthening when 4 divides m. The m=2 exception is retained.

Prime powers q=9,25,27,49,81 use exact coset Burnside sums and divisor-period recursion; the smaller ranges are independently cross-checked by exhaustive enumeration. The q=27 odd-degree control has equal nonzero Frobenius-fixed counts 16/16. The even-degree cases have no source fixed type but positive target fixed counts.

Hostile control: the q=9 even-root graph is independently constructed in F_3[t]/(t²+1) and identified with the 18-edge rook graph on nine vertices. A verified graph involution swaps the endpoints of exactly three edges, disproving full sign-evenness despite even primitive-root weight.

The first-order complexity statement is not a finite test result. Its proof imports only the classical nonuniform constant-depth parity lower bound, with primary citations in the main note.

No repository checkout files were edited by the producer. `probe.py` and `counts.py` are scratch experiments and are not needed for publication; publish the main note, final verifier, this receipt, JSON/output evidence, and SHA256SUMS only.
