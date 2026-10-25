# Validation and boundaries

`check.py` compares two representations: the parity-constraint compiler
in `explore.py`, and a direct Boolean adjacency matrix evaluator which
does not use its anti or K-move routine. All test gates raise explicit
exceptions and remain active under `python -O`.

Normal and optimized runs are required to have identical JSON receipts.
The result is **FINITE-EXACT**, with the following explicit universes:

- At orders2..5, every labelled tournament having the standard directed
  blue path and every second Hamiltonian word:2030 direct comparisons.
  The152 pairs of path words are compared with exact affine-family sizes.
- At order3, all8 labelled sources and all216 triples of words:1728
  literal target/legality controls, of which36 are legal ledgers.
- Every2438 source cell in the new guarded theorem at orders11..36;
  the solver proves legality and all three blackness claims simultaneously
  for **all** its free arc completions. The receipt counts the cells using
  the new overlap range separately.
  The population sum counts a tournament repeatedly if it belongs to
  more than one cell; it is not a census of distinct tournament classes.
- The441 overlap cycle tables for every even L=8..300 and j=1,2,3.
- All2137 possible second positional reversals through order10, with the
  first reversal fixed, comparing simultaneous-anti feasibility with the
  odd-order criterion. Prescribed path arcs are deliberately absent from
  this particular theorem and test.
- The14300 long-cut queries at orders15..40, with only the basic legality
  exclusions. The four recorded feasible blue final systems are genuine
  limitations, not negative connectivity witnesses.

There is no SCC, mixedness, quotient or isomorphism filter in the source
generation. Every endpoint port is retained. The equations are necessary
and sufficient for the stated affine source families. Positive controls
retain unconstrained simultaneous anti systems; hostile controls retain
an odd signed cycle, a correct edge-XOR with an illegal reversed join,
and the four feasible long-cut blue failures.

The unlimited mathematical claims depend on the explicit proofs in
`README.md`, not on extending these finite ranges. In particular neither
the14300-query experiment nor the generic compiler proves that the
current menu covers all orders. The tool avoids enumerating completions
of a cell; it does not prove a bounded number of cells suffices.
