# Session130 independent audit

**PASS — general symbolic claims audited; finite claims independently rebuilt.**
Source baseline: `5277b4505b6426f2381a01eb9186f0779566ce4e`.
No producer imports or data reads occur in either independent checker. No publication was performed by this auditor. Concurrent THM-4854 remains a source-qualified PREAUDIT claim, not a proved dependency.

## Results first

- **Literal n5 metagraph:** 64 tilings, 10 merged nodes, 32 lines = 8 blue + 24 black. There are 12 black twin pairs, including two loop pairs at nodes 6 and 8. The full graph is connected with integral cycle rank 23; its reflection quotient has rank 11. Loops contribute two to degree.
- **Integral symmetry:** quotient-cycle lifts plus 12 twin differences form a unimodular integral basis. Black-coordinate restriction has binary rank 10; the invariant plus anti-invariant sublattices therefore have index **1024**, with ranks 11 and 12. At orders 3 and 4 the corresponding indices are 1 and 1. The generic integral splitting and mixed-color component-incidence exact sequences are valid, including disconnected graphs.
- **Colored deletion:** each lower literal line has `2^(n-2)` lifts. Exactly two are blue precisely when the marked central board is reflection-fixed. Independently exhausted n4–n7. Every n5 raw/blue/primitive deletion vector and black parent-pair graph agrees with the colored-incidence packet.
- **Hamiltonian refinement:** independently checked `H ≡ 1+2q+4i (mod 8)`, full low-order OCF, length-refined fixed-packing parity, odd-center deletion, exact weighted cut, and fixed-cycle/mirror-arc parity. Coverage is all 598 reflection-fixed masks through order 7 and 128 deterministic order-8 masks, not an order-8 census.
- **Arithmetic:** 270 admissible scales, 17,280 binary codewords and 270 signed vectors independently verify the physical first-enabled cycle `GBDEFA`, scale gcd, integral decoder, and square-reflection-norm criterion. The primitive-increment-lattice proof and binary-only square-converse boundary are sound. Raw 2/6 periodic-state unions fail by the blue/black orbit-count obstruction; the Boolean decoration repair is typed correctly.
- **Hostiles:** the two H9 tournaments have isomorphic length/reflection-decorated conflict graph K4 but blue degrees 1 and 3. `(H,blue)` distinguishes all ten n5 nodes but fails at n6: masks 1 and 77 have `(17,5)`, and masks 50 and 43 have `(29,7)`, in distinct classes with distinct deletion profiles. The n6 audit independently enumerates **all 32,768 labelled tournaments**, giving 56 ordinary, 12 self-converse and 34 merged classes.
- **Gain boundary:** endpoint-ratio gains telescope, whereas the PRIMEGAME block triangle has exact holonomy `2025/4802`. A node codec cannot supply that missing edge gain.

## Independent methods and controls

`check_independent.py` first partitions the complete labelled tournament universe into explicit permutation orbits, then locates the fixed-path atlas inside it. This reverses the producers' primary fixed-path-first construction. Automorphism orders come from labelled orbit sizes; Hamiltonian and symmetric-path counts are separately enumerated. Reflection is computed both as converse plus vertex reversal and as a tile-coordinate permutation. Integral cycle bases use a separate spanning-forest construction and exact Bareiss determinant.

`check_structural.py` uses a different lexicographic free-edge gauge. Hamiltonian counts and ordered two-path covers use subset dynamic programming; directed odd cycles use DFS. Length-refined coefficients count literal cycle collections, retaining distinct cyclic words with identical support. Arithmetic uses direct valuations, integer square roots, and the original first-enabled fraction order.

Five mutation controls reject swapped color, a dropped literal line, counting loops only once, omission of the automorphism denominator, and omission of the converse sheet. Both programs also expose `--test-failing-gate`, which exits 1 normally and under `-O`. Normal and optimized result files agree byte for byte.

## Symbolic packet review

| Packet | Audit |
|---|---|
| `inheritance/INHERITANCE_AND_THEORY.md` | PASS: multivariate fixed-packing quotient; optional one fixed-center cycle; hostile; normalized cut composition. OCF and historical prior claims remain inherited/source-qualified, not independently reproved literature. |
| `colored_incidence/REPORT.md` | PASS: all-order extension count, n5 independence, raw incidence factorization, quotient failures, n5/n6 classifier boundary. |
| `primegame_lines/THEORY.md` | PASS: all seven scoped theorems; allowed-scale unipotent inverse; primitive square converse; raw-cycle obstruction; phase-shift non-descent. |
| `hamiltonian/PROOFS.md` | PASS: mod8 mirror-overlap formula, correct cut direction, division to mod4, automorphism normalization and strongest hostile. |
| `gain_cohomology/THEORY.md` | PASS: integral splitting, elementary-abelian parity quotient, mixed-color exact sequence and holonomy boundary. |

One wording repair was requested and accepted: in the gain proof the black-restriction kernel is the blue-supported subspace **on the binary cycle space**. The integral parity map also kills doubled cycles, so its kernel should not be called merely blue-supported. No theorem changed. Loop-twin orientations must be transported by reflection rather than chosen independently by smaller mask.

The finite tests do not imply all-order injectivity of profiles or `(H,blue)`, connectivity at every order, or a PRIMEGAME dynamical model of merged nodes. The proof scopes in the packets retain these boundaries.

## Reproduction

```sh
python3 check_independent.py > normal.json
python3 -O check_independent.py > optimized.json
cmp normal.json optimized.json
python3 check_structural.py > structural.json
python3 -O check_structural.py > structural.optimized.json
cmp structural.json structural.optimized.json
# Each of the following must fail with exit 1:
python3 check_independent.py --test-failing-gate
python3 -O check_independent.py --test-failing-gate
python3 check_structural.py --test-failing-gate
python3 -O check_structural.py --test-failing-gate
sha256sum -c SHA256SUMS
```
