# Validation and scope of the attachment packet

23 October 2026. **PROVED scoped / FINITE-EXACT**, as separated below.
All scripts use exact integer, rational, permutation or finite-set arithmetic.
No numerical fit is used as an infinite theorem. No scarce theorem ID has
been reserved or promoted by this packet.

## Mathematical review

- The lattice lane independently reviewed [KERNELS.md](KERNELS.md) and
  [MOVES.md](MOVES.md); see its [proof audit](lattice/AUDIT-kernels-moves.md).
  Common kernels are normal in the allowed projections, not necessarily in
  the entire cell group. Source kernels act on transporter fibres on the
  right. The nested-cell group and marking formulas are restricted to the
  stated equality-prefix blocks. The wide-move lower bound concerns the
  one-hot encoding; it is not a complexity lower bound.
- The complement lane independently reviewed the cyclic lattice theorem;
  see [its audit](complement/AUDIT-lattice.md). The any-degree elimination
  proof uses the Helly property of cosets in one cyclic group. Unary domains,
  every pairwise fill-in, the common hidden variable and the gcd witness
  multiplicity are retained. The noncyclic three-line hostile is essential.
  Extra path, electrical or nonconstant-weight guards are not silently
  absorbed into this binary-coset solver.
- The square-sum lane independently reviewed the observer obstruction,
  incidence identity and insertion laws; see
  [its audit](square_sum/AUDIT-complement.md). A graph/complement pair is not
  assigned one intrinsic tournament class by forgetting its observer.
  The positive order-averaged distribution has a different scope.
- The primary agent reviewed [the incidence flow](complement/INCIDENCE_FLOW.md)
  against the incoming connectedness proof: row normalization, detailed
  balance and the conditional-variance identity follow by expanding the
  exact two-sided incidence weights. The explicit star/empty/matching
  collision refutes lossless decoding of arbitrary source mixtures.
  The square-sum agent subsequently performed an independent proof review
  and normal/optimized replay, appended to its complement audit.
- A separate lattice agent reviewed [the simultaneous Gaussian decoder](INCOMING.md),
  including its axis, diagonal and zero boundary stabilizers. The complement
  agent independently tested all 2,401 ordered pairs in `[-3,3]^2`, finding
  325 joint D4 orbits and no six-invariant collision. Our frozen exact
  checker extends that finite control to `[-4,4]^2`; the all-real theorem
  uses the case proof. Invariant separation is not ring generation.
- The primary agent reviewed the square-sum founder proof, its exact finite
  base and the pinned recoloring witness in [the report](square_sum/REPORT.md).
  An alternative edge generator and descending-path search check the finite
  base independently inside its checker. This is a same-agent proof review
  of that lane, not an independent full proof audit. The all-founder tail
  starts at 71, while 70 is a sharp hostile. The three-site recoloring is
  minimal only with the stated outside pins; without them a one-site repair
  is possible.
- The lattice agent independently reviewed the
  [signed shared-clock compiler](square_sum/AUDIT-shared-clocks.md) and tested
  1,000 additional networks by literal vertex-tuple enumeration. It verified
  loops, cycle pins, absolute marks and distinctness filters. Its three-domain
  hostile demonstrates why filtered integer domains do not inherit cyclic
  coset Helly-2. Optional cross-component inequalities can already encode
  ordinary graph coloring, so no unrestricted polynomial counting claim is
  made for that extended input.

The [actual nine-cell hostile](lattice/NON_NORMAL_KERNEL.md) separately tests
normality, left/right transporter action and observation loss. It exhausts
all 9! permutations to recover the cell group and all 81 squared local pairs
for the source/target relations. This is a realizable tournament example,
not merely an abstract subgroup warning.

## Declared computational universes

| Executable | Universe or independent control |
|---|---|
| `check_kernels_moves.py` | Intrinsic prefix profiles through depth 4; 27 three-cell depth-labelled networks; all 81 cubed nine-cell frames through exact buckets; mixed-depth six-mark witnesses; native triangle chains through 64 cells; finite primitive-move box; literal Fibonacci summaries through depth 12 and exact depth-100 prefix queries |
| `check_joint_alignment.py` | All 6,561 ordered vector pairs from the integer square `[-4,4]^2`; literal simultaneous D4 orbits versus six invariants; boundary stabilizer and additive-loss controls |
| `lattice/check.py` | All cyclic subgroups and cosets for endpoint orders through 9; exact series compositions; cyclic coset triples through 12; 600 seeded finite networks; actual composite-cell and nonabelian controls |
| `lattice/audit_non_normal_kernel.py` | Actual nine-cell full automorphism census, nonnormal independent kernel, translated source/target relation and owner enumerator |
| `complement/check.py` | Complete graph and tournament isomorphism classes through order 7; every order of each graph isotype; merged incidence distributions and triangle histograms |
| `complement/cut_response.py` | Every graph, observer order, insertion cut and new neighborhood through four old vertices; direct triangle recount against score and cut formulas |
| `complement/check_flow.py` | Direct backward-word incidence through order 4; exact rational transition, conductance and quadratic-form identities; three-source mixture collision |
| `complement/audit_lattice.py` | Subgroups generated directly by pairs through endpoint order 9; all 163,309 triples of cyclic cosets through modulus 18; noncyclic pairwise/joint hostile |
| `square_sum/check.py` | Independent edge and ancestry enumeration through 120; all founder-copy colorings through 35; 95 declared arithmetic color recipes; pinned component repair; exact Pythagorean controls |
| `square_sum/audit_complement.py` | Independent literal adjacency-matrix orientation and canonicalization through order 5; both incidence directions and adjacent-swap delta |
| `square_sum/check_shared_clocks.py` | Exhaustive finite local and shared-port networks compared with literal start tuples; cycle pins, excluded starts and retained cross-component factors |
| `square_sum/audit_shared_clocks_independent.py` | 1,000 seeded networks checked by literal vertex-tuple enumeration, including equalities within one path, absolute marks and cross-component distinctness |

The order-seven histogram separation is **FINITE-EXACT**, not a proof of
all-order reconstruction. The colored ancestry search counts colorings,
not choices of parent forests. Generic S3 finite-state letters in the
Fibonacci adapter are not asserted to be automorphisms of a directed
triangle. The native C3-chain tests are separate.

## Reproduction and integrity

From the repository root run:

```text
python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/reproduce.py
```

The wrapper first verifies an LF-normalized SHA-256 inventory of every
Markdown, Python and JSON file in this packet. Then it runs every listed
checker normally and under `-O`, with exceptions rather than removable
assertions enforcing the checks. The packet has no external computational
inputs; linked older results remain mathematical dependencies with their
own scopes. `--write-manifest` explicitly freezes files and does not certify
them. The manifest is an integrity record, not a substitute for the proofs.

Prepublication repairs: a tuple/list receipt serialization mismatch was
normalized before default replay; a missing full-block shortcut in the
Fibonacci prefix evaluator was repaired before the depth-100 test completed.
These draft defects did not produce published mathematical claims. The
proof reviews also clarified converse-merged indicator lifting and cyclic
elimination's extra-guard boundary. No known demonstrated false theorem is
left live by this packet.

Unrestricted source-only metagraph wiring, arbitrary large-cell compact
canonicalization, the proposed all-size square-sum coloring rule and universal Collatz
coverage remain **OPEN**. A supplied-model solver is not their completion.
