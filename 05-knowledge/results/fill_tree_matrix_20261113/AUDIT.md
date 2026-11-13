# Proof reviews and reproducible computation boundaries

**ACCEPTED at the stated scopes. Chan–Pak and ordinary Collatz remain OPEN.**

## Proof review

The coordinator reviewed the square-zero path matrix, column-Laplacian
out-tree cofactor, endpoint coefficient, forest falling-factorial formula,
positive cycle birth factor, forced book cycles and missing-arc cut formula.
The unmarked forest formula supplies P2/P3; root-specific S2/E2/S3 need
additional marks and are not silently claimed as recovered by that scalar
polynomial. The underlying matrix-tree, OCF and marked-walk identities are
credited as existing mathematics.

The coordinator and Collatz lane reviewed the exact neutral SCC-cut
classification, its orthant quotient, literal multiplicity g_i+1, the d!
history count, uniform Pólya law, and rankwise difference formula. The
independent written review is
[collatz/MULTIPLICITY_REVIEW.md](collatz/MULTIPLICITY_REVIEW.md).
The empty order-zero convention was made explicit before freeze.

The coordinator reviewed every guard in the literal Collatz relaxation,
the dyadic-detour elimination, the forced even-K family, and the
rank/arborescence equivalence. The proof gives an unavoidable lower peak;
only the twenty saved family members are claimed to attain it. The
all-start convergence test is finite through100000, not a universal law.

The matrix and multiplicity lanes independently reviewed the contracted
arithmetic tree. Its physical sufficient condition is correctly marked
OPEN and its pure-birth form is REFUTED. The coordinator reviewed the
all-order target<=73 induced-union reduction and source conventions.
The independent matrix audit enumerates simple directed paths instead of
using endpoint DP, and constructs factor multisets instead of the
producer's binary divisor recurrence. It agrees on all target unions.

## Recovery before conclusions

The actual fixed-H multiplicity theory is inherited from the reviewed
October21 packet. The prime-source13→27,29→57,31→63 exclusions already
appear in the November2 metagraph/PRIMEGAME packet. Their recovery and
new cone interpretation are labelled accordingly. The current all-target
union also treats composite source counts25 and35. No priority claim is
made for the underlying combinatorial identities.

The most important inherited hostile remains the previous session's
all-order unit-ready hole at53084231. The new ordinary containment table
uses a different predicate and does not supersede or weaken that hole.
The earlier scalar-arrow nontransitivity47→73→117 but47 not→117 is
also retained: pairwise target witnesses do not automatically compose.

## Exact universes

| Check | Universe and independent comparison |
|---|---|
| Multiplicity and four numbers | All33,867 labelled tournaments through6; independent permutation checks of the positive atoms; integer Busch cutoff |
| Neutral growth | All1,098 attachments through parent order4; actual named histories;5,004 rational urn cells |
| Matrix and forest identities | All4,165 loopless labelled digraphs through4; Hamilton permutations, marked matrices, cofactor/forest expansion, path covers and labelled bins |
| Tournament subset fill | All1,099 labelled tournaments through5;33,866 subset cells;1,098 birth factors |
| Induced count union |148 inherited strong representatives through8,23,512 subsets; independent87,388 simple-path DFS; all37 target unions through73 |
| Arithmetic tree |65,534 allowed labels through131071, with the four exact splice edges |
| Literal Collatz | Native and odd quotient box comparisons at100,1000,10000,100000; all-start policy through100000; twenty exact-peak paths at even K<=40 |

The complete isomorphism bank is inherited, hashed and previously
reviewed. This session independently recounts its relevant subset data;
it does not claim a fresh generation of the complete bank or new
unrestricted A000568 values.

## Reproduce

From the repository root run each of:

```text
python -B 05-knowledge/results/fill_tree_matrix_20261113/tree_blueprint.py
python -B 05-knowledge/results/fill_tree_matrix_20261113/containment_edges.py
python -B 05-knowledge/results/fill_tree_matrix_20261113/matrix/check.py
python -B 05-knowledge/results/fill_tree_matrix_20261113/matrix/containment_audit.py
python -B 05-knowledge/results/fill_tree_matrix_20261113/multiplicity/check.py
python -B 05-knowledge/results/fill_tree_matrix_20261113/collatz/check.py
```

Each has a matching optimized receipt from `python -O -B`. Conditions
use explicit runtime checks, not removable assertions. The packet-level
`packet_audit.py` checks receipt agreement, UTF-8, local links and all
declared source hashes; it does not repeat proof reviews or inherited
isomorphism generation.
