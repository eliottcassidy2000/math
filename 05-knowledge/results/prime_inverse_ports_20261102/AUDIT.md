# Audit, status and reproduction

Scope: same-agent mathematical review with independent counting algorithms.
No independent person, subagent, Lean build or external novelty review is
claimed. No universal spectrum or Collatz theorem is claimed.

## Mathematical gates

| claim | status | critical gate / hostile |
|---|---|---|
| weighted singular-block valuation | PROVED | coordinate DAG allows at most two unpaid transitions; initial coordinates integral; sharp C1 C1 C0^r endpoint r=0 excluded |
| adjacent/separated prime-cycle form | PROVED | actual cyclic symmetry, uniform endpoint loads, cuts and intervening E product retained; five-block separation can remove rather than add cancellation |
| one-zero terminal h<=N/3 | PROVED for nonempty cut | disjoint legal extensions with recoverable parent; not a bound on arbitrary realizers |
| symmetric7p list for prime p<4375 | PROVED from inherited reductions + FINITE-EXACT table | two moving orbits need at least30625 paths; marked-eight module minimum32445; all marks checked, not one representative mark |
| regular-seven exchange / mod18 | PROVED | ordered run kernel with homogeneous ports; nonhomogeneous cut3 gives18 rather than14; p-residue is necessary only |
| regular cycle balance / six shared high-run responses | PROVED | literal C5=tr(A^5)/5; triangle-pair inclusion-exclusion; n=1,3 boundary; regularity essential |
| two-coordinate regular-nine exchange | PROVED + FINITE-EXACT controls | deficit-six sidecar retained;50 comparisons against full adjacency |
| all six-vertex local trap repairs | FINITE-EXACT complete | all32768 labels; common original path; no quotient-distance inference;288 three-edit hostiles |
| ordinal-context extension of repair | PROVED | outside arcs fixed; error optimality only for edits inside the six-vertex block |

The inherited prime-asymmetry,593 threshold, p-singular r-2 theorem,
singleton/D cost inverse, dyadic5/19 family and finite-bank sparsity are
credited in RECOVERY.md and the mathematical files. They are not relabelled
as new results. The known five-cycle extremal bound is cited separately in
REGULAR_BALANCE.md. No previously published repo claim was refuted in this
session; attempted extensions are recorded with their exact hostiles rather
than promoted and later retracted.

The small same-score hostile is order5: score list(1,2,2,2,3) admits H11
with C5+2M33=1 and H13 with C5+2M33=2. Thus the regular identity cannot be
extended just by fixing the degree multiset. Rows are in regular_balance.json.

## Input universe and independent paths

The frozen nauty2.9.3 inputs contain every unlabelled tournament at orders2
through8; their class counts and SHA-256 values are in GENERATOR.json.
All marks of every strong quotient are used. The strong orders3..8 have
counts1,1,6,35,353,6008, hence50782 marked strong inputs. This relies on
one trusted complete class generator, not two independent generators.

- `research.py`: forward endpoint DP for skeleton words and subset-cover
  convolution for inner modules.
- `module_audit.cpp`: separate remaining-vertex recursion, with covers
  counted by permutations and legal break sets. It agrees on every marked
  minimum and all12 prime targets; compiled at O0 and O2.
- `exchange.py`:85 small homogeneous contexts and12 complete adjacency
  prime witnesses, plus the nonhomogeneous/order9 hostiles.
- `valuation.py`:314 actual profiles,12500 valuation tests,5376 two-event
  tests,122 full adjacency comparisons and54 terminal tests. The order28
  weighted example is exact transport, not an order28 adjacency census.
- `local_escape.cpp`: full six-vertex adjacency DP, BFS, minimax Dijkstra.
- `escape_audit.py`: independent720-path/1024-completion subcube count,
  sublevel flood fill, complete relabeling orbits, retained-path verification,
  and737280 path-labelled native-complement incidences.
- `regular_balance.py`:86 regular controls (complete classes through7,
  every circulant sign choice at9/11, seeded degree-preserving samples at13/15),
  literal cycles and disjoint pairs versus skew-matrix multiplication;
  three nine-module cover polynomials and50 full adjacency exchanges.

Every Python check uses an explicit exception, so optimization does not
remove it. Normal and `-O` outputs agree. Both C++ audits agree at `-O0`
and `-O2` with `-Wall -Wextra -Werror`. REPLAY.json records the successful
release run plus the subsequent normal/-O replay of the added regular law.

## Exact reproduction

From the repository root (Python3 and a C++17 compiler required):

```sh
python3 05-knowledge/results/prime_inverse_ports_20261102/verify.py
python3 agents/check_docs.py
```

The first command uses frozen inputs and needs no nauty installation. For
a separate regeneration of the inputs, the generator command for each n is
`gentourng -q n`; its executable/version/hash and the exact stored banks are
in GENERATOR.json. The stored code convention is independently consumed by
the Python and C++ routines, with complete class totals as controls.

For the new all-order law alone:

```sh
python3 05-knowledge/results/prime_inverse_ports_20261102/regular_balance.py
python3 -O 05-knowledge/results/prime_inverse_ports_20261102/regular_balance.py
```

SOURCE_PINS.json fixes the inherited truth sources at public0d2b365955.
MANIFEST.sha256 seals this packet's files except itself and bytecode caches.
No shared checkout changes were included.
