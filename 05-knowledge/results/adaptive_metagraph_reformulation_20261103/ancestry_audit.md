# Independent ancestry audit

**PASS.** Read and audited [ANCESTRY.md](ANCESTRY.md), A1–A7. The all-order
proofs are sound with the stated strong-parent and ordinary-isomorphism scope.
The primitive-core inventory remains explicit and unresolved.

A1's degree argument is decisive: three indegree-one vertices would form a
dominating C3 in a larger tournament; with two, deleting their internal tail
makes their internal head a source. Consequently at most one deletion can
leave a strong parent. Isomorphisms therefore transport the newborn and its
unique predecessor. Repeating this makes each chronological newborn fixed
under every final automorphism.

A2's exact group is the core automorphism group fixing the core letters used
as roots, pointwise. Earlier newborn roots add no further core condition.
This proves the root-word orbit description in both directions. For A3,
a core automorphism g has fix(g)+i choices at birth index i; Burnside applies
to the product. Empty histories need the empty-product convention. The
triangle formula (n−1)!/6 applies only for n≥4; A3=1 separately. Nonfree core
examples confirm why dividing every root-word count by |Aut(K)| would fail.

A4 counts newborn name, root name and named parent exactly once by A1.
Dividing by n!(n−1)! makes the recurrence telescope, proving its stated
expansion. A5's tree representation retains the chronological labels and the
initial triangle gauge: erasing chronology is not an authorized quotient.
The indegree update in A7 follows directly from the birth cut.

## Independent computation

`ancestry_audit.py` imports no producer or repository modules. It builds every
labelled tournament at orders three through five, checks strongness by literal
reachability, and audits every strong source:

| Order | Strong labelled sources | Irreducible labelled cores |
|---|---:|---:|
| 3 | 2 | 2 |
| 4 | 24 | 0 |
| 5 | 544 | 64 |

For all 570 sources it replays named histories and compares the **full**
automorphism group from permutations with the predicted pointwise stabilizer.
It separately constructs every root word of lengths one and two above two
primitive five-vertex cores and canonicalizes the actual six/seven-vertex
objects by all permutations:

| Core | Aut order | Depth-one ordinary descendants | Depth-two ordinary descendants |
|---|---:|---:|---:|
| Regular R5 | 5 | 1 | 6 |
| C3[C3,1,1] | 3 | 3 | 14 |

For the second core, the descendant automorphism histograms are
{1:1,3:2} and {1:8,3:6}, respectively. Each ordinary descendant's root-word
multiplicity times its automorphism order equals the core automorphism order.
These checks retain source weights and do not treat all root words as free.

## Quotient and rewrite boundaries

The independent converse hostile has history (0,1,2), source rows
(10,20,33,6,13,27), and converse rows (52,41,26,49,34,4). Their oriented core
orders are three and six. Every C3-root history at orders four and five was
checked without such an exit. Ordinary ancestry therefore needs an explicit
additional gauge when applied to converse-merged nodes.

The separately proved [two-sided obstruction](TWO_SIDED_PEELING.md) is stronger
than this gauge warning: allowing both deletion orientations is nonconfluent
at every source order at least six, with terminal cores of different orders.
Neither boundary contradicts the one-sided ancestry theorem.

Reproduction:

    python3 ancestry_audit.py > ancestry_audit.json

Normal and optimized Python runs produce identical JSON. The finite audit
checks implementation and sharp examples; the all-order statements rest on
the proof above and in ANCESTRY.md.
