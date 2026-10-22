# Independent audit: colored rooted-forest profiles

Status: **ACCEPTED / PROVED**, with **FINITE-EXACT** independent controls on
the explicitly stated small orders. No mathematical counterexample or
unresolved proof gap was found. This is a mathematical and implementation
audit, not a literature-priority determination.

Audited submission: `BLUE_FOREST_CONSTRAINTS.md`, `check_blue_forests.py`,
`check_wiring.py`, and their helper `probe.py` in the blue-forests lane.
Exact source hashes are in `independent_receipt.json`. No Git mutation was
performed in this audit lane.

## 1. Conventions and inherited canon

The forward fixed-path convention is `1 -> ... -> n`; the free edge `(i,j)`
is sent to `(n+1-j,n+1-i)` **without negating its bit**. This is precisely
the `sigma` convention of
`01-canon/theorems/THM-643-blue-black-line-parity-structure.md`, the
`(sigma,tau)` frame and T5. It is converse followed by reversing labels,
not ordinary label reversal alone. The note does not inherit the refuted
pure-blue-leaf or blue-fiber-cap extrapolations in that theorem's historical
sections.

With `D=(n-1)(n-2)/2`, the reflection has `f=floor((n-1)/2)` fixed free
edges and `o=(D+f)/2=floor((n-1)^2/4)` orbits. Hence the presentation masses
are `M=2^D`, `B=2^o`. For `n>=3`, complement has no fixed point and
complement-reflection has none because the free pair `(1,n)` is fixed by
reflection but toggled by complement. Therefore Klein orbits have sizes
two on blue and four on black. The three counts are correctly distinguished:

- Blue complement lines: `B/2`.
- Literal black complement lines: `(M-B)/2`.
- Fundamental paired-black edges: `(M-B)/4`.

The last convention agrees with
`01-canon/theorems/THM-849-n8-black-self-line-obstruction.md`, section 2,
and `01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md`,
statement 3. The exceptional orders `n=2` and `n=3` are properly excluded
from line and conditional-black formulas where needed.

The following are **CITED prior canon, not new results**:

- The cyclic-triangle/energy identity and endpoint flux:
  `01-canon/theorems/THM-785-cyclic-triangle-flow-on-blue-black-metagraph.md`,
  sections 1-3.
- The blue parity and leg laws:
  `01-canon/theorems/THM-790-blue-parity-law-proved.md`, sections 1-4.
- Odd-order absence of blue loops, and the necessary but insufficient
  even-order endpoint-score condition:
  `01-canon/theorems/THM-648-blue-selfloops-only-at-even-n.md`.

THM-648 uses the opposite path labeling: its `s(1)=(n-2)/2` becomes the
current `s(n)=n/2-1`, and hence `s(1)=n/2`, after reversal. There is no
one-unit contradiction. The submission explicitly preserves these
inheritance labels and distinguishes the colored forest refinement.

## 2. Independent energy and degree-two derivation

Write `v_ij=e_i-e_j`, `b=e_1-e_n`. Path arcs telescope, giving
`d=b+sum x_ij v_ij`. Independent mean-zero free signs imply

    E_all W = ||b||^2 + sum ||v_ij||^2 = 2+2D = n^2-3n+4.

On blue presentations one sign is shared by every reflection orbit, so
replace the individual incidence vectors by their orbit sums. Two disjoint
edges retain squared norm four. The only overlapping two-edge orbit is
`(i,c),(c,n+1-i)` at odd `n`, with `c=(n+1)/2`; its vectors have inner
product minus one, reducing the norm contribution by two. There are
`(n-3)/2` such free-edge orbits. Thus `E_blue W=mu-delta`, with the stated
`delta`, and the black conditional mean follows by subtracting the blue
mass from the total. No uniform-isomorphism-class averaging is used.

For `L=diag(s)-A`, the second determinant coefficient is the sum of its
two-by-two principal minors. Each is `s_i s_j-A_ij A_ji=s_i s_j`, since a
tournament has no digon. Substitution of `d_i=2s_i-(n-1)` gives

    F_2 = (4N^2-n(n-1)^2-W)/8.

The source-count proof of
`c3=binom(n,3)-sum_i binom(s_i,2)` gives the claimed affine relation
between `F_2` and `c3`. Therefore the degree-two coordinate is explicitly
an inherited triangle invariant; its color-resolved mean is the additional
calculation. The formula for `q_n` and both colored shifts check exactly.

## 3. Full compatible-orbit formula

The row/outdegree Laplacian counts forests with at most one outgoing arc
per vertex, whose arcs point toward their roots. This direction convention
matters. In a fixed forest, backward path arcs are impossible. Every used
free arc imposes a bit on its reflection orbit; conflicting imposed bits
make the forest impossible on the blue cube. If there are no conflicts,
exactly `q(F)` independent orbit bits are fixed, leaving `2^(o-q(F))`
blue presentations containing the forest. Summing over forests proves both
the mean and total formulas coefficientwise.

This also explains why replacing blue arcs by their marginal probabilities
does not work: touching an orbit twice is not the same as touching two
independent bits. At `n=5`, the mistaken marginal substitution misses
`1/4` in the degree-two coefficient. The all-presentation law correctly
uses the number of free *edges*, rather than touched reflection orbits.

Independent control: the audit assigns one outgoing neighbor or a root to
each vertex, rejecting an undirected cycle by union-find as soon as it
appears. In an outdegree-at-most-one digraph this is exactly the rooted-forest
condition. It then counts all and blue occurrence totals from edge/orbit
constraints, without using a determinant. These totals agree coefficient by
coefficient with an independently implemented principal-minor determinant
calculation for every `n=2,...,6`.

## 4. Degree-three proof and first even-order mismatch

For a blue tournament, reversal satisfies `RAR=A^T` and `Rd=-d`.
Consequently `K=diag(d)-(A-A^T)` satisfies `K1=0` and `RKR=-K`.
The constant line is invariant under both operators, so they descend to the
quotient. On that `(n-1)`-dimensional space, `L=nI/2+K/2` because the
all-ones matrix vanishes in the quotient. Similarity of `K` to `-K` forces
the quotient determinant to have only even powers of its `K` parameter.
No diagonalizability assumption is needed.

Also `tr K=0` and `tr K^2=W-n(n-1)`, so the first nonconstant quotient term
is `c_1=(n(n-1)-W)/8`. Its blue mean is `(n-2)/4+delta/8`. In the previously
derived all-presentation matching formula, the corresponding coefficient
is `(n-2)/4`. Expanding the factor `(1+nt/2)^(n-3)` therefore gives exactly

    [t^3](Q_n-P_n) = delta*n*(n-3)/16.

Higher paired terms begin at degree four, so they cannot affect this
argument. Exact independent cube averages give

    Q_5-P_5 = t^2/4 + 5t^3/4 + 27t^4/16,
    Q_6-P_6 = 5t^4/16 + 15t^5/16.

The second line is an explicit hostile to the implication “equal colored
score-energy means imply equal full forest means.” It does not prove an
all-order sign or positivity law, and the submission claims none.

## 5. Merged inventories and witness

On a distinct converse pair, the two ordinary fixed-path fibers have equal
size: reversing a Hamiltonian path gives a bijection with paths in the
converse, and the automorphism groups agree. Thus their combined mass is
`2h`; their full polynomial contribution is `2h Fbar`, not `2h F` for an
arbitrary representative. On a self-converse class the forest polynomial
is unchanged by converse because the tournaments are isomorphic.

Dividing black presentation mass by two gives fundamental-black **degree**
mass, not edge mass; loops count twice. This proves both full-polynomial
inventory equations, and the energy/degree-two equations follow by taking
their relevant coefficients. Independent ordinary-class and converse-pair
enumeration through order five checks the identities. It also finds
non-self-converse pairs with different full forest polynomials already at
order four, confirming that the averaging instruction is necessary.

The explicit six-vertex adjacency list in the submission was entered
literally, without using its claimed orbit mask. It is a tournament,
contains the fixed path, and is reflection-fixed. Its free-edge complement
is likewise blue. Both have score multiset `{1,2,2,3,3,4}`, energy `22`,
and six cyclic triangles, but their forest polynomials are respectively

    (1,15,91,279,433,273),
    (1,15,91,279,432,270).

These differ by `t^4+3t^5`. Distinct forest polynomials prove
nonisomorphism, and both endpoints are self-converse via label reversal.
An additional independent permutation test found neither an isomorphism
nor a converse-isomorphism between the endpoints. Hence they remain
distinct on the converse-merged graph.

Minimality is accepted in the stated nontrivial complement-line setting
`n>=3`: odd orders cannot have zero blue flux by the inherited parity law;
at order four both zero-flux blue presentations are loop endpoints. At
order six there are 24 zero-flux blue presentations and only four loop
endpoints. These last two pools were independently checked by permutations,
not inferred from the inherited census.

## 6. Implementation review and scope

The submission's Newton-identity polynomial routine is standard and
correct; its exact integer divisibility assertion is appropriate. Its
Bareiss evaluation checks supply a second determinant route. This audit
instead uses signed column-subset dynamic programming for every principal
minor, importing no source-worker module.

One non-mathematical test-quality nit was reported to the author and
**repaired before finalization**: the original `check_wiring.py` accumulated
a `blue` Boolean without asserting or returning it, and inferred the blue
mode from `len(orbs)<D`, which skipped that mode's endpoint check at `n=3`.
The final implementation has an explicit `expected_blue` flag, passes it
for every blue run including `n=3`, checks that route's endpoint identity,
and asserts the computed reflection flag. The change was inspected and
the independent checker rerun against its final source hashes. No proof,
numerical claim, or enumeration bound changed. This audit also explicitly
verifies reflection and blue endpoint identities in every blue presentation
of the small cube.

Reproduction from the working root:

    python3 session31-blue-audit/independent_checks.py

The checker also accepts `--source-dir PATH` and recognizes sibling
`blue-forests` in a published packet. It requires only Python's standard
library. Do not use `-O`; the checker deliberately refuses to run with its
assertions disabled. The receipt records exact source hashes, computed
polynomial totals, the literal witness, and all tested ranges.

Independent finite controls cover all fixed-path presentations through
`n=6` (at most 1,024), incidence-vector formulas through `n=31` without
presentation enumeration, rooted-forest occurrence sums through `n=6`,
and merged-node inventories through `n=5`. The submission's larger
`n=7,8,9` rows were not independently recomputed here. General claims rest
on the proofs above, not extrapolation from these tests. No claim of
universal coefficient monotonicity, complete graph wiring, or a new closed
product for the blue forest mean is accepted or made.
