# Independent alphabet-quotient and information-bound audit

October 23, 2026. This addendum is separate from the original order-eight
wiring and electrical audit. It does not alter their mathematical claims.

## Exact ternary forest and binary boundary

`audit_ternary_quotients.py` reads only the inherited, independently certified
list of all 6,880 canonical order-eight classes. It reconstructs every
tournament, counts cyclic triangles directly, and constructs both forward
and converse signatures by reversing the actual orientation. It does not
import a producer checker, use producer converse IDs, or consult the
quotient-search transcript.

For the color rule

    chi_T(u -> v) = (c_v + c_uv) mod 3,

it independently finds 6,864 ordinary signatures and 6,880 distinct ordered
forward/converse signature pairs. The bipartite signature-incidence
multigraph has 13,728 vertices, 6,880 edges, and 6,848 components:

* 6,820 isolated edges;
* 24 paths of length two;
* four paths of length three.

Its maximum degree is two, its cycle rank is zero, and simultaneous leaf
peeling removes 6,876 classes followed by the remaining four. Those four
middle-edge classes have canonical codes 8393858, 8654866, 25694466, and
34621760. Thus the singleton-or-converse-singleton shortcut is insufficient
here, but the full forest argument succeeds.

For completeness, the independent audit directly evaluates every one of
the 255 partitions of the nine E33 residue symbols into two nonempty
unlabeled blocks. No pruning is used. None produces a forest; their cycle
ranks range from 10 to 3,591. The single-color quotient has cycle rank
6,713. The 255 binary results, including an explicit cycle-closing class
and full graph counts for each partition, are retained in the receipt.

Consequently three is the exact minimum for the **forest criterion among
all set-partition quotients of the fixed nine E33 symbols**. This argument
alone does not imply that binary quotients fail the actual constrained
port-table problem. A cycle can be obstructed by capacities or support.
The full producer search's 117 ternary successes are not independently
recomputed by this audit: the selected success and every smaller quotient
are sufficient to establish the stated minimum.

### Stronger boundary: actual literal port counterfeits for all 255 binary quotients

The separate `audit_binary_literal_switches.py` closes the capacity gap.
It checks all 255 explicit four-port matching certificates against the
independent FIRST replay's complete raw mask map, ordinary class ledger,
native port ledger, and directed merged ledger. The producer supplies only
the proposed witness matching; fresh direct cyclic-triangle counts supply
the quotient signatures and triangle potential.

Every witness passes all of the following checks:

* four distinct actual closed ports map bijectively to exactly their four
  native open ports, with endpoint and black-color compatibility;
* each target signature agrees with its native target, including at both
  reflected partners, and the literal bijection commutes with reflection;
* all native source endpoint/color/signature margins and all actual target
  endpoint/color capacities are unchanged;
* the induced ordinary-cell delta agrees exactly with the reported delta,
  every resulting count is nonnegative, and all reflection equalities hold;
* literal black twin multiplicities remain even, reflection-fixed ordinary
  black cells remain even, and every blue entry and loop is untouched;
* all degrees and cyclic-triangle currents are unchanged, while both the
  directed and undirected merged wiring projections genuinely change.

There are 61 distinct literal switches across the 255 binary quotients.
The full certificate is pinned at SHA256
`9ac70b69fba35e90b36eae1e417e31fdcd989088dddd6a4ea80fbff93f0c87df`.
Any binary witness survives the final coarsening to one color.

Together with the ternary forest, this proves the stronger exact minimum:
**three colors among all set-partition quotients of the fixed nine E33
symbols for the specified reflected endpoint/color/signature port-matching
model**. It is no longer only a forest-certification minimum. It is not a
minimum over arbitrary color rules, an all-order statement, or evidence
that the alternate matchings arise from another canonical-complement
tournament realization. They are genuine alternate matchings of actual
ports satisfying the declared observations.

### Why cycle failure survives quotienting

A color quotient merges signature vertices within each of the two shores,
but preserves the individual ordinary-class edges, including parallel
edges. If B is the oriented incidence matrix, the new matrix is RB for a
row-merging matrix R. Every nonzero vector in ker B remains in ker(RB).
Therefore a nonforest cannot become a forest under further quotienting.
Equivalently, merging vertices in the same component increases cycle rank
by one, while merging vertices from distinct components leaves it fixed.

Every successful partition consequently has only successful refinements.
A chain of single-block merges from the identity reaches it through
successful parents. This independently validates the producer's pruning
principle; parallel class edges must never be deleted.

## Aligned local decoder, not anonymous reconstruction

Write r_v and r*_v for the weighted sums of the forward and converse color
histograms at the **same aligned vertex** v. Cyclic-triangle incidence gives

    sum_(v -> w) c_vw = c_v = sum_(w -> v) c_vw,
    sum_v c_v = 3 c3(T).

Hence r=(I+A)c and r*=(I+A^t)c. Since A+A^t=J-I,

    r+r* = (I+J)c = c  mod 3.

If each actual edge's color and head identity are also available, then
c_uv=chi_T(u->v)-c_v mod 3. The sign convention for the paired-edge
gradient is explicit:

    chi_T(u->v) - chi_Top(v->u) = c_v-c_u mod 3.

The reversed subtraction has the negative gradient. The independent
checker verifies these identities with 440,320 vertex/edge checks across
all order-eight classes. The displayed proof establishes every order.

For the general linear palette alpha*c_v+beta*c_uv, the aligned sum is
(2 beta-alpha)c mod 3. Thus this particular sum decoder is invertible
exactly when 2 beta-alpha is nonzero. The plus palette decodes; the minus
palette's sum vanishes. An edge-only palette can also decode aligned c,
so aligned decoding alone cannot explain anonymous completeness.

Neither independent outer multiset provides the vertex alignment needed
above, and a histogram does not identify the head of each colored edge.
The identities do not recover those missing associations, establish
all-size class separation, or certify realizability of arbitrary colored
arrays. The gradient lives on tournament vertices, not metagraph classes.

## Growing-alphabet information obstruction

The independently reviewed model assigns L colors to arcs, retains at
each vertex its outgoing color histogram, and retains only the multiset
of those n histograms. There is no additional side information.

Every signature has total histogram mass N=n(n-1)/2. Ordering its rows
overcounts the possible outputs, giving at most

    binom(N+nL-1,nL-1)

signatures. Conversely, fixing one directed spanning cycle gives
2^(N-n) labeled strong tournaments, hence at least 2^(N-n)/n! strong
ordinary classes. Pigeonhole comparisons between these bounds are valid
even if the coloring rule itself is arbitrarily sophisticated.

Let f(x)=(1+x)log_2(1+x)-x log_2 x. The binomial entropy bound implies
that eventual separation of all strong classes requires

    liminf_(n -> infinity) L(n)/n >= alpha,
    f(2 alpha)=1,
    alpha = 0.146907686670207746658300826652... .

The liminf conclusion does not assume that L/n converges: any subsequence
bounded above by a fixed constant less than alpha has a fixed entropy
gap against the strong-class lower bound. This is necessary, not sufficient.

The explicit uniform obstruction for n>=64 and 1<=L<=n/16 also passes
independent proof review. The integer inequalities 17^17<2^70 and
65^128<2^15*64^128, with the exact n=64 base case, establish the all-order
comparison by induction. It is not inferred from a finite test alone.
`audit_alphabet_counting.py` separately verifies 744 stars-and-bars counts
by dynamic programming, exact entropy-root rational brackets, and both
row-count and fixed-mass pigeonhole inequalities for all n=64 through 512.

These are bounds on one unaugmented class signature. An authenticated
full inventory, endpoint margins, and reflected port constraints provide
additional information. The counting argument does not rule out a fixed
small alphabet forcing a full wiring table from that larger input, prove
unequal Hamiltonian counts, or identify a smallest collision order.

## Reproduction

From the published packet directory, using Python 3 standard library only:

    python -B wiring-audit/audit_ternary_quotients.py
    python -O -B wiring-audit/audit_ternary_quotients.py
    python -B wiring-audit/audit_alphabet_counting.py
    python -O -B wiring-audit/audit_alphabet_counting.py
    python -B wiring-audit/audit_binary_literal_switches.py
    python -O -B wiring-audit/audit_binary_literal_switches.py

Each command compares a frozen receipt by default. `--write` intentionally
regenerates that receipt. Normal and optimized outputs agree. The ternary
checker needs only `classes_input_n8.tsv`, already in the first published
audit; the information-bound checker has no external input files.
The literal-switch checker also needs the independently generated
intermediates from `FIRST_README.md` and the sibling
`binary-counterfeits/receipt.json`. An explicit `--witness PATH` overrides
that location. No producer Python module is imported by any checker.
