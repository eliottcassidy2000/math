# Independent review: colored relative fibres and marking efficiency

**Verdict: PASS** for `triples/COLORED_FIBERS.md`, relative to the explicitly
cited cyclic-triple fibre classification. Reviewed by the manuscript lane.
The elementary rigidity proof was reviewed separately in `rigidity_reviewer.md`.

## Exact labelled fibre grammar

Full port incidence vectors induce genuine vertex colors. For any flipped
interior edge, a port detects the flip in a cyclic triple iff its endpoints
have different port bits. Therefore fixing all cross-color orientations is
exactly the additional restriction imposed by the fixed-port augmented-C3
observer. Theorems9–11 make this the complete eliminated boundary fibre,
not just a necessary filter.

At a prime node, quotient reversal changes every cross-child edge. Preserving
all cross-color arcs forces every cross-child pair to have one color. Since
there are at least two nonempty children, this is equivalent to all leaves
of the node having the same color. The allowed choice count is2 or1.

At a linear node, an inversion changes all arcs between the two involved
child modules. It is permissible iff both children are pure and have the
same color. Mixed children and different-colored children cannot be crossed;
allowed permutations are exactly those inside consecutive pure-same-color
runs. Least-common-ancestor control of each edge forbids compensation between
levels. The inherited faithful independent node parametrization consequently
gives `2^p_mono product_R |R|!` without collisions.

This uses neither the retracted automorphism-order invariant nor ordinary
component-isotype invariance. It counts labelled orientations in a relative
observer fibre, not automorphisms, unlabelled classes, or new scalar counts.

## Added-port refinement criterion

The criterion is exact when checked at every old monochromatic prime node
and every old pure linear run with at least two children. Constancy preserves
all formerly allowed operations. Nonconstancy destroys at least one: a prime
reversal, or a permutation moving a pair of children distinguished by the
new bit. Constraints below a node cannot undo the changed cross-child edge.
The revised statement omits redundant singleton-run conditions, avoiding an
unneeded recursive necessity argument. Refinement can never introduce a new
allowed old orientation, so destroyed operations cannot be compensated by
new choices elsewhere.

## One engineered port

The recursive coloring proof is valid. Nontrivial children are mixed by
induction. Therefore the only pure children of any linear node are singletons;
alternating them along each consecutive singleton run removes every nontrivial
factorial. Prime nodes either contain a mixed child or are colored explicitly
with both bits. No prime factor2 survives. This gives a genuine incidence
vector for one physical new vertex and singleton relative fibre.

For n>=2, zero ports cannot suffice: T and its converse are distinct labelled
orientations but have identical all-subset H by path reversal. This gives the
lower bound independently of the cited fibre count. For n1 the no-port fibre
is already singleton. The construction makes no claim of primeness of the
extended tournament or preservation of H after adding its discriminator.

## Random ports and expected information

With q=2^k, independent fair port bits give independent uniform vertex colors.
Nontrivial fibre is exactly the union of monochromatic-prime events and
monochromatic unions of adjacent children at linear nodes. Each event involves
at least two leaves and has probability at most1/q. There are at most n-1
such events by the rooted-tree identity `sum(arity-1)=n-1`. No independence
between these events is used. This proves the uniform bound.

For a transitive interior, counting proper colorings of its ordered path
proves the exact `(1-1/q)^(n-1)` uniqueness probability. Thus the stated
fixed-success logarithmic scale and the unbounded additive margin for
worst-case success tending to1 are correct. At n1 the probability is1,
including the k0 convention.

For expected log2 fibre size, expand each maximal pure run's log factorial
as a sum over its intervals with weight `log2(ell/(ell-1))`. Each interval
has monochromatic probability q^(1-leafcount); each prime node contributes
one bit with probability q^(1-leafcount). Linearity proves the formula without
assuming event independence. The exact identity was additionally tested in
the free abelian group of prime valuations, avoiding floating-point logarithms.

## Reproducible controls

- `colored_fibers_checks.py`: all248,832 labelled three-color/tournament cases
  at order5, direct grouping by cyclic triples and cross-color edge bits.
  Also16,932 exact added-bit refinement cases through order4.
- `random_marking_checks.py`: all225 tournament/palette cases through order4
  for q=1,2,3;99 transitive cases and450 exact prime-valuation expectation
  identities. All arithmetic is rational.
- `engineered_port_checks.py`: independently implemented recursive color
  assignment on every1,099 labelled tournament through order5, then direct
  grouping confirms its relative fibre is singleton.

These reviewer scripts reuse the producer's literal strong-module tree builder
but use independent cyclic-triple detection, direct fibre grouping, recursive
coloring, and rational expectation calculations. They do not claim independent
reimplementation of the entire modular decomposition algorithm. Normal and
optimized runs agree; their JSON receipts give the exact universes.

Run each script by its repository-relative path and repeat with `python3 -O`.
The coordinator's final manifest should pin the final producer and reviewer
sources. No mathematical discrepancy was found in this extension.

## Final extension: sparse quadratic guards

**PASS, with a scope wording clarification requested.** A selected quadratic
port-start query is either vacuous or fixes the orientation of exactly its
queried pair, because its change is `(r_u-r_v) delta A_uv`. At the pair's
LCA this either locks a prime quotient's dualization or imposes a comparison
on a linear quotient. Consistency with the original linear order ensures an
acyclic comparison poset. Its total extensions count the allowed orders.
The inherited faithful modular-fibre grammar makes these choices independent
across nodes, proving the proposed product of linear-extension counts and
unlocked binary factors. This remains a conditional use of the cited fibre
classification, not a new audit of that classification.

The engineered binary port supplies a split pair across each required pair
of adjacent linear children and some split pair across prime children. A
mixed node cannot have all its cross-child pairs equal-colored: that would
make all its children one common color. Hence the selection always exists.
The tree identity gives at most n-1 queries. Locking adjacent comparisons
forces the complete linear order, while one prime comparison locks the prime
choice. In the transitive case every original adjacent pair is necessary,
since swapping just that pair preserves all other pair comparisons. Thus
this designed-query bound and its transitive worst-case sharpness are valid.

Scope clarification: retaining the complete internal C3 inventory still
preserves **interior H**. Sparse guards can change the **augmented or attached
system's H** by changing its unobserved boundary response. The coordinator was
asked to make that distinction explicit in the sentence about final H.

`sparse_guards_checks.py` compares the formula to direct labelled C3/pair-bit
fibres for all pair-constraint subsets through order4 and 51 deterministic
subsets at order5 (56,389 total tournament/constraint cases). It also directly
checks singleton fibres for the designed certificates on all1,099 labelled
tournaments through order5. Linear extensions are counted by literal child
permutations; C3 is tested by induced outdegrees. The producer module-tree
and engineered-color routines are reused. Normal and optimized runs agree;
`sparse_guards_checks.json` records the exact finite universes. No formula or
certificate discrepancy was found. These are finite controls of the proof,
not claims of arbitrary-size exhaustive computation or Chan–Pak coverage.
