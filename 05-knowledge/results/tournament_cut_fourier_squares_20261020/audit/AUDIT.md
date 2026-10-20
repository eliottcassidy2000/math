# Independent audit: Eulerian Fourier modes and commuting squares

2026-10-20. **PASS after a finite-check scope clarification. No mathematical correction is required.**

The final audited [proof note](../FOURIER_MODES_AND_COMMUTING_SQUARES.md) has SHA-256
`d37f5e05f5d62de6c4ffa5c9828399a588b64b584edaee711c09e0e6f2dfa3b9`.
The inherited [main cut-aggregate packet](../../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md)
and [order-five distance-four packet](../../tournament_graph_distance_four_20261020/DISTANCE_FOUR_SIX_CUBE.md)
remain frozen and unchanged. This audit covers the weighted cut aggregate, not a
black-edge, regular-return, ordinary-isomorphism, or converse-merged metagraph.

## 1. Result and correction lineage

The complete Fourier transform, all-signing forest expansion, triangle formula,
nullity lower bounds, peripheral spectrum, moment identities, and commuting-square
proof all pass. The original labelled graph's transitive leaves are a valid
counterexample to transferring bridge-freeness from the aggregate.

The handoff initially named hash `e31ccedf...`. Before the audit pinned a version,
the producer incorporated the parent's explicit `n>=3` qualifier for the peripheral
multiplicity equivalence, giving `a6ad6a90...`; the positive-apex cycle statement
already had its necessary `n>=2` restriction. The final audit requested one further
wording clarification: the producer verifier checks forest expansions for every
**Eulerian** signing through order five, not every arbitrary signing. Section 6
was narrowed accordingly, giving the final hash above. The theorem, verifier,
and receipts were unchanged by that final edit.

The independent audit additionally checks the formula for all 1,099 arbitrary
signings at orders one through five. This supports the more general theorem
without overstating the producer receipt's universe. No remaining required edit
was found.

## 2. Character transform and forest expansion

The annihilator of the cut space consists exactly of edge sets of even degree
at every vertex: pairing with a vertex cut returns that vertex degree modulo two.
Therefore Eulerian edge sets, including the empty one and isolated vertices,
index all quotient characters. The character count equals the quotient size.
Inserting each character into the incidence convolution proves the full spectrum,
with the factor `(-1)^|Z|` correctly converting path intersection to path-complement
intersection. This is not raw tournament H and does not impose Eulerian degrees
on a tournament.

For an arbitrary signing Z, expansion of the product of `1-2[e in P]` counts
all subsets F of Z contained in a Hamiltonian path. Such F must be acyclic with
maximum degree two. Its nontrivial path components have two orientations each;
its isolated vertices have only one. Contracting the nontrivial components leaves
exactly `n-r(F)` blocks, giving

    (-2)^r(F) 2^c_+(F) (n-r(F))!.

The inverse counting map is unique even when separate blocks happen to become
adjacent in the final word: the prescribed subset F still determines where the
blocks begin and end. The empty forest term, connected spanning-path case,
isolated vertices, and small orders one and two all satisfy the formula.

The positive-apex connection uses a graph of order `n+1` and unoriented cycles.
For `n>=2`, deleting the apex gives an undirected Hamiltonian path with exactly
two ordered orientations, so the factor two is correct. It is not a same-order
cycle identity, and no theorem about signed-cycle switching invariance is imported.
A direct hostile is already available at order three: the empty signing has
signed path sum 6, whereas switching its signs at one vertex gives signed path
sum -2. Thus the explicit warning against that transfer is substantive.

## 3. Triangle modes and full/component multiplicities

For a triangle signing, the only allowed subsets are the empty set, its three
single edges, and its three two-edge paths. The coefficients give

    S_n(C3) = n! - 12(n-1)! + 24(n-2)!
            = (n-2)!(n-4)(n-9).

The triangle has odd size, so its weighted eigenvalue is the negative of this
expression. Every displayed value and sign boundary in the note follows,
including the two zeros at orders four and nine.

At odd order the full edge set E is Eulerian, and complementing a character
multiplies its eigenvalue by the common generator parity sign. At order nine
this sign is positive. The 84 triangles and their 84 complements are pairwise
distinct Eulerian characters, giving 168 independent zero eigenfunctions on the
full quotient. On the even component, restriction identifies exactly the pairs
`Z` and `E xor Z`, because its annihilator in the full quotient dual is `{0,E}`.
The 84 resulting character classes are distinct, giving 84 independent zero
modes per component.

For precision, on the odd affine component the two restricted functions are
scalar negatives, rather than literally equal functions. After choosing an
origin, or simply passing to their common one-dimensional eigenspace, they have
the same identification and multiplicity. This does not change any stated bound.
The audit explicitly checks this sign distinction.

The order-nine figures are lower bounds, not an exact nullity classification.
No order-nine quotient spectrum or `2^28`-vertex graph is enumerated. At order
four, four zero modes are exact by the inherited cross-polytope spectrum. At
order five, the full spectrum has eigenvalue 24 with multiplicity 20, pairing
to multiplicity 10 on each component; division by edge weight 8 gives the
previously established support eigenvalue 3 with multiplicity 10.

## 4. Peripheral spectrum and normalization

An eigenvalue is a sum of `n!` signs. Its absolute value is `n!` precisely when
all signs agree, equivalently when its character annihilates every generator
difference. The inherited difference-span theorem gives a zero-dimensional
annihilator at even `n>=4` and the two characters `{0,E}` at odd `n>=3`.
This proves all peripheral multiplicities, not just the existence of the
constant character or a bipartition character.

The explicit `n>=3` qualifier on the equivalence is necessary. At order one the
matrix is `[1]`, not a two-component matrix despite `1=1 mod 4`. At order two it
is `[2]`. The note now handles these separately and correctly.

For `n>=3`, looplessness makes the spectral first moment zero. Squaring the
weighted convolution gives diagonal entry `sum_g C[0,g]^2`. Since the common
positive weight is w_n and the row sum is n!, this is `n! w_n`. Thus the second
moment is `|Q| n! w_n`; from order eight onward its character average is `2 n!`.
The support spectrum divides by w_n, and its squared moment divides by w_n².
No incidence/support or full/component normalization was lost.

## 5. Commuting squares, bridges, and the nonlifting hostile

For distinct nonzero generators g,h over F2, the four vertices

    q, q+g, q+g+h, q+h

are distinct. Their successive displacements g,h,g,h are support generators,
so they form a four-cycle. Additional diagonals do not invalidate it. The
inherited generator-span dimension is at least two for every `n>=4`, ensuring
that h exists for every chosen g. The argument works within either component
in the exceptional parity-disconnected cases.

At order three the support has only one generator, so it is K2 and its edge
is a bridge. Its positive entry is 6. Interpreting that weight as six parallel
incidence edges changes the individual-edge bridge question: deleting one of
six parallels leaves the adjacency, whereas deleting the support edge deletes
all six. The note keeps these conventions separate. At orders one and two,
the diagonal entries count identity-return incidences, not twice an undirected
loop degree; removing loops leaves no support edges.

A transitive tournament has one directed Hamiltonian word. Unique path witnesses
make its degree one in the original labelled graph, and for `n>=3` its move
changes at least one edge. Thus that edge is a genuine bridge. At `n>=4` its
aggregate image nevertheless lies on a commuting square. This is a precise
obstruction to lifting the quotient conclusion, not a claim about a differently
defined metagraph. No extra switching move is implicitly allowed upstairs.

## 6. Exact computational and provenance checks

Producer normal and optimized replays reproduce all **2,744** checks and the
existing receipt byte-for-byte:

- verifier SHA-256:
  `663222a67adbc27004a1e583faa73d69f49c9a08d3ee32332ba5c3c1af646292`;
- receipt SHA-256:
  `3a51b766683a4e521e15a6084c122ad5dd0e34a280164e2af24d0c6f2960104d`.

The [independent checker](independent_controls.py) imports no producer code and
uses no Walsh transform. It passes **4,661** live checks in normal and optimized
modes with identical [receipts](independent_controls.json). Its controls include:

- all Eulerian signings and complete direct signed-path spectra at orders 1–6;
- an independently enumerated cut kernel, compared character-by-character;
- all 1,099 arbitrary signings at orders 1–5, with all their forest expansions;
- each forest/block coefficient compared against a literal containing-word count;
- the positive-apex orientation factor for every signing at orders 2–5;
- first and second spectral moments and peripheral multiplicities;
- every generator's square certificate at orders 4–6 and the transitive-leaf
  hostile at orders 3–6;
- triangle values at orders 3–11 by an independent marked-position count:
  choose the three triangle-vertex positions, sign by consecutive marked
  positions, and multiply by `6(n-3)!`;
- all 168 order-nine masks and their 84 distinct restriction classes, without
  enumerating quotient vertices.

The full n=6 weighted spectrum independently reproduces

    (-240)^6, (-80)^72, (-48)^255, 16^540, 48^120, 144^30, 720^1.

This remains a bounded control, not an all-order spectral classification.

All six producer `SHA256SUMS` entries and all three inherited-source hashes
match. The cited signed-cycle theorem is used only to identify a distinct
nearby source, not as an unverified dependency of the path proof. No repository
or remote writes were performed. This is an elementary proof and exact-computation
audit, not a formal proof-assistant certification or priority claim.

## 7. Portable reproduction

From the published packet directory:

    python3 check_fourier_squares.py > /tmp/cut-fourier-replay.json
    cmp verification.json /tmp/cut-fourier-replay.json
    python3 -O check_fourier_squares.py > /tmp/cut-fourier-replay-optimized.json
    cmp verification.optimized.json /tmp/cut-fourier-replay-optimized.json
    python3 audit/independent_controls.py > /tmp/cut-fourier-independent.json
    cmp audit/independent_controls.json /tmp/cut-fourier-independent.json
    python3 -O audit/independent_controls.py > /tmp/cut-fourier-independent-optimized.json
    cmp audit/independent_controls.optimized.json /tmp/cut-fourier-independent-optimized.json

The accompanying publication manifest preserves the final producer bytes,
includes the independent evidence, and declares the two co-published inherited
packet links. It requires no scratch-directory paths when run after publication.
