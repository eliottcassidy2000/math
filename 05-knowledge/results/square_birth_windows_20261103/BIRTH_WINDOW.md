# Complete finite-horizon birth controllers

Status: **PROVED general reduction; FINITE-EXACT minima at the stated sizes.**
This extends the inherited [native-26 proof](../eleven_wedge_honeycomb_20261102/NATIVE_26.md)
and uses the inherited fixed-order 33/71/637 controllers from the
[three-window packet](../twenty_curvature_ports_20261102/THREE_WINDOW_CARRIER.md).

## Object and permitted operations

Q_n has vertices 1,...,n and undirected edges {i,j} with i+j square, i!=j.
A guard source is an orientation x of ALL these edges. Edge bits are ordered
by larger endpoint, then smaller endpoint, so new edges append. Bit 1 means
smaller -> larger. A path name is its literal directed vertex sequence;
names are lexicographically sorted independently at each stage.

The immediate output is the set L_n(x) of legal named Hamilton paths.
For P in L_n(x), command K_P retains P's edges and reverses every other edge:

    K_P(x) = x xor (all_edge_bits xor path_edge_mask(P)).

This is the restriction of the actual tournament path complement, not global
converse. It is a partial involution, and the selected P stays legal.
A birth appends n+1 and supplies the orientation bit of EVERY new guard edge.
The source orientation is retained; there is no class-valued reset or sampling.
Arbitrarily many legal commands may occur before, between and after births.
The stage clock and new-edge names are retained. Global converse is allowed
as well and is verified to descend without refining the reported quotients.

An observer is minimal here if no two of its states can be identified while
preserving all future outputs for this exact operation language. Vertex-name
quotients, unrestricted re-marking and averaged joint-source laws are separate
problems.

## A finite reduction for every finite horizon

Let A_n be the active sources L_n(x)!=empty. Fix a terminal size N, and define

    S_N = A_N,
    S_n = A_n union projection_n(S_(n+1)).

**Proposition (all finite nested guards).** Explicitly retaining S_n and one
residual state of weight 2^|E_n|-|S_n|, when nonempty, is sufficient for every
word through horizon N. Backward refinement with the next-stage quotient
gives the coarsest observation congruence at each stage.

Proof. A legal K_P keeps a path and hence stays within A_n. A source outside
S_n has no legal old command and EVERY new-edge assignment lies outside
S_(n+1). The whole excluded cube therefore has the same present and future
observations. This is a residual block, not a selected representative whose
arcs are silently re-used. At N, start with legal-path outputs and refine
by the labelled successor blocks. At n<N, initialize with legal-path
outputs and the vector of next-stage classes for ALL new-edge assignments,
then refine by old commands. Stable classes preserve every permitted word
by induction on the word and on the remaining births. Every split has a
distinguishing observation: either an immediate output, a birth followed by
a distinguishing future word, or one old command prepended to such a word.
Thus no split can be undone by an exact observer. This proves minimality.

When converse permutes each path bank, the S_n are converse-invariant.
Converse conjugates K_P to K_(reverse P) and complements the new-bit word.
The output language is closed under these changes, so converse descends.
The program additionally checks this on every explicit source and birth.

The proposition is a finite-horizon compiler, not a bound on how many states
an infinite horizon needs. It generalizes the inherited one-birth residual
reduction. In arbitrary guards an empty residual should simply be omitted;
the present three graphs all have nonempty residuals.

Unrolling gives S_n = union_(j=n)^N projection_n(A_j). Before its first active
stage a source has no legal path-complement command, so this direct-birth
support also describes relevance when arbitrary intermediate commands are
allowed. Each projected path event is a cube on the old edges: delete the
future vertices from the named path and retain its ordered spanning path
forest. If that forest has r components, it fixes n-r old edge bits. In a
connected old guard with cycle rank beta, its cube has dimension
|E_n|-n+r=beta+r-1. This is an elementary way to construct support directly
from future path banks without enumerating their free future edge bits.
Overlap and subsequent command refinement are still required.

## Exact universes and minima

Q25,Q26,Q27 have respectively 32,34,36 edges and 20,24,70 directed paths.
Leaf 18 is their only degree-one vertex, so every Hamilton path has 18 as an
endpoint. Exhaustive DFS from 18 followed by reversal finds the complete
bank. Each path fixes n-1 edge bits. The union of its completion cubes is
exactly A_n; it has sizes 4,864, 11,232 and 62,952.

| Stage and horizon | Explicit S_n | Residual weight | Active-source classes | All classes |
|---|---:|---:|---:|---:|
| 26 through 27 | 68,040 | 17,179,801,144 | 214 | **761** |
| 25 through 27 | 68,760 | 4,294,898,536 | 84 | **903** |

The complete refinement sequences, including the repeated stable count, are
673,745,761,761 and 869,901,903,903. Each starts from immediate legality plus
the full successor-class vector. A separate inverse-image block-splitting
algorithm gives exactly the same partition, up to class renaming.

Every source in a block is checked for equal outputs, command targets,
birth targets and converse targets. Every class weight is its exact number
of actual orientations. Legal complements preserve class weights, and
summing over ALL source classes and ALL new-bit assignments recovers every
target class weight exactly. This includes the large symbolic residual.
The complete tables, paths, edge orders and witnesses are in `birth.json`.

At 27, precisely 6,144 active orientations restrict to A26 and 56,808 do not.
The latter can arise from spanning two-path forests at 26. Deleting an
internal 27 from a new Hamilton path does not leave an old Hamilton path.
The two neighbors 9 and 22 do not form a square-sum edge. This explains why
the old active cube is inadequate without an appeal to cardinalities alone.

Deleting BOTH 26 and 27 from the complete 70-path bank at 27 leaves three
old25 path components in 62 cases and two in 8 cases. Their projected cubes
have respectively 10 and 9 free old-edge bits, compared with 8 for a single
old25 Hamilton path. This gives literal multi-run attachment data behind the
growth from 33 to 903 states; scalar path totals cannot supply it. The incoming
[Session148 insertion recurrence](../session148_square_sum_20261103/oriented/REPORT.md)
counts the one-birth two-path contribution. The present controller additionally
closes the named joint observations under repeated actual complements and
two successive births, on every source including dormant ones.

## Literal one-bit hostile

In the birth edge order, take old 26-source masks

    x=5392403568,    y=5392403576=x xor 8.

Only edge {2,7} changes. Both lie in fixed-26 class 6, so no old legal
complement word with named-path observations separates them. Supply the
same new-bit word 2 on edges (9,27),(22,27), namely the ear 22 -> 27 -> 9.
Their immediate 27-path names are respectively {26,27,53} and {53}.
Literal directed-arc membership independently verifies these masks.
The full sequences for each path name are stored in the fixed-27 path bank.

This is an actual omitted source bit exposed by a literal birth. It is not
an inference from equal Hamilton counts or from an arbitrary abstract probe.

## Relation to labelled tournaments

Every guard orientation at 25 has 2^268 labelled tournament completions.
Births 26 and 27 introduce respectively 23 and 24 additional nonguard arcs;
the final free factor is 2^315 because binom(27,2)-36=315. All these arcs
remain outside the square-sum guard forever along this prefix family.
Their actual reversals do not affect ANY consumer in this specified language.
This is an exact sufficient quotient of the tournament operations for this
target, not a reconstruction of the full tournament or its automorphism group.
In particular it supplies no new A000568 coefficient and does not replace the
81/27 transporter control for unrestricted tournament composition.
