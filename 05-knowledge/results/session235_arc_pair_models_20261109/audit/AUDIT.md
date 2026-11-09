# Session235 independent arc-pair audit

Status: PROVED (encoding implications, one-hole constructions and padding);
FINITE-EXACT (full oriented census through six and independent controls).
Producer census, theory PB1–PB7, core transport and portable replay: PASS. This report does not claim a complete oriented-core grammar.

## Inheritance and scope

Closest mechanism: Session234 DG3 ordered-sum multiplication, DG6 resolution
weights, and DG7 endpoint increments, in
`05-knowledge/results/session234_digraphs_20261109/PROOFS.md`.
Hostile example: the directed triangle with an unused added reverse arc.
Corrected near miss: no absent pair means semicomplete only under ordinary
adjacency semantics. Least-used sidecar here: the direction-specific number of
paths using a deleted arc. Concept board: raw pair decoder; tournament
completion; deletion response; SCC boundary; fixed-order versus all-order range.
Anchor: audit the interpretation. Niche: exact one-arc deletion. Wildcard:
directed cycles realize every integer >=3 outside the tournament category.
No LRC implication is asserted.

## 1. Raw bits must be decoded before counting

Retain the unequal-pair convention 10 = i->j and 01 = j->i, and set 11 = no
arc. The unresolved 00 needs an explicit choice:

1. If 00 also means no arc, the decoded category is oriented graphs. The
   decoder is not injective. An oriented graph with z missing pairs has
   exactly 2^z raw off-diagonal encodings. A no-raw-00 restriction selects
   exactly one encoding per oriented graph.
2. If 00 is forbidden, the allowed three states encode oriented graphs
   bijectively. The same no-raw-00 sector results.
3. If 00 means both arcs, all four states encode arbitrary digraphs
   bijectively. Restricting away 00 again gives oriented graphs.

Consequently, under each interpretation, absence of off-diagonal raw 00 does
not force a tournament: 11 can still be a missing pair. A tournament is the
unequal-only sector. Diagonal bits are irrelevant to vertex-distinct paths
when the decoder assigns them only self-loops.

If one instead multiplies raw entries along permutations, one has silently
returned to ordinary adjacency semantics, where 11 is reciprocal. The path
functional must be applied to the decoded adjacency matrix.

## 2. One missing pair already restores both omitted counts

All matrices in this section are ordinary adjacency matrices, NOT raw codes.
The accompanying JSON records the alternate raw codes separately.

H=7, order five:

    00001
    10000
    01000
    11100
    01110

The sole missing pair is {0,2}. The seven paths are
04321,10432,21043,30421,31042,32104,43210.
Both tournament completions have H=9. Deleting either chosen arc removes
exactly two paths. To encode this with 11=empty, change the entries (0,2)
and (2,0) to one. Every other off-diagonal pair is unequal.

H=21, order six:

    000011
    100000
    110000
    111000
    011100
    001110

The sole missing pair is {1,5}. Its two tournament completions have counts
41 and 23. The 21 actual paths are listed explicitly in `independent.json`.
The alternative raw matrix changes entries (1,5) and (5,1) to one. Again
there is no off-diagonal 00 under any of the three interpretations above.

For a tournament T and existing arc e,
H(T-e)=H(T)-number of Hamiltonian paths in T containing e.
This is a literal partition of the path set. The finite omission theorem
says nothing that prevents the difference from being 7 or 21. Likewise,
if D is missing exactly {u,v}, its two completion counts are H(D)+N_uv and
H(D)+N_vu, where N counts newly enabled paths. Their odd parity implies
N_uv and N_vu are even when H(D) is odd; it does not prohibit that H(D).

Ordered sums T_a -> D -> T_b preserve H(D). All added cross pairs have
exactly one arc. Thus these witnesses generate one-missing-pair examples
at EVERY order n>=5 for 7 and n>=6 for 21. This is an existence family,
not an exhaustive classification.

Every oriented graph embeds in a tournament on the same vertices, and
adding arcs cannot decrease H. The maximum tournament H on four vertices
is 5, and on five vertices is 15 (independent exhaustive checks here).
Therefore those thresholds are minimal even among all oriented graphs.
By the inherited all-order tournament omission theorem, zero missing pairs
cannot realize 7 or 21 at any order. Exactly one missing pair is therefore
the sharp minimum at every attainable order.

## 3. Independent numerical coverage

`independent.py` does not import producer code or its endpoint recurrence.
It scatters each vertex permutation into all compatible tournament or
oriented pair codes. For tournaments it separately accumulates which pair
arcs each path uses, then subtracts those counts. Witnesses are subsequently
checked by a separate literal enumeration of all permutations.

Run from repository root:

    python3 05-knowledge/results/session235_arc_pair_models_20261109/audit/independent.py

Exact one-missing-pair counts, labelled and loopless:

| n | H=7 | H=21 |
|---|---:|---:|
|2,3,4|0|0|
|5|600|0|
|6|10080|8640|

The number of tournament/arc deletion presentations is twice each number:
each missing-pair graph has precisely two tournament completions. No
unproved isomorphism quotient or integer division by group orders is used.

All oriented graphs are independently scanned through order five, universe
3^(n choose 2). At n=5 there are 1440 labelled H=7 graphs, versus 600 with
exactly one missing pair; H=21 is absent. Maximum H for n=1,...,5 is
1,1,3,5,15. A subsequent independent direct-scatter comparison checks every oriented
entry at order six, as documented below.

## 4. All-order boundary

The no-missing-pair theorem is genuinely all-order. The sharp minimal order
and minimal defect statements above combine finite maxima, explicit witnesses,
and padding. They do not require an all-order oriented-graph classification.

The directed cycle on h vertices has exactly h Hamiltonian paths for every
h>=3: its path is forced after choosing the starting vertex. Together with
an isolated two-vertex graph (H=0), a directed path (H=1), and the oriented
four-vertex graph with rows 0100/0010/1001/0100 (H=2), this shows the oriented
graph spectrum over all orders is all nonnegative integers. This elementary
range fact does not solve the tournament spectrum conjecture.

Do not copy Session234's finite-core proof into the oriented category:
strong oriented graphs need not have a Hamiltonian cycle. For example the
arcs i->i+1 (0<=i<n-1), (n-1)->1 and (n-2)->0 form a strong oriented graph
for every n>=4, with exactly two Hamiltonian paths: the increasing path,
and (n-1),1,2,...,(n-2),0. Indeed, if vertex 0 is first, the increasing path
is forced; otherwise its only predecessor is n-2, and because its only
successor 1 then precedes n-2 in any completed path, it must be last. The
remaining path is forced to start at n-1. Thus fixed H does not bound the
order of strong oriented cores, already at H=2. No direct extension of
DG4 or its strong-floor bound is valid here.

## 5. Source attribution and exact implication

Primary source checked on 2026-11-09:
[Chan and Pak, Computational complexity of counting coincidences,
arXiv:2308.10214v2, Section 6.5, Conjecture 6.2](https://arxiv.org/html/2308.10214v2#S6.SS5).
It conjectures that the tournament path-count range is precisely the positive
odd integers excluding 7 and 21, and expressly credits bof and Gordon Royle
(March 2016). It does not claim an efficient construction or prescribe a
vertex bound. The conjectural part is attainment of all remaining odd values;
changing the underlying graph category does not settle that part.

The [original MathOverflow discussion](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament)
contains bof's strong-component argument for 7, Royle's finite data, and bof's
March 11, 2016 comments reducing the 21 exclusion to the strong-order floor.
Do not claim priority for the exclusions or call them new Chan-Pak results.
The in-repo all-order proof is
`01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`.
Its historical source-audit link to
`05-knowledge/results/rotation_owners_spectrum_20261020e/spectrum.md` is absent
in this checkout; the primary links above were inspected directly instead.


## 6. Final producer audit: PASS

Reviewed `census/census.cpp`, `census/verify.py`, `census/REPORT.md`, and
`theory/PROOFS.md` PB1–PB7. The pair indexing, zeta transform, reciprocal
coordinate identity, and full permutation-orbit quotient are correct. The
integer range is safe at n<=6 (at most 720 paths). Replaying producer code
in `/tmp` reproduced all 24 persistent CSV files byte-for-byte.

The independent `scatter_compare.cpp` explicitly enumerates all compatible
oriented supergraphs of EACH permutation, with no zeta transform and no subset
DP. Its output equals EVERY oriented table entry at every order through six,
including all 14,348,907 order-six values. In that order it scatters
42,515,280 path/supergraph incidences and confirms H7=255600 and H21=29520.
`scatter_compare.txt` and `producer_replay.txt` retain the results. This upgrades
the oriented sweep from a sampled second-algorithm check to a full independent
algorithmic check. It does not independently rescan all semicomplete entries;
the inherited marginal verification and endpoint samples retain their stated
scope.

To reproduce the full cross-check, compile producer `census.cpp` into a temporary
directory and run it there for orders 1..6, then compile this audit's
`scatter_compare.cpp` and give that temporary directory as its sole argument.
Binary arrays are scratch only; their recreation is necessary after cleanup.

The all-order identities PB2/PB3 are literal path-set double counts and
inclusion-exclusion. PB4's spanning-linear-forest expansion is valid, including
the c! factor because components are unordered and isolated vertices count.
Only c=1 survives modulo 2 for n>=1. Reverse-transpose contributes no sign.
The second determinant proof also respects the square-zero coefficient ring.

The cut-vertex lobe identity PB6 is valid: both nonempty lobes must appear as
single blocks on opposite sides of v. For its cycle+hub lobe, paths ending at v
are indexed precisely by I and paths starting at v precisely by O. The pendant
cycle contributes one path of each endpoint type. Therefore the total is h,
with n=h+m arbitrarily large; H2 is supplied separately by two cycle lobes.
This decisively refutes a finite alphabet of actual strong oriented cores for
7 and 21 themselves, while leaving other compressed grammar notions open.

`theory_controls.py` independently checks 2460 cycle/hub/pendant cases, all
nontrivial partitions for h=3..8 and pendant lengths 3..7; 36 two-cycle
bouquets; and all 4165 digraphs of orders 1..4 for complement-reverse parity.
All pass. These are controls for the proofs, not their replacement.

The producer census's conditional statements about same-raw decoding are finite
unless explicitly transported by the inherited complete semicomplete SCC
classification. That transport is valid even when deleting reciprocal arcs
splits an original strong core: ordered-sum multiplicativity only requires
all cross arcs to point forward, not that each prescribed block remain strong.
Accordingly, checking the finite complete source core list suffices for an
all-order image statement. Here the map is the relabel-equivariant involution τ, so source and decoded
automorphism groups are identical. This justifies using the same orbit weights
on both sides; no arbitrary division by a presumed group order occurs.


## 7. Core transport and generating functions: PASS

Independently recomputed all Hamiltonian path counts and full automorphism
groups of the 37 inherited source cores (two H3, three H7, thirty-two H21)
from the frozen matrices in `census/core_transport.json`, without importing
producer counting code. Checked each decoded matrix equals the explicit τ
transform. The group counts and all target image values match.

Combined with the inherited ALL-ORDER complete source grammar, this proves:

- Every semicomplete H7 source decodes to H1.
- Semicomplete H21 sources decode precisely to {1,3,5,7,15} over all orders.
- All five H21 image values already occur at order six and every larger order
  by singleton padding. At order five the image is {1,3,5,7}.

Mixed core words have decoded product 1*1 or 1*3; this creates no extra value.
The new exact joint generating functions retain source isomorphism classes;
τ being an equivariant bijection also identifies them with decoded image
classes. `core_transport_check.py` independently enumerates singleton gap
placements for source words and verifies every displayed GF coefficient at
orders 1..12. This is a finite control of the exact word argument, not a
finite-to-infinite extrapolation. All 24 checks pass.

The pair-model verifier was rerun independently and reproduces its validation
JSON byte-for-byte. All requested producer components now pass this audit.


## 8. Portable full replay

`python3 audit/replay_all.py` can be run from any working directory using the
appropriate path to the script. It uses a temporary sibling layout, copying
only source programs and the inherited Session234 generator/marginal inputs.
All census binaries and intermediate outputs stay outside the source packet.
It regenerates and verifies all persistent census CSVs, the full independent
oriented scatter, the models and theory checks, conditional histograms, core
transport, and every independent audit control. Successful replay compares
36 persistent artifacts and verifies the reference files were not changed.

Optional `--receipt PATH` writes only the requested JSON receipt; inside the
packet it permits only `audit/replay_all.json`. Without the option it writes
no packet files. The frozen receipt is `replay_all.json`. It contains stable
artifact hashes and scope statements, but no times, machine names, temporary
paths, or other nondeterministic fields. Compiler and scratch products are
removed on exit. Requirements are Python 3 and a C++17 compiler (`CXX`, default
`g++`). The inherited source classification is an input, not independently
re-proved by the replay runner.


## 9. Recovery verification against current main

The earlier finite checkpoint and completed final audit files were recovered
before publication. A fresh full replay reproduced the frozen deterministic
receipt exactly: all 36 artifact comparisons passed. An additional read-only
proof review checked PB1-PB7 and the complete joint generating functions,
and independently rechecked the 37 cores and 24 coefficient controls through
order twelve. No mathematical correction was needed. The provenance sentence
in PB4 now explicitly credits Berge's complement-parity theorem and links the
concurrent THM-5071 weighted-cycle refinement. No novelty or tournament
completeness claim is added.
