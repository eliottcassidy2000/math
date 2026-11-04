# Twenty actual reader repairs retain the original names

**FINITE-EXACT:** both missing24/44 returns at each clock t=1,…,10 have
explicit native programs from the specified original source graphs to the
actual reader witnesses. **PROVED:** a two-helper bridge fallback makes
each required old-arc repair constructive without an unbounded Hamilton-word
search. **OPEN:** all-clock target acceptance, a reusable common return
word, and regeneration of two output branches.

## Actual source, actual target

The sources are the previous sibling families on prefix names10,…,10+t−1,
suffix names0,…,4, and root99:

    family24 source: H=−13+21·2^t,
    family44 source: H=−11+21·2^t.

The new target has 3H+4 paths in each case. The records in
[`reader/witnesses.json`](../../reader/witnesses.json) provide its complete
named matrix, selected seed and reader word. Every old physical name occurs
injectively, and root99 remains the outer root. Some old prefix names move
into the target component; their old incidences can change. These are actual
source-preserving-name rewrites, not replacements of a source by another
graph with equal H.

The closest mechanism is the inherited one-helper compiler in
[`METAGRAPH.md`, §5](../../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md),
implemented in `metagraph_probe.py`; its finite seed-switching application is
[`recursive_seed_switching_20261103/transport/NOTE.md`](../../../recursive_seed_switching_20261103/transport/NOTE.md).
The canonical hostile is the transitive triangle's long arc, which need not
appear consecutively in any original Hamilton word. The corrected near miss
is to regard a count-formula hit as a physical native operation. The retained
sidecar here is the exact actual birth cut and two guarded words for each
arc, plus the complete named endpoint matrices. The live concepts are reader
acceptance / source names / changing prefix roles / guarded complements /
joint endpoint and phase readouts.

## A guaranteed bridge fallback

For an actual arc a→b in any tournament T, construct a Hamilton path P of
T−{a,b} by the elementary insertion algorithm. Add a new bridge z such that
b→z and z beats every other old vertex. Then

    Q=(a,b,z,P)

is a directed Hamilton word of T+z, including when P is empty. In particular
b is not last. Add the inherited append-helper h, which beats only a. The
two native guards are

    Q1=(a,b,z,P,h),
    Q2=(a,h,b,z,P).

After the first complement, a→h and h→b are directed; all other consecutive
arcs of Q2 were consecutive in Q1 and were retained. The two path-edge sets
differ on the original old pairs only at {a,b}. Delete h and z. Exactly the
specified old pair reverses; every other old arc and every old name is
unchanged. This needs two helper births, two native complements, two helper
deletions, and peak order |T|+2. The insertion-built guard takes polynomial
time in the graph order. It is an alternative implementation budget, not an
improvement on the inherited one-helper peak-order theorem.

The actual compiler first tries the inherited one-helper method with bounded
forward/inverse path searches (1,000 states per attempt). A failed bounded
search triggers the bridge formula. Thus it never treats search failure as
nonexistence or waits for an unbounded Hamilton-language enumeration.

Neither helper/native intermediate graphs nor completed individual arc
repairs are required to preserve H, strongness, root stars, phase, or common
path frames. Only the declared source and final target are asserted strong.
This is not fixed-order native or black-core connectivity.

## Complete literal programs

For each target, first birth its fresh permanent vertices with their final
incidences to all present names. The remaining differences are old-old arcs.
Reverse those pairs one at a time; the actual graph after each two-native
gadget is checked against an independent literal single-arc reversal. An
elementary pairwise native interpreter validates both guards in the actual
intermediate matrices. Final equality is checked after transporting matrix
indices by physical names, not by canonical class or scalar count.

| Clock | Family24 source→target | Arc gadgets | Family44 source→target | Arc gadgets |
|---:|---|---:|---|---:|
| 1 | 29→91 | 9 | 31→97 | 4 |
| 2 | 71→217 | 11 | 73→223 | 5 |
| 3 | 155→469 | 12 | 157→475 | 6 |
| 4 | 323→973 | 19 | 325→979 | 15 |
| 5 | 659→1981 | 25 | 661→1987 | 25 |
| 6 | 1331→3997 | 23 | 1333→4003 | 26 |
| 7 | 2675→8029 | 31 | 2677→8035 | 30 |
| 8 | 5363→16093 | 26 | 5365→16099 | 35 |
| 9 | 10739→32221 | 28 | 10741→32227 | 47 |
| 10 | 21491→64477 | 41 | 21493→64483 | 52 |

There are470 arc gadgets /940 native complements in total. Of these,
449 use one helper and21 use the guaranteed two-helper bridge. The20
programs use13 permanent births and491 temporary births/deletions; no old
name is deleted. Maximum peak order is19. These costs describe the specified
programs; no optimality claim is made.

## Independent counts and freshly decoded sidecars

The inherited compiled gap evaluator computes H directly from the source
and endpoint matrices, independently of the reader's arithmetic recurrence.
Its required feedback deletion is found by an exact triangle-hitting search
and checked by the engine to leave a transitive tournament. Both source and
target counts are checked against the stated 3H+4 relation.

At every endpoint, ports99,0,1,2,3,4 have fresh E_r(2), golden triples,
F49-plus-charge coordinates and F343-vector coordinates. Six ordered joint
J queries are freshly computed in both directions on pairs{99,0},{99,2},
{0,4}. All names and full rows remain available for higher joint responses;
these selected queries are not asserted to be a complete small interface.
The reader sidecar supplies additional promoted-prefix-port observations.

The chart formulas are the inherited ones in
[`PHASE_CONTROLS.md`](../../../multiport_defect_grammar_20261103/transport/PHASE_CONTROLS.md).
Every changed graph is decoded anew. No physical Frobenius action, phase
invariance, or transfer of the unrelated81/27 source symmetry labels is
asserted. The actual81/27 controls were unchanged and are inherited, not
rerun or identified with these new targets.

Reproduction from repository root:

    python 05-knowledge/results/c4_port_returns_20261103/native/reader_repairs/check_repairs.py
    python -O 05-knowledge/results/c4_port_returns_20261103/native/reader_repairs/check_repairs.py --output 05-knowledge/results/c4_port_returns_20261103/native/reader_repairs/receipt_optimized.json

The normal and optimized receipts must be byte-identical. In addition to the
20 programs, every oriented arc of every labeled tournament at orders2–4
is independently replayed through the bridge fallback:410 controls. The
compiler and counter binaries live only in temporary directories.

## What this does and does not overcome

The fixed homogeneous-clone no-go and these positive finite repairs are
compatible. The targets change the exterior component/reader and can promote
old prefix vertices into different incidence roles. They therefore leave
the fixed marked-quotient hypothesis in the exact way the obstruction
requires. Twenty accepted and physically realized targets establish those
twenty cases. They do not prove that the bank accepts every later target,
that infinitely many changes are necessary, or that the output admits both
next returns.
