# Shared-pair flags exactly characterize local partner gluing

23 October2026. **PROVED**, elementary all-size full-record compatibility
lemma. **FINITE-EXACT** exhaustive controls and minimal examples are in
`partner_gluing.py` and `partner_gluing_receipt.json`. No claim of literature
novelty is made.

## Statement

Let C_i be vertex subsets of one finite actual owner set V. On each C_i
supply a tournament S_i and a Hamiltonian path P_i, with literal owner
identifications retained. Assume the source tournaments are coherent:
whenever two cells contain the same pair{u,v}, their source orientations
of that pair agree. Equivalently they are induced restrictions of at least
one source host tournament; uncovered source pairs may be oriented freely.

For a pair in C_i, put

    c_i(u,v)=1 if u and v are consecutive in P_i, and0 otherwise.

Define the local partner K_i by preserving the source arc on consecutive
path pairs and reversing it on every other pair.

**Gluing theorem.** There exists a tournament K on V with
K[C_i]=K_i for every i if and only if

    c_i(u,v)=c_j(u,v)

for every pair{u,v} contained in both C_i and C_j. If this holds and e pairs
of V are uncovered by all cells, exactly2^e labelled partner hosts extend
the specified local records.

**Proof.** Fix a covered owner pair. Its source orientation is the same in
every cell containing it. The local partner reverses that common orientation
exactly when c_i=0. Hence the prescribed partner orientations agree exactly
when the flags agree. Necessity follows by restricting a putative host K.
For sufficiency assign every covered pair its common partner orientation.
Each uncovered pair may independently receive either orientation. A
tournament requires one orientation per distinct pair and no additional
consistency condition on triples, so these assignments are precisely the
2^e extensions. QED.

If source coherence is absent, flag equality alone is not the criterion.
The general edgewise criterion is equality of the resulting oriented
partner arcs; opposite source orientations and opposite flags can even
cancel. Source coherence is an explicit hypothesis, not a conclusion
from matching owner names.

## Small exact obstructions

Take three ordered owners a<c<b. Cell{a,b} has source/path a->b, so its
partner preserves a->b. Cell{a,c,b} has the transitive source and its unique
path a->c->b, so its partner reverses a->b. The sources agree on their
shared pair and the owner identifications are literal, but no partner
host exists. Three owners are minimal if unequal cell sizes are allowed:
on at most two owners every shared pair lies in a two-vertex cell and is
consecutive in its unique Hamiltonian path.

For equal three-vertex cell sizes, use the coherent source host ordered

    a<d<b<c.

The cells{a,b,c} and{a,d,b} have paths a->b->c and a->d->b. Their common pair
{a,b} is consecutive only in the first path, so the partner arcs disagree.
Four distinct owners are minimal for two distinct three-vertex cell sets.
This is a minimal distinct-cell example, not a classification of multiple
different observers on one identical cell.

## What a global partner host does not supply

Gluing does not supply a global Hamiltonian observer whose consecutive-pair
flags induce all the local flags. For example, take a transitive source on
four owners a<b<c<d and the three two-vertex cells{a,b},{a,c},{a,d}. Their
path flags all equal1. There is no repeated shared pair, so the partner
constraints glue. But any global observer meeting all three flags would
have three consecutive neighbours at a, impossible in a linear path.
Every local path is legal; the missing global coordinate is the path's
degree-two and connectedness structure. The present theorem does not
claim a global observer or a common-host realization of arbitrary local
records whose source arcs were never authenticated.

## Exact compatibility factor and converse

For binary endpoint records at each cell, pairwise compatibility is decided
by comparing the full shared-pair flag vectors, after checking source
coherence and owner embeddings. This supplies a literal pairwise predicate
to `../compatibility/NEUTRAL_FIBRE_LIFTING.md`. If the all-A and all-B old
assemblies are compatible, forbidden mixed choices compile to that theorem's
implication poset. If either endpoint fails, its implication-only interface
does not apply without retaining the additional restrictions.

Converse transport replaces every S_i by its converse and every P_i by its
reversal. Consecutiveness is unchanged, all partner orientations reverse,
and the same flag-agreement condition holds. Literal clock positions are
reflected cell by cell. Matching a source owner's coordinates without its
shared-pair flags loses precisely the information detected by the examples.

Source: actual coherent induced source cells with marked Hamiltonian paths.
Target: compatible local partner restrictions of a single tournament host.
Preserved predicate: exact existence of a partner host, including the
2^e extension count. Forgotten by that host alone: which observers and
clock coordinates generated its local arcs. Required sidecar: the owner
embeddings, local full records and path-consecutiveness flags. Neither
electrical guards nor normalized-resource capacity follow from gluing;
they remain separate factors of a native packet compiler.

## Reproduction

The checker exhausts all64 labelled source tournaments on four owners,
all pairs of distinct three-owner cells, and all their legal local
Hamiltonian observers. It compares the flag criterion with direct enumeration
of all64 candidate partner hosts, checks extension counts and converse, and
includes both incompatible examples and the compatible/global-observer
hostile above.

    python 05-knowledge/results/compatible_charged_helpers_20261023/controller/partner_gluing.py
    python -O 05-knowledge/results/compatible_charged_helpers_20261023/controller/partner_gluing.py
