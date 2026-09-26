# 233: shared scales, height certificates, and Hamiltonian path coordinates

**Status: PROVED scoped results / INDEPENDENTLY AUDITED / FINITE-EXACT
controls / OPEN Collatz.** Session 2026-09-26, starting at802d4893e.
The user's corrected seed is233. The earlier223 results retain their
stated scope; this session follows different recovered mechanisms.

The strongest advances are a rigorous separation between lower and ordinary
pairing density, an elementary pointwise descent-or-rank alternative, and an
exact source-height filtration of a marked path's comparison graph.

## Inheritance and concept board

Anchor: globally consistent descent certificates and their integer height.
Niche: the depth233 carry collision and incident gcd magnitudes.
Wildcard: the user's smaller-tournament/Hamiltonian-path construction.

Closest proved mechanisms: [THM-4491's tree costs](../../01-canon/theorems/THM-4491-pairing-two-step-tree-extra-density.md),
[THM-4490's prime banks](../../01-canon/theorems/THM-4490-collatz-affine-word-sunit-specialization.md),
and [THM-354's component cuts](../../01-canon/theorems/THM-354-good-cut-scc-count.md).
Canonical hostile: equal time, equal odd count, and every prefix slope
greater than one still permit the depth233 source merger.
Corrected near miss: lower-density attainment does not imply ordinary-density
attainment. Least-used sidecars: common assignments across intermediate
cutoffs, the magnitudes of incident carries, and labelled height thresholds.

| Concept | Operation | What changed |
|---|---|---|
| Shared pair bits | Optimize several cutoffs under one assignment | Positive-scale incompatibility forces density oscillation |
| Lower versus upper density | Retain both liminf and limsup | The old lower-density optimum is not an ordinary density |
| Prime valuations | Apply triangle inequality at a maximal relation coefficient | Explicit gcd dominance and source-height certificates |
| Depth233 carries | Compose blocks and decode the last differing digit | This collision gadget retains capacity two, rather than multiplying fibres |
| Hamiltonian spine | Transfer reduced edges to labelled chords | Exact edge coordinates survive; incidence and cycles do not |
| Affine source height | Keep the intercept and vary the source within one cylinder | Components split at exact rational contact thresholds |

## What233 actually represents in the recovered work

| Occurrence | Exact role | Transfer or boundary |
|---|---|---|
| [Old flow section8](crossroads_20260926_flow.md) | Merger time233 with153 odd steps, between N and N-4 | Ordered carries permit collisions despite equal slope; no minimality claim |
| [Mirror and clock work](collatz_procgen_20260922_q1_mirror.md) | Upper approximation233/147 to log_2(3), with2^233>3^147 | A slow-descent clock, distinct from the153-odd-step merger |
| [Procedural order laws](collatz_procgen_20260922_order_laws.md) | The segment231 to233 has17 odd steps and27 halvings | Its slope contracts but its positive carry raises the endpoint |
| [THM-316 staircase](../../01-canon/theorems/THM-316-staircase-antipalindrome.md) |233 Hamiltonian paths on a particular eight-vertex staircase | Independently reproduced; neighboring counts1,5,29,233,2489 are not Fibonacci |
| Fibonacci addresses | F_13=233 and144*377-233^2=-1 | Cassini is an identity; no descent-preserving map to the merger was found |

The theorem number, a path count, a source value, and a time index are
different types. The successful transfers below use explicit maps instead
of identifying those types because their displayed integers agree.

## 1. Common cutoffs force density oscillation

**[THM-4492](../../01-canon/theorems/THM-4492-pairing-two-cutoff-density-separation.md):
PROVED, independently audited.**

For global two-step pairings let A_F(X) count flipped pairs through X.
THM-4491's attained minimum lower density alpha satisfies

    alpha <=0.2953650955221956... .

The previously open ordinary-density attainment target is now REFUTED.
Every global two-step member has upper density

    u >=821510388809/2677850419968
       =0.30677978974599224... .

This is also a lower bound on any existing ordinary density. A member
attaining lower density alpha must instead satisfy the stronger bounds

    u >=0.3116986076593383...,
    u-alpha >=0.016333512137142698... .

The first finite obstruction is tiny: C(2)=0 and C(5)=1, but one common
assignment has minimum A_F(2)+A_F(5)=2. The two optima require incompatible
choices. Weighted tree costs turn this into a positive-scale obstruction.

The general method uses any finite nonnegative rational kernel a_k on
cutoffs floor((2/3)^k X). The full-depth root weight is the cumulative sum
of a_k through that depth. Periodic local tolls give an exact convergent
series, normalized by sum a_k(2/3)^k. Eight adjacent cutoffs with

    a=[128,192,288,432,648,972,1458,2187]

give the strongest retained certificate. Independent full-root-cost arrays
through2^20 residue types reproduce its exact numerator3286041555236.
The simpler two-cutoff kernel gives the stronger conditional oscillation
bound for alpha-attainers.

The mechanism retains ONE common assignment across scales. It does not
determine the minimum ordinary or upper density. The two-step clauses fail
at horizon four, so these bounds do not transfer to longer horizons without
rebuilding their legal path constraints. [Full proof and experiment](crossroads233_20260926_flow.md).

## 2. An elementary height/rank alternative

**[THM-4493](../../01-canon/theorems/THM-4493-collatz-gcd-height-rank-certificate.md):
PROVED, independently audited.**

The user's triangle idea has a valid arithmetic realization. Given a
nonzero relation product a_i^c_i=1, choose i of maximal |c_i|. The triangle
inequality in each prime-valuation equation proves

    a_i divides product_(j!=i) gcd(a_i,a_j).

Thus strict dominance over that gcd product at EVERY node suffices for
multiplicative independence. This is a sufficient criterion, not an
equivalence.

For odd Collatz nodes m_0=n,...,m_(N-1), pairwise gcds divide exact subword
carries. If every node stays at least n, n>=N, and

    n>H_N=2^(N-1)(N-1)!3^(N(N-2))/2^((N-1)(N-2)/2),

then3m_0,...,3m_(N-1) are multiplicatively independent.
Here H_2=2, H_3=108, H_4=39366, H_5=86093442.

Equivalently, above this explicit height, the first N odd nodes must
either dip below their start or have full rank. Every sufficiently large
no-dip source has rank through

    N=floor(c sqrt(log_2 n)), for any fixed c<0.960047311978290... .

A proved horizon sidecar applies this to the usual no-descent set observed
through floor(log_2 n) shortcut steps. Any hypothetical infinite injective
positive orbit consequently has arbitrarily long consecutive independent
blocks and infinite total first-slot rank: use its future tail minima.
This does not prove independence of every prefix (G2) or force descent.

This is pointwise, whereas THM-4490 reaches a longer growing prefix for
almost all starts uniformly over intervals. The fixed-word cutoff is
elementary; the general subject of eventual independence under translation
already exists in the literature, linked in the [carry note](crossroads233_20260926_carry.md).
No publication-priority claim is made.

Both depth233 colliding paths have153 independent first slots, each proved
by private prime factors without factoring the large nodes. Thus rank
does not restore injectivity. Repeating only the two collision blocks
retains capacity two, by exact rational-base suffix decoding.

## 3. The precise tournament and metric reconstruction

[Graph proof and controls](crossroads233_20260926_graph.md):
**PROVED elementary reconstruction and fixed-word threshold theorem,
root-audited; FINITE-EXACT controls.**

Fix a marked path v_1 to ... to v_N. Map each edge{a,b} of K_(N-1),
a<b, to the chord{v_a,v_(b+1)}. This is an exact bijection of orientation
coordinates and gives

    binom(N,2)=(N-1)+binom(N-1,2).

It is not a vertex map: triangles and incidence change. A reduced backward
K2 lifts to a directed triangle. A complete weighted triangle with edge
lengths1,1,3 has unique edges and no loops yet violates triangle inequality.
The strongest line-metric repair requires one COMMON ordering for all
additive betweenness equalities; separate degenerate triangles do not suffice.

For a fixed Collatz word, retain chronological spine arcs and orient every
nonadjacent chord by its actual height-descent sign. The normalized affine
intercepts strictly increase. As the positive source height grows, backward
chords can only disappear; component boundaries open at exact maxima of
contact thresholds. There are at most N-1 component-splitting heights.

For the inherited231-to233 word, the only finite cut threshold is

    1441352971/5077565,

at cut seven. At source231 all18 positions form one strong component.
Every larger source in the exact cylinder231+2^28 t, t>=1, has exactly the
two components[0,...,6] and[7,...,17]. The endpoint changes from a rise to
a descent. This segment already dips to31 at position seven; it is not a
no-descent counterexample.

These are comparison cycles, not dynamical periodic orbits. The observer
forgets adjacent-step descent signs unless they are retained separately.

Concurrent work supplied an [insertion-slot interpretation](tournament_insertion_slots_20260926.md):
its slot count is exact, but the proposed Hamiltonian-path recursion fails
already on C3. The metric typing and slot-location wording were repaired
during integration. The two interpretations are compatible but different.

## Tests, synthesis, and remaining targets

Four retained programs run normally and with -O and require identical output.
The [verification runner](../../04-computation/experiments/crossroads233_20260926_verify.py)
writes the [artifact manifest](crossroads233_20260926_manifest.json).
The universes include33,867 fixed-spine tournaments;1,364 valuation words
and8,324 rational chambers;8,835 small integer tuples;66,405 no-dip
prefixes;12,007 finite-horizon checks; and independent depth20 tree costs.
Exact controls check the written proofs; the proofs supply their infinite
quantifiers.

The strongest remaining targets are:

- Determine the minimum upper or ordinary density in the two-step family,
  using finite positive cutoff functionals and matching global constructions.
  The best certificate in a bounded kernel search is not an optimum theorem.
- Rebuild the legal longer-horizon compatibility object for HYP-9140.
  Separate private or single-cutoff minima cannot be silently glued.
- Find a descent consequence of chronological rank that survives the
  depth233 merger and the infinite-rank alternative for divergent orbits.
- Seek a uniform height-sensitive statement across words. The present
  affine component filtration concerns one fixed word at a time.

The concurrent AMM constant improvement in THM-4494 was read as a stopping
signal for the numerical log_2(3) window analogy, not used as a dependency.
Collatz, HYP-9140, the uniform Robin comparison, and G2 remain OPEN.
