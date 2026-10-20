# Tournament metagraph growth: exact recursions, multiplicative strata, and the missing state

**20 October 2026 — PROVED scoped constructions / FINITE-EXACT independent
census / OPEN compressed all-size class-count recurrence.**

The user's counting identity is correct. The substantial remaining question
is which information lets the structures grow without repeatedly solving
the whole isomorphism problem. This packet recovers the existing exact
additive recursion, adds a checked multiplicative decomposition and its
line-transport correction, and tests several smaller proposed states.
It does not claim a faster general formula for A000568 or a new proof of
any unrelated open problem.

## 1. What is counted, and what is already known

Let a_n count ordinary tournament isomorphism classes, s_n self-converse
classes, M_n converse-merged nodes, and P_n merged nodes containing a pair.
Then, for every n,

    a_n = M_n + P_n = 2M_n - s_n.

Fixing a directed Hamiltonian path leaves m_n=binom(n-1,2) free arcs.
The full tiling set X_n has 2^m_n elements. An individual ordinary class T
has fibre H(T)/|Aut(T)|, rather than 2^m_n elements. Summing those fibres
gives 2^m_n. A merged non-self-converse node combines two such fibres.
These are inherited from
[THM-1430, the tiling/class dictionary](../../../01-canon/theorems/THM-1430-the-tiling-class-metagraph-dictionary-and-which-tricks-pay.md);
no novelty is claimed for them.

Two graphs must be distinguished:

- The arc-flip graph connects classes related by ONE arc reversal.
- The blue/black line multigraph joins the nodes containing t and tau(t),
  where tau complements ALL free tiling bits. It retains loops and
  parallel line instances. Grid reflection sigma represents tournament
  converse after relabelling; blue means sigma(t)=t.

The exact common object is the span M_n <- X_n -> L_n, where
L_n=X_n/<tau>. The two involutions commute on X_n, but tau is not a
function on M_n. A tournament class admits several path presentations,
whose complements can lie in different classes.

The independent [census](enumeration.md) recovers the following without
using Burnside's formula in class generation:

| n | a_n | s_n | M_n | P_n | total tilings |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 0 | 2 |
| 4 | 4 | 2 | 3 | 1 | 8 |
| 5 | 12 | 8 | 10 | 2 | 64 |
| 6 | 56 | 12 | 34 | 22 | 1024 |
| 7 | 456 | 88 | 272 | 184 | 32768 |

The a_n values agree with [OEIS A000568](https://oeis.org/A000568).
The sequence entry is a comparison source, not an input to generation.

## 2. The smaller graph is exact, but its labels matter

[THM-830, the gap-tournament suspension](../../../01-canon/theorems/THM-830-b3-deletion-deck-mirror-current-calculus.md)
already gives a precise form of the user's small-graph idea: a size-n
tiling is a labelled tournament on the n-1 gaps of its distinguished path.
There are exactly binom(n-1,2) gap-pair bits.

On that gap tournament tau is genuine converse. On the original
n-vertex tournament it reverses only non-path arcs. These are different
operations on different objects. Forgetting the order of the gaps loses
the suspension information: isomorphic small gap tournaments need not
produce isomorphic large tournaments. Thus the useful carrier is a
gap tournament with its ordered path context.

## 3. Add or remove one or two vertices: recover the existing theorem

[THM-796, three-sorted recursive incidence](../../../01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md)
already proves the exact all-size construction

    X_n = (X_(n-1) x_[X_(n-2)] X_(n-1)) x {apex bit}.

The two faces must agree on their common interior. For n>=5, replacing
the two tiling faces by bare complement lines loses exactly one coherent
endpoint-phase bit. It is therefore incorrect to glue independently
chosen node or line labels without overlap information.

Appending one path endpoint gives 2^(n-1) lifts from X_n to X_(n+1).
Reflection exchanges left and right deletion. Adding one endpoint at
EACH side gives a particularly clean n -> n+2 law:

| parent line, n>=3 | blue child lines | black child lines |
|---|---:|---:|
| blue | 2^n | 2^(2n-1)-2^n |
| black | 0 | 2^(2n-1) |

This follows because the 2n-1 new boundary bits consist of n-1 reflected
pairs and one reflection-fixed apex. A nonzero interior reflection defect
survives every extension. These colour totals are already proved in
THM-796; this packet independently checks them. They do not determine
which upper isomorphism nodes receive the child lines.

The same theorem already proves why multiplying node-only transition
matrices fails: at n=5, lift multiplicities vary within 11 of 19 nonzero
node blocks. A uniform distribution within a node is not preserved.
Later finite codecs must retain their scope: THM-809 succeeds at n=8,
but [THM-828](../../../01-canon/theorems/THM-828-n9-s2-parity-syndrome-join-reduction.md)
finds 58 reflection doubletons at n=9 for its raw count address.

## 4. Multiplication of sizes: a canonical integer scale

Write D_k(T)=T[TT_k]: replace every vertex by an ordered transitive
k-chain and replace each base arc by all k^2 cross-arcs in its direction.
The construction itself is standard in the repository. The useful
additional coordinate is the following canonical scale.

Join two vertices when they have identical relations to every OTHER
vertex. The resulting homogeneous-pair graph is a disjoint union of
directed paths; its components are precisely the maximal transitive
modules. This decomposition is inherited from
[THM-4463, homogeneous contacts](../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md).

Let their lengths be l_1,...,l_b and define g(T)=gcd(l_1,...,l_b).
The [growth proof](growth.md) establishes:

    U is D_k(T) for some T  iff  k divides every l_i(U).

When it exists, the root T is unique up to isomorphism. The same holds
after converse merging. D_k preserves the reduced quotient and multiplies
every chain length by k; division recovers the unique root.

If p_n counts classes with g(T)=1, the resulting exact divisor law is

    a_n = sum_(d|n) p_d,
    p_n = sum_(d|n) mu(d) a_(n/d).

The checked values p_1,...,p_6 are 1,0,1,3,11,54. "Primitive" here means
only chain-gcd one. It does not mean strongly connected or modularly prime.
This is a structural decomposition, not yet a recurrence computing p_n
from smaller orders: new primitive classes remain to be classified.

### The exact correction needed for black lines

In the concatenated block-path presentation, D_k does not commute with tau.
There is, however, a deterministic mask E_(n,k) such that

    tau_(kn) D_k(t) = D_k(tau_n(t)) XOR E_(n,k).

Its support consists exactly of:

1. Non-path arcs internal to a clone block.
2. Cross-arcs between consecutive base blocks, except their path connector.

Consequently

    |E_(n,k)| = n*binom(k-1,2) + (n-1)*(k^2-1)
              = (k-1)*(3nk-2k-2)/2.

The mask is reflection-invariant. Correcting tau by this mask gives an
exact line transport on the block-respecting presentation subcube.
Ordinary tau can leave the inflated family: D_2(TT_2)=TT_4 has chain
length 4, whereas its complemented free-bit presentation has chain
lengths 1,1,2 and no doubling root.

Even the corrected map does not cover an inflated node's entire tiling
fibre. The directed triangle has one presentation modulo automorphisms;
D_2(C3) has 15. Interleavings of the clone blocks supply the other
presentations. A block partition plus its path order is indispensable.

Even a path's block word and run lengths can lose too much. In D_2(C3),
the paths A1 B1 B2 C1 C2 A2 and A2 B1 B2 C1 C2 A1 have the same word
ABBCCA. Their ordinary tau partners have degree multisets
(2,2,2,3,3,3) and (1,2,2,3,3,4), so they are different even after
converse merging. Which ordered clone occupies each separated visit is
an additional load-bearing coordinate.

## 5. A small hostile makes the missing growth state visible

At n=3 the transitive and cyclic tournaments are the two ends of one blue
line. Both merged fibres have size one. The unmarked coloured multigraph
has an automorphism swapping them, and eight natural local counts agree.

Yet the actual one-end extension counts to the three merged n=4 nodes
(transitive, non-self-converse pair, other self-converse) are

    transitive parent: (1,1,2),
    cyclic parent:     (0,1,3).

Thus the unmarked coloured metagraph alone does not specify growth.
Marking the transitive seed, or retaining H and |Aut|, separates this
particular pair. This witness does not rule out every enriched state.

The stronger finite probe at n=7 finds nine pairs with the same eight
local coordinates AND the same H and |Aut| but different line-partner
fibre profiles. One such pair has fibre 38, no blue tilings, H=19,
|Aut|=1, line-support degree 16, and arc-flip degree 20. One has two
black incidences to fibre-41 nodes; the other has none.
The full exact codes, coordinates and profiles are in
[enumeration.out](enumeration.out).

The actual closure test for a proposed observer O and operation F is

    O(x)=O(y)  implies  O(F(x))=O(F(y)).

This is the inherited operation-specific congruence test from
[THM-840](../../../01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md).
Matching a few counts or distinguishing today's nodes is weaker than
passing this test for future extension, deletion and line transport.

## 6. What can already replace Burnside, and what remains

Class BFS is already an exact alternative: start at the transitive class,
reverse each arc, canonicalize, and continue until no class is new.
Connectedness of the labelled orientation cube proves completeness.
Our independent implementation visits 456 classes at n=7 using 9576 arc
trials; it need not construct the full labelled cube or retain metagraph
edges. This revalidates THM-1430's algorithm, rather than inventing it.

For the count alone this is not a demonstrated improvement over Davis's
formula, whose terms run over odd partitions. Enumeration still visits
every class. The desired advance would be a smaller structural grammar
with proved completeness and operation-compatible states.

The strongest concrete next targets are:

- Use the unique chain root and its weight vector to separate multiplicative
  copies from genuinely new primitive classes; retain the block-path
  interleavings and clone allocation when transporting lines.
- Classify the primitive weighted gap cores and their extension profiles.
  Test closure before replacing them by node degrees, colours or H.
- Exploit already successful restricted grammars. For example,
  [THM-2454's centre law](../../../01-canon/theorems/THM-2454-pure-blue-alphabet-and-center-law.md)
  factors blue fibres of palindromic strong-component stacks. Its proposed
  complete atom alphabet remains an open classification beyond its finite
  censuses. A proved complete alphabet would yield an actual counting
  grammar for that stratum, rather than only a total-mass identity.

This is also the precise methodological connection to the earlier Collatz
rank work: retain the information required by the next operation. Here
lost phase and path incidence are exposed by finite counterexamples.
No Collatz conclusion is inferred from the tournament constructions.

## 7. Research and verification record

Anchor: exact cross-size generation. Niche: canonical transitive-chain
scale. Wildcard: the gap-tournament suspension and its operation congruence.
The concept board was class/converse sheet, path fibre, coherent phase,
reflection defect, chain scale, and extension response.

Closest proved mechanism: THM-796's pullback. Canonical hostile: tau and
endpoint deletion fail to descend already at n=4. Corrected near miss:
a complete static fingerprint need not be an update state. Least-used
sidecar: canonical chain lengths together with block-interleaving paths.
Meta-patterns used: search the statement before the method; retain local
incidence through counting; correct the object before sharpening technique.
No new universal meta-pattern is promoted from this one thread.

The source/target maps and their losses are explicit above. Independent
normal and optimized replays, full-permutation canonicalization controls,
literal finite universes and artifact hashes are recorded in
[VALIDATION.md](VALIDATION.md). The old
[metagraph explainer](../../../03-artifacts/explainers/09-metagraph.md)
was repaired: incorrect class counts, a false mixed-triangle prohibition,
and an unsupported greedy claim are withdrawn with lineage in MISTAKES.md.

Concurrent THM-4729/4730/4731 reservations were not used as proved
dependencies. This packet reserves no scarce theorem IDs.
