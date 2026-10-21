# Order twelve has no strong self-converse deficit-six center

October 21, 2026. **FINITE-EXACT complete exclusion**, with a structural
normalization and a no-search literal certificate consumer. Every necessary
normalized order-twelve source is either not self-converse or has at least
eight black Hamilton-path orbits. Thus

\[
[z^{12}]R_6(z)=0.
\]

Combined with the [proved cutoff and order-ten exclusion](../deficit_six_centers_20261021b/COMBINED_CUTOFF.md),
this packet alone leaves the nine orders

\[
\boxed{14,15,18,20,21,24,25,28,35}.
\]

Those coefficients remain **OPEN in this packet**. Later independent results
may remove further orders. The previously reported list had ten orders; here
the reduction to nine is a fresh complete exclusion of twelve, not an assumed
correction to that list. In particular

\[
R_6(z)=z^5+z^6+z^7+O(z^{14}),\qquad\deg R_6\le35.
\]

The minimal eventual recurrence order eighteen and guaranteed onset forty-eight
are unchanged. This is not a census of ordinary isomorphism classes, and it
does not compute the exact deficit of every normalized source.

## 1. Inheritance and the certificate being retained

For a literal tournament $T$, let $H$ count directed Hamilton words, let
$A=\operatorname{Aut}(T)$, and put $h=H/|A|$. A word is **blue** when its
positional reversal is an anti-automorphism of the same source. Let $b$ count
blue word orbits. For a self-converse source, the deficit is $\Delta=h-b$.
The target is a strong self-converse source with $\Delta=6$.

The nearest proved mechanism is the marked-cycle index and the incoming
score-run obstruction. The exact cycle-index argument is inherited from the
[cycle-cut source](../paired_black_leaf_reduction_20261020/CYCLE_CUT_NEAR_PURITY_OBSTRUCTION.md)
and restated below. The [odd-shore proof](../odd_shore_deficit_centers_20261021/session22-centers/PROOF.md)
shows how full-score strings distinguish new prefix orbits; this packet uses
the same invariant on complete literal words rather than a single run.
The [order-ten producer and receipt](../deficit_six_centers_20261021b/README.md)
are the closest finite normalization. Their full automorphism enumeration is
not needed at order twelve.

The canonical hostile is Paley7: every vertex has score three, so every path
has the same score string, although the source has black paths and deficit
six. The corrected near miss is to identify equal score strings with equal
path orbits, or to identify a chosen period with the exact cut index. Neither
implication is used. The least-used sidecar is the entire Hamilton word with
its source-owned complement target.

The live concepts are the actual Hamilton cycle, its odd rotation stabilizer,
the full-score word, positional reversal, and the literal target. The source
binary relation is the tournament arc itself; score ties remain ties and are
not oriented into a second tournament.

The map used by the certificate is

\[
P=(v_0,\ldots,v_{n-1})\longmapsto
\sigma_T(P)=(d_T^+(v_0),\ldots,d_T^+(v_{n-1})).             \tag{1}
\]

It preserves a sufficient test for distinct source-automorphism orbits and
for blackness. It loses adjacency, vertex labels, exact orbit equality and
the target graph. The bank therefore retains each actual source and word;
the target is reconstructed losslessly and compared with the independent
producer's complete target stream.

## 2. Why all order-twelve candidates occur in exactly $2^{18}$ frames

Fix a Hamilton cycle and label it

\[
0\to1\to\cdots\to11\to0,
\]

with successor $R$. Let $m=|A\cap\langle R\rangle|$ and $q=12/m$.
The automorphism group of a tournament has odd order: an involution would
reverse the arc joining a transposed pair. Consequently $m$ is an odd divisor
of twelve, hence $m=1$ or $3$.

The twelve cuts of the cycle meet exactly $q$ path orbits. Indeed the unique
vertex map taking one ordered cut to another is a power of $R$. If one cut
is blue, the product of its positional reflection with another cut's
reflection is $R^{2t}$. Thus the blue cut orbits solve one congruence
$2t=0\pmod q$ and number $\gcd(2,q)\le2$. If none is blue, all cut orbits
are black. Therefore

\[
q-2\le\Delta=6,\qquad q\le8.
\]

This rules out $m=1,q=12$, leaving precisely

\[
\boxed{m=3,\quad q=4,\quad R^4\in A}.                    \tag{2}
\]

This order-twelve reduction needs neither the odd-order closing-cycle theorem
nor a score-regularity assumption. Strongness supplies an actual Hamilton
cycle; all subsequent coordinates refer to that cycle.

The action of $R^4$ on the 66 unordered vertex pairs has 22 orbits, each of
size three. No orbit can reverse its own orientation, since the acting group
has odd order. Its two orientation choices are independent of those for the
other orbits. The twelve prescribed cycle arcs occupy exactly four such
orbits, forcing their orientations. The remaining eighteen orbits give

\[
2^{22-4}=262\,144                                             \tag{3}
\]

literal normalized frames. Each has the displayed directed Hamilton cycle
and is therefore strong. Each is invariant under $R^4$. Distinct masks give
different labelled adjacency matrices. Different masks may still be
isomorphic; no quotient or class count is asserted. Every candidate from
(2) is represented after relabelling its marked cycle.

The precise mask gauge is reproducible. For each unordered pair, take the
lexicographically least unordered pair in its $R^4$ orbit. Order these 22
representatives lexicographically; skip the four forced orbits; read free
bits in that order. Bit one directs the representative's smaller label to
its larger label, and the other two arcs are its rotations. The independent
consumer derives each pair's canonical representative, rather than using
the producer's visited-edge orbit sweep.

## 3. A sufficient black-orbit certificate without automorphism search

Every automorphism preserves full outdegrees. Thus two Hamilton words with
different strings (1) belong to different $A$-orbits. If $P$ is blue, its
positional reversal is an anti-automorphism, so

\[
\sigma_T(P)_i+\sigma_T(P)_{n-1-i}=n-1
\quad\hbox{for every }i.                                  \tag{4}
\]

The contrapositive is sufficient: failure of (4) at even one coordinate
certifies that this actual Hamilton word is black. Consequently eight actual
words with eight distinct strings, each failing (4), certify

\[
\Delta\ge8.                                               \tag{5}
\]

The certificate itself works for any source; self-converseness is only needed
to place it in the center classification. No equality case or necessity is
claimed for the score observation.

Before reading paths, there is a cheaper self-converse guard. A source with
an anti-automorphism has equally many vertices of scores $d$ and $n-1-d$ for
every $d$. A histogram violating this symmetry is therefore not self-converse.
Passing this guard does not prove self-converseness; all passing frames are
over-included in the path bank.

The complete partition is

| Exact branch | Normalized frames | Authenticated conclusion |
|---|---:|---|
| Non-complementary full-score histogram | 111,104 | Not self-converse |
| Complementary histogram, eight supplied black score words | 151,040 | At least eight black path orbits |
| Unresolved | 0 | None |

This proves the order-twelve exclusion. It does not merely filter a random
sample or trust a bounded search failure: every admitted source has a finite
positive exclusion certificate.

## 4. Producer, source binding, targets and independent consumer

`census.cpp` searches the 151,040 score-balanced frames. It branches on score
strings. For a fixed prefix string it retains reachable states
`(used_vertex_mask,last_vertex)` with a literal path witness. Merging two
states with the same used vertices and last vertex is valid because their
remaining extension problem is identical. A completed nonpalindromic string
is saved, and search stops after eight distinct strings. This search used
9,135,120 prefix nodes and about four seconds in the recorded local run.
Neither timing nor prefix-node count is a premise of the proof.

For every saved word, the producer also writes the complete target graph:
retain each of the eleven path arcs and reverse every other source arc. The
target uses the same twelve vertex labels. The expanded stream has
1,208,320 words and targets and 18,124,808 bytes, with SHA-256

```text
e80a2fde831fdca7cc9ad43d31a024d58b267cab1b053b53d64cf85f9c755132
```

`pack_bank.py` keeps the six-byte word representation and a commitment to
the complete expanded stream. Since the target is a deterministic function
of the retained source and word, this is lossless compression. The stored
`inputs.json` is about 2.35 MB. It is ordered by the actual frame masks,
skipping only the independently rechecked histogram failures.

`check.py` performs **no Hamilton-path search, automorphism enumeration or
source replacement**. It reconstructs all frames, authenticates the actual
cycle and rotation, validates every supplied permutation and each successive
arc, recomputes the score tests, and reconstructs all targets. The complete
expanded stream must match the independent producer's hash. The positive
exclusion argument rests on the individually replayed paths and score
predicates, not on a hash standing in for mathematical checking.

The public `audit_record` binds a record to its exact mask and immutable
source rows. The generic `score_ports` consumes supplied `(word,target)`
pairs and exports a sufficient lower bound. It never discovers paths or
accepts a caller's asserted score or blue status. Exact integer/type checks
precede indexing or comparisons; floats, bools, negative indices, duplicate
vertices, altered sources, duplicate score words and wrong targets reject.

The smallest positive control used here is the known six-vertex center with
rows `(36,17,10,19,37,14)`. An independent enumeration of all 720 permutations
gives

\[
(H,|A|,H_{\rm blue},\Delta)=(45,3,27,6).
\]

It has exactly six distinct nonpalindromic score strings, accepted as a
six-port certificate; the eight-port request rejects. Paley7 and its literal
black word `(0,1,2,4,5,6,3)` show that a constant score string can hide real
blackness. The displayed order-twelve sample also replays its positional
reversal directly and constructs targets by a separate pairwise method.

## 5. Reproduction and remaining obligation

Immutable verification needs only the Python standard library and the saved
bank:

```text
python -B 05-knowledge/results/deficit_twelve_20261021c/check.py --verify
python -B -O 05-knowledge/results/deficit_twelve_20261021c/check.py --verify
```

To reproduce the separate producer, place the executable and temporary raw
bank in a scratch directory, then run:

```text
g++ -std=c++17 -O3 05-knowledge/results/deficit_twelve_20261021c/census.cpp -o <scratch>/census
<scratch>/census 262144 <scratch>/bank.bin
python -B 05-knowledge/results/deficit_twelve_20261021c/pack_bank.py <scratch>/bank.bin
python -B 05-knowledge/results/deficit_twelve_20261021c/check.py --verify
```

The final receipt records 5,357,643 always-active gates, the complete frame
partition, hashes of all literal sources and ports/targets, the positive and
hostile controls, and LF pins for the inherited proof/receipt inputs and the
producer. Normal and optimized immutable replays both pass and exactly match
the saved receipt, whose LF SHA-256 is
`fe54716a7b6e6fac30a41e118e9a59e2e80cdea0da6e7602932cfa0697779da2`.
Independent family proof/API review passed: complete frame normalization,
score-histogram direction, eight-black-orbit implication, source-bound target
consumer and positive/hostile boundaries. That review did not duplicate the
full bank replay; the normal/optimized replays above are the producer lane's
recorded runs.

The new proof step is the complete order-twelve exclusion through an
automorphism-invariant, cheaply distinguishable observation with literal
path witnesses. It leaves the nine orders above, including any stronger
structural obligations supplied elsewhere for 25 and 35. It does not infer
vertex transitivity merely from the old regularity condition. The bound 35
and recurrence onset 48 do not improve when the coefficient at twelve is
removed.
