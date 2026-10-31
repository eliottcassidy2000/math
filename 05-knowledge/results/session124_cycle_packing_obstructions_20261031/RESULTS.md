# Session124: odd-cycle realization has a genuine packing barrier

## Status and scope

This packet proves a **restricted-family obstruction**, not a counterexample to
the ordinary Hamiltonian-path spectrum conjecture. In particular,29 is attained
by a tournament, but cannot be attained by any tournament whose directed cycles
are pairwise vertex-intersecting. The standard conjecture that only7 and21 are
missing remains open. No novelty priority claim is made for these elementary
structural consequences or the finite forbidden-pattern technique.

The motivating unrestricted identity is the proved odd-cycle collection formula

\[
 H(T)=\sum_{\mathcal C}2^{|\mathcal C|},
\]

where the sum ranges over collections of **literal directed odd cycles modulo
cyclic rotation**, pairwise vertex-disjoint. Distinct cycles with the same vertex
set remain distinct. Primary source: Irving–Omar, *A Survey of Rédei–Berge
Symmetric Functions via Matrix Algebra*, Corollary19,
<https://arxiv.org/html/2412.10572v1>, attributed there to Grinberg–Stanley,
<https://arxiv.org/abs/2307.05569>. The local canon is
`01-canon/theorems/THM-002-ocf.md`; the corrected header of
`THM-209-independence-polynomial.md` explicitly forbids collapsing supports.

## Theorem1: an exact relative spectrum through73, at every order

Let \(\nu(T)\) be the maximum number of vertex-disjoint directed triangles.
The all-order set of values \(H(T)\le73\) with \(\nu(T)\le1\) is exactly

\[
\{1,3,5,9,11,13,15,17,19,23,25,27,31,33,37,39,
45,47,49,51,53,61,63,65,69,73\}.
\]

Consequently the nine ordinary attainable values

\[
\boxed{29,35,41,43,55,57,59,67,71}
\]

cannot be attained within this class. Every tournament attaining any of these
nine numbers has at least two vertex-disjoint directed triangles. Equivalently,
its odd-cycle conflict graph has a nonedge and its OCF coefficient of \(4\) is
positive. It is therefore impossible to prove the full ordinary spectrum
conjecture solely by realizing complete odd-cycle conflict graphs.

### Proof and finite certificate

1. A tournament has a directed cycle iff it contains a directed triangle.
   Applying this to the vertex support of each cycle shows that the largest
   number of disjoint cycles, disjoint odd cycles, and disjoint triangles are
   equal. Thus \(\nu\le1\) implies every pair of literal odd cycles intersects,
   and \(H=1+2c_{\rm odd}\).
2. If there were two nonsingleton strongly connected components, each would
   contain a triangle, contradicting \(\nu\le1\). Removing singleton SCCs
   preserves \(H\) by ordinal-join multiplicativity. A nontrivial value therefore
   has a strong representative with the same packing restriction.
3. The cited Busch strong-tournament bound implies
   \(H^3\ge5^{n-1}\). Since \(5^8>73^3\), such a representative has \(n\le8\).
   This is an exact integer cutoff, not extrapolation from finite examples.
   Source: A. Busch, *A note on the number of Hamiltonian paths in strong
   tournaments*, Electronic Journal of Combinatorics13(2006), N3,
   <https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf>.
4. Choose any triangle \(a\to b\to c\to a\). Its complement is transitive:
   another triangle there would be disjoint. Order that complement
   \(v_1\to\cdots\to v_m\). All remaining edges are encoded by a word of
   length \(m\) over eight letters: bit \(d\) of letter \(i\) is1 iff core
   vertex \(d\) beats \(v_i\). Conversely every such word defines an actual
   tournament. Enumerate all words for \(m=0,\ldots,5\), rejecting those with
   two disjoint triangles. This is **37,449 presentations**, not an assumption
   about a canonical-isomorphism census. `check.py` independently counts
   endpoint-DP Hamiltonian paths and literal odd cycles on every accepted
   presentation, checks their OCF agreement, and obtains precisely the set above.
5. `global_witnesses.json` contains actual eight-vertex adjacency witnesses for
   each of the nine relative holes. `check_witnesses.py` counts every vertex
   permutation and separately verifies endpoint DP and OCF. These are positive
   ordinary-realizability certificates, not assumed spectrum completeness.

For29 there is also the six-vertex out-neighbor bitmask witness
`[50,16,3,7,12,30]` on vertices0,...,5. It has6 directed triangles,4 directed
5-cycles, and2 disjoint triangle pairs, so
\[
29=1+2(6+4)+4(2).
\]
There is no way to replace it by a tournament with exactly14 pairwise-intersecting
odd cycles, at any order.

## Theorem2: exact subsequence classification for the triangle-core family

Fix the directed triangle core and the eight-letter encoding above. There is a
finite minimal forbidden **subsequence** basis consisting of240 words of length3
and35 of length4, explicitly frozen in `patterns.json`. A word defines a
packing-at-most-one tournament iff it avoids every member of this basis.

**Proof.** Every triangle meets the core. Two disjoint triangles therefore use
at most four spine vertices. They cannot use fewer than three spine vertices,
since the core has only three vertices. Their induced presentation is a bad
word of length3 or4. Deleting further spine vertices until minimal gives the
stated form. The275 words are obtained by exhaustive enumeration of the512
length3 and4096 length4 words, dropping every length4 word with a bad
length3 subsequence. Conversely, a forbidden subsequence is an induced witness
of two disjoint triangles in the full tournament. `pattern_check.py` regenerates
the entire basis and independently checks the resulting recognizer against
actual triangle packing on all37,449 words through length5.

This is a complete all-order representation by a finite forbidden language,
not merely a census through a fixed order. It is generally not a unique
representation up to isomorphism.

### Exact freely repeatable alphabets

An alphabet \(A\subseteq\{0,\ldots,7\}\) satisfies \(A^*\) contained in the
legal language iff it is contained in one of these nine binary alphabets:

```
{0,1} {0,2} {0,4} {1,5} {2,3} {3,7} {4,6} {5,7} {6,7}.
```

Here is a direct proof, independently checked against the forbidden list.
If two distinct core coordinates vary across \(A\), choose separate pairs of
spine vertices showing a1 followed by a0 in each coordinate. This gives two
vertex-disjoint triangles, each using its own core vertex. Thus at most one
coordinate varies, so \(A\) is a singleton or one of the twelve edges of the
3-cube. Suppose the variable coordinate is \(a\), and the other core vertices
satisfy \(b\to c\). If their constant bits are \(b=0,c=1\), a spine vertex
forms a triangle with \(b,c\); two additional spine vertices realizing the
variable coordinate's1-to0 transition form a disjoint triangle through \(a\).
This excludes exactly three cube edges. Otherwise every triangle either uses
\(a\), or is the original core triangle (which also uses \(a\)); hence no two
are disjoint. The nine remaining edges are precisely those listed.

### Global finite-phase restriction and entropy bound

Store which proper prefixes of the275 forbidden words have appeared as
subsequences of the current input. There are84 distinct proper prefixes,
including the empty word, hence83 nonempty Boolean state coordinates. The
state only increases. Reject if appending a letter completes a forbidden word.
There are at most83 strict state updates in any accepted word.

At a fixed state, every self-loop letter can repeat arbitrarily in any order.
Its loop alphabet must therefore be contained in one of the nine binary
alphabets. It follows that every legal word consists of **at most84 binary
free phases separated by at most83 exceptional single letters**. This does
not say there are only84 unseparated blocks, nor that arbitrary independent
choices of the phases are legal.

The number \(L_m\) of accepted length-\(m\) presentations consequently obeys

\[
L_m\le 2^m\sum_{j=0}^{\min(83,m)}\binom mj4^j
=O(m^{83}2^m).
\]

To count, specify the positions of the \(j\) strict updates; their letters have
at most8 choices, and every remaining position has at most2 choices. The
automaton is deterministic, so this overcounts rather than misses words.

**Boundary:** this is a restriction on packing-one constructions. It does not
bound the unrestricted tournament spectrum and does not prove a density-zero
claim for the set of attainable \(H\)-values.

## Theorem3: bounded-packing entropy principle

For each fixed positive integer \(k\), tournaments with \(\nu(T)\le k\)
have a finite core-and-word representation whose legal languages have
\(O(m^{P_k}2^{km})\) words of spine length \(m\), for a finite effective
constant \(P_k\). In particular their number of isomorphism classes on \(n\)
vertices is at most \(O(n^{P_k}2^{kn})\). This is only an upper bound; word
representations can repeat isomorphism classes.

**Proof.** Choose a maximum packing of \(r\le k\) triangles. Their \(3r\)
vertices form a core; the complement is transitive. The alphabet size is
\(2^{3r}\le2^{3k}\). For each fixed core, forbidden packings of \(k+1\)
triangles use at most \(2k+2\) spine vertices, because every triangle meets the
core. Thus legality is characterized by finitely many forbidden subsequences
of bounded length. A monotone proper-prefix automaton has at most
\[
P_k\le\sum_{i=1}^{2k+1}(2^{3k})^i
\]
nonempty coordinates (this loose universal bound suffices).

If a freely repeatable alphabet varies in \(k+1\) core coordinates, choose for
each one two distinct spine occurrences with bit1 then bit0. They form
\(k+1\) vertex-disjoint triangles. Therefore any recurrent alphabet varies
in at most \(k\) coordinates, so has at most \(2^k\) letters. The same strict
update counting gives
\[
L_m\le2^{km}\sum_{j=0}^{\min(P_k,m)}\binom mj2^{2kj}.
\]
There are only finitely many possible cores for fixed \(k\), and summing over
\(r\le k\) proves the isomorphism-class bound. The entirely transitive case
\(r=0\) adds one class per order. OCF in this class has degree at most \(k\),
but **not every degree-\(k\) independence polynomial is thereby realizable**.

## Verification / provenance

All code uses the Python standard library and explicit exception checks,
not assertions that vanish under optimization. Run from any directory:

```
python /path/to/package/check.py
python /path/to/package/pattern_check.py
python /path/to/package/check_witnesses.py
```

Each has a matching `python -O` receipt. The finite spectrum theorem does not
rely on external data. Positive nine-value witnesses were extracted from
baseline `f5b6c85b8c1d905ac51b40f82649321769a1a8d8`, path
`05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/n8_classes.csv`,
and are rechecked from their actual adjacency, without trusting its H column.
The exploratory canonical census files were removed from this packet.

The earlier repo result
`05-knowledge/results/h63_n8_isoclass_census_s10.out` already showed that both
order8 witnesses with H63 have complete conflict graph K31. Thus63 is not a
packing barrier. The first relative barrier identified here is29, not63.

## Reconciliation with the older “poison-factor” exploration

The historical scripts `04-computation/k3_poison_complete.py` and
`04-computation/k3_poison_fast.py` explored the factor \(1+3x\), including
\(I(P_4,x)=(1+x)(1+3x)=I(K_1\sqcup K_3,x)\). Those experiments concern arbitrary
simple graphs and must not be read as a realization theorem for tournament
odd-cycle conflict graphs. The complete-script output is empty. The fast-script output records only
1,099 graphs on at most5 vertices, but its historical conclusion extrapolates
to every graph. That extrapolation is false: \(K_{10}\) has value21 but
polynomial \(1+10x\), not divisible by \(1+3x\). Its phrase that \(K_3\) cannot
appear in a conflict graph is also invalid if read as forbidding a subgraph: an
actual conflict graph \(K_7\), for example, contains many such subgraphs. Only
being the *entire* conflict graph \(K_3\) is ruled out by H7. The unrelated
`04-computation/h21_poisoning_graph.py` uses “poisoning” for a vertex-deletion
relation in triangle-packing-two tournaments; that is a different construction.

Our result gives a particularly direct obstruction to a poison-only criterion:
\[
I(K_{14},x)=1+14x
\]
has no factor \(1+3x\), and its value at2 is the ordinarily attainable number29,
yet **\(K_{14}\) cannot be the literal odd-cycle conflict graph of any
tournament**. If it were, that tournament would have packing number1 and H29,
contradicting Theorem1. Thus neither avoiding the old factor nor obtaining an
allowed numerical value certifies conflict-graph realizability. In contrast,
the actual six-vertex H29 witness has polynomial \(1+10x+2x^2\).
