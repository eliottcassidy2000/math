# Exact infinite-source readiness, and an old-arc repair that changes the architecture

**PROVED using complete finite target-capped automata.** The scalar window is
finite; source word lengths, neutral padding and suffix lengths are unbounded.
This distinction is the point of the experiment. No all-h readiness theorem or
Chan–Pak closure is claimed.

## 1. The full target language supplies an exact inverse constructor

`selector.py` closes and minimizes the entire language H=N, then searches the
finite machine for an accepting path. If none exists, it certifies absence in
the entire two-apex grammar, not merely absence below a chosen order. If one
exists, its word produces literal tournament rows and is independently counted
by subset DP. `SELECTOR63.json` contains the full39-state minimal transition
table and the executable accepting path

```
word 100122,
states 0 -> 1 -> 4 -> 14 -> 24 -> 32 -> 37,
outgoing masks (190,252,248,176,96,64,137,48), H=63.
```

The same procedure returns the one-state empty languages for7 and21. This is
an exact per-target selector, not a uniform positive answer for arbitrary N.

The comparator retains the intended target, not only scalar count. For example

```
H(120)=H(212)=13,
H(12011)=41,       H(21211)=59.
```

The suffix11 is a shortest suffix distinguishing those prefixes for target41.
The words are therefore distinct exact target-continuation types even though
both scalar H and their smallest realization order agree.

## 2. Complete one-ended binary-menu classification through source count63

An append-only binary fan from a source word w of count h consists of two
arbitrary suffixes x,y such that

```
H(wx)=2h-1,        H(wy)=2h+1.                         (1)
```

The two suffixes need not have the same length. They may depend on all the
source information. They retain the named apex pair and old spine order/arcs.

**Theorem R.** Among attained odd source counts `3<=h<=63`, a source supporting
(1) exists exactly at h=5 and h=9. There is no depth-two such fan from any
attained `3<=h<=31`: both first-stage physical children cannot themselves
support their respective two branches.

This is not a theorem for all source h. It permits arbitrarily long source
words and arbitrary suffix lengths. It does not constrain insertion inside
the old spine, changing marks, old-arc repair, or rebuilding a new source.

**Complete finite proof.** Close the capped36-state transfer at target127.
It has28016 live capped states. Minimize with the entire scalar output; the
result has1080 states including overflow. Its non-loop transition graph is
acyclic. For each state compute its full descendant-state set, including all
possible suffixes, by reverse topological propagation. Group states by their
actual h. A state is ready precisely when its descendant set contains both
output values in (1). For depth two require those two descendant states to be
ready in turn. Every intermediate count is bounded by its terminal count by
induced monotonicity, so the cap discards no successful path. All first-stage
targets for h<=63 lie at most127; every second-stage target for h<=31 does too.
No assertion is made about depth two for h>31 (`null` in the receipt).
`readiness.py` implements exactly this proof; `READINESS.json` records every
source count's entire set of branch masks and explicit positive controls.

Positive common-source controls are

```
12 (H5)  -> 121 (H9),   122 (H11);
102 (H9) -> 1021 (H17), 1022 (H19).
```

A nontrivial loss-of-menu example occurs at h=23, with both target values
individually possible from some source of that same h:

```
1321 (H23) -> 13211 (H45),
1222 (H23) -> 12222 (H47).
```

Nevertheless NO single H23 source, at ANY word length, admits both branches
by arbitrary appending. This is stricter than checking the minimal sources
and is a concrete reason why replacing witnesses solely by scalar equality
breaks an induction. The target-specific machines expose exactly that loss.

## 3. Interior insertion alone does not repair H13 -> H27

The complete H13 normal form in `PROOF.md` and inherited Session191 is

```
0* c 3*,  c in {120,212,230,302,310,321}.              (2)
```

There is no two-apex word of count27 containing ANY of these six core words
as a subsequence. Therefore arbitrary interior insertions, not just appending,
cannot grow an H13 source to H27 while retaining the apex pair and old arcs.

**Exact automaton certificate.** Take the47-state minimal exact-H27 machine and
its product with the four-state subsequence recognizer for each three-letter
core c. The latter tracks how many letters of c have been greedily matched;
it accepts at3. The product's reachable state counts for the six cores are
respectively56,57,58,59,59,57. No reachable state accepts both predicates.
This finite closure proves the assertion at every length. Every source in(2)
contains its core, so retaining the source as a subsequence would retain that
core as well. `insertion.py` records the certificate in `INSERTION.json`.

The conclusion is deliberately narrower than an arbitrary induced-subgraph
claim for all tournament families. The certificate retains the fixed two
apexes and word embedding. It also does not say H27 itself is absent.

## 4. One old-arc repair opens the missing branch

The failed architecture points to a concrete change, not just a negative
conclusion. Start from the actual H13 word120. Its lower binary target is
obtained by a single physical prepend:

```
120 -> 1120,      13 ->25=2*13-1.
```

For the upper target add a neutral source first. The word0120 has the rows

```
(62,44,24,48,34,4), H=13.
```

Here vertex0 is new and vertices1,2,3,4,5 are old. Reverse the OLD arc1->2,
then reverse the arc0->2 incident with the new vertex:

```
(62,40,26,48,34,4), H=15;
(58,40,27,48,34,4), H=27=2*13+1.                       (3)
```

The final spine is still transitive after reordering its names as(2,0,1,3),
with the same apexes4,5. It is the two-apex word2010. Thus this is an explicit
small architecture-changing binary repair: one old-arc change is what the
all-order insertion-only obstruction forbids us to omit. Every intermediate
count and the final reframe are independently checked by literal graph DP.

This repaired fork is NOT automatically renewable. A cheap hostile extension
repeats the initial1 block in the source to1120 (H25); the analogous two-edge
repair of01120, explicitly reversing(2,3) and(0,3) so that the
last initial1 vertex and the new source reverse their arcs to the letter2
vertex, givesH71, not51. Keeping the original numerical arc indices instead
is a different operation. Growing a physical pattern changes the full
response, not just the scalar input. A uniform old-arc repair law with its
preserved interfaces remains the positive next obligation.

## 5. Why this is not merely a larger finite spectrum check

The outputs classify infinite-source continuation languages, prove that all
fixed-target recurrence is neutral padding, give exact per-target inverse
witnesses, and locate precisely which information a binary induction loses.
The core positivity problem remains unbounded in N. The next constructive
route should search old-arc rewrites or changing apex interfaces that preserve
these successful-target languages, rather than adding more suffix depth to
an already complete failed append-only menu.
