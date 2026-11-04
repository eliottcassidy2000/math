# A name-preserving, all-strong native +4 repair of the first +4 sibling

2026-11-03. **PROVED — independent analytic/source/replay audit passed.**
The inherited guard compiler is PROVED; the displayed programs and
local strongness certificates are FINITE-EXACT. This note gives one all-clock
map. Its bulk output has the separately audited continuation in
[PREPARED_FOUR.md](PREPARED_FOUR.md); two-branch closure remains open.

## 1. Exact theorem and inherited mechanism

For every integer `t>=1`, let the source have physical names

```text
p0,...,p(t-1), b0,b1,b2,b3,b4, r.
```

The `p` vertices form a transitive prefix in that order and beat all `b`
vertices. The `b` vertices form `TT5`; `b0,b2,b4` beat the same physical root
`r`, while `r` beats every other vertex. This is the inherited first `+4`
sibling, with

```text
H(source_t)=21*2^t-11,       F_r-H=12.
```

**Theorem.** There is an explicit algorithm for every such `t` giving a strong
named target on exactly the same old vertices, with

```text
H(target_t)=3H(source_t)+4=63*2^t-29.
```

The algorithm uses only temporary helper births, legally guarded native
Hamilton-path complements, and deletion of those helpers. Every old name,
including `r`, survives unchanged. **Every intermediate tournament is strong.**
There are no permanent births or old-vertex deletions. At `t=1` the budget is
9 helpers and18 complements; at every `t>=2` it is11 helpers and22 complements.
At most one helper is present at a time, and peak order is exactly `t+7`.

The target graph recipe has a uniform bulk form for `t>=3` and two explicit
finite patches. This piecewise algorithm has the stated all-clock count
identity; it is not a single unchanged suffix representation at the two
boundary clocks.

The direct guard construction is inherited from
[`recursive_seed_switching_20261103/transport/FOLLOWUP.md`](../../recursive_seed_switching_20261103/transport/FOLLOWUP.md),
section “A direct guard avoids subset search.” The two-native-move helper
gadget is inherited from
[`adaptive_metagraph_reformulation_20261103/METAGRAPH.md`](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md),
§5. Both apply to **any** named tournament and pair; they require neither a
Hamilton count nor strongness. The new strongness proof below concerns this
particular family and schedule. It is not attributed to the general compiler.

The adjacency target is the new order-six repair found in this session. Its
scalar lower boundary is the independently audited [SHIFT_PROOF.md](SHIFT_PROOF.md):
SCCs through order five cannot implement the bulk clock-uniform `+4` return
from this sibling. Finite exceptions at individual clocks are explicitly
outside that obstruction. The two finite targets here are recovered from the
audited `spectrum_acceptance.json` receipt, not asserted as new scalar counts.

## 2. Bulk adjacency return, with all old names retained

Let `t>=3`, and put `k=t-2>=1`. Change the frame, without changing any arc:

```text
new prefix: p0,...,p(k-1),
new suffix positions0,...,6: p(k),p(k+1),b0,b1,b2,b3,b4.
```

The reframed source suffix is `TT7` with port mask84, that is positions
`{2,4,6}`. Its actual signature is

```text
(h,p,f,g,j,e,s)=(1,128,84,21,12,1,0).
```

Change its first six vertices from `TT6` to the strong tournament

```text
core = (2,4,57,3,11,27).
```

All six still beat the last suffix vertex. The final seven-vertex suffix is

```text
B' = (66,68,121,67,75,91,0),     W'={4,6},
sigma(B',W')=(17,448,252,105,46,17,4).
```

Its two SCCs have orders6 and1; the order-six core is literally strong.
The old source vertex `b0`, now position2, is changed from a root port to a
nonport. The physical root remains `r`. Thus

```text
H(target_t)=-29+252*2^(t-2)=3(21*2^t-11)+4,
F_r(target_t)-H(target_t)=46.
```

The discrepancy schedule is the following ten pairs in suffix positions,
in the displayed order, then the root-port pair `(2,r)`:

```text
(0,2),(0,3),(0,4),(0,5),
(1,3),(1,4),(1,5),
(3,4),(3,5),(4,5),
(2,r).
```

At every completed adjacency stage, suffix vertex6 is a terminal singleton
and remains a root port. The remaining prefix is nonempty. The exact
terminal-SCC criterion from S182 therefore proves that **all twelve completed
adjacency states are strong for every `t>=3`**, including the initial state.
No claim about these states depends on the partial six-vertex core being strong.

## 3. Two finite patches and a sharp empty-prefix hostile

At `t=1`, keep `p0` as the one-vertex prefix and use the five `b` names, in
their existing order, for the suffix

```text
rows=(18,8,3,5,14), ports=19,
sigma=(13,104,55,58,26,8,9).
```

This is the inherited `V5_29_13` scalar witness. Its count is97, so it gives
`31 -> 97`. The final root defect is26. There are9 old-pair discrepancies.

At `t=2`, keep only `p0` as the prefix; use
`(p1,b0,b1,b2,b3,b4)` as the suffix, preserving those physical names:

```text
rows=(56,33,35,38,46,0), ports=33,
sigma=(11,200,113,37,14,11,1).
```

This is the inherited `V5_c_1 -> T1` witness. It gives `73 -> 223`, with final
root defect14. There are11 old-pair discrepancies. In each patch, reversing
the discrepancies in lexicographic order on the fixed physical vertex list
keeps every completed adjacency state strong; this is checked on the literal
finite graphs. Every helper/native intermediate is checked separately below.

The unpatched bulk recipe at `t=2` fails. It would leave an empty outer prefix.
The actual root-only target has

```text
H=17+252-46-4=219,
```

whereas `3*73+4=223`. Its literal count is independently verified. This is why
the bulk proof cannot simply substitute zero into the nonempty-prefix formula.
The finite patch supplies the missing actual graph; it does not erase this
failure boundary.

## 4. The inherited legal compiler instantiated on every discrepancy

For a current named tournament and selected pair `{a,b}`, construct a directed
Hamilton word on the other vertices by insertion immediately before the first
out-neighbor, or at the end if no such neighbor exists. Insert the pair block
before the first remaining vertex beaten by either endpoint, choosing its final
endpoint to beat that vertex. If there is none, append the block. All other
path arcs are directed in the current tournament. Choose the pair orientation
to agree with the word; the result is legal either in the current tournament
or its single-pair reversal.

For a legal traced pair `u->v` on a path `P`, use the inherited recipe:

- If `v` is not last in `P`, birth helper `z` beating only `u`, use guard
  `Pz`, then the guard obtained by inserting `z` between `u,v` in `P`.
- If `v` is last, birth `z` losing only to `v`, use guard `zP`, then the
  same internal insertion guard.

The first native complement makes the second guard legal. The two protected
edge sets, restricted to old vertices, differ only at `{u,v}`. Deleting the
helper therefore reverses exactly that old pair. If the initial pair word was
legal in the reversed tournament, traverse these two involutions backward.
The helper's birth cut is then the formerly final cut.

Both temporary cuts have either one or two in-neighbors, or one or two
out-neighbors. They are nontrivial in this family. All old names and every other
old arc survive each gadget. This constructs actual legal guards at every clock
in polynomial time; no Hamilton-path search is hidden. With a constant number
of discrepancies, the direct implementation takes `O(t^2)` adjacency tests.

The independent interpreter in `native_six.py` uses labeled directed adjacency
sets. It reconstructs each helper cut, checks every full guard word, performs
the off-path complement, and deletes only the named helper. The compiler's bit
matrix states and the interpreter agree after both native moves and deletion.
Repeated-name and reverse-directed illegal guards are explicitly rejected.

## 5. Why every native intermediate is strong for all bulk clocks

The general helper compiler alone does **not** preserve strongness. This
section supplies the additional family-specific argument, without extrapolating
from a few clock values.

At each of the eleven bulk discrepancy stages, the old graph has a prefix block
`P=TT_k`, followed by the same seven named suffix vertices, and the same root.
Every prefix vertex beats every suffix vertex; the root beats every prefix
vertex. The selected pair is always outside `P`. The outside adjacency depends
on the discrepancy stage, but not on `k`.

**Guard block lemma.** The direct guard algorithm keeps `P` as one consecutive
ordered block, and its word after replacing that block by one symbol, pair
orientation, append/prepend case, and both helper cuts are independent of
`k>=1`.

To prove this, first delete the selected outside endpoints just for constructing
the word. The prefix vertices are processed first, and append in transitive
order. Every subsequently processed suffix vertex loses to every prefix vertex,
so its insertion occurs after the entire prefix. If the root remains, it is
processed last and inserted immediately before the prefix, since it beats every
prefix vertex and no suffix vertex has been inserted before the prefix.

If both selected endpoints are suffix vertices, they lose to the whole prefix.
Their final block therefore goes before the root, among suffix vertices, or at
the end; it never splits the prefix. If the root is selected, it beats the first
prefix vertex, so the pair block goes immediately before the prefix. Every
choice depends only on the fixed outside arcs and these uniform comparisons.
At least five unselected suffix vertices remain after the prefix. Hence both
ends of the full pair word lie outside the prefix. Inserting the helper at the
selected outside pair, or appending/prepending it, also leaves the prefix block
intact. The forward helper cuts refer only to the selected endpoints, and the
reverse cuts only to those endpoints and the first/last word vertices. These
are all outside `P`. The helper cuts are therefore uniform on the prefix and
independent of `k`. This proves the lemma, including the reversed-gadget case.

Immediately after the first native complement, the induced prefix graph is

```text
R_k: i -> i+1, and j -> i whenever j>=i+2.
```

For `k>=3`, it contains the directed spanning cycle
`0 -> 1 -> ... -> k-1 -> 0`, so it is strong. The graph on the nine outside
vertices (seven suffix vertices, root, helper) is independent of `k`. Each
outside-to-prefix incidence is uniform except possibly at the first or last
prefix vertex: only the path's entry and exit edges escape complementation.
The uniform interior incidence and the two boundary exceptions are fixed by the
compressed word. For every `k>=3`, first, last, and interior classes are all
nonempty.

Contract the strong block `R_k` to one vertex, retaining an edge in each
direction that occurs. This is a directed quotient, not necessarily a
tournament. Its reachability graph is identical for all `k>=3`, and the full
graph is strong if and only if this quotient is strong. The eleven literal
`k=3` quotients are checked exactly and are all strong. Thus **every first-native
intermediate is strong for every `k>=3`**. The cases `k=1,2` are the separately
saved literal programs at `t=3,4` and are also strong.

Before the first complement and after the second, the restriction to old
vertices is one of the completed strong adjacency states. The helper has both
an old in-neighbor and an old out-neighbor by the explicit nontrivial one/two
exception cut. Adding it to a strong old graph therefore preserves strongness.
Deleting it gives the next already-proved strong adjacency state.

This proves all-intermediate strongness for every bulk clock. The two finite
patches are individually replayed and every one of their intermediate graphs
is also strong. Their combination proves the theorem for every `t>=1`.

## 6. Changing physical defect and the observer boundary

The root remains the same physical `r`. The exact source defect is12. At the
final target it is26 at `t=1`,14 at `t=2`, and46 at every bulk clock. The finite
receipt recomputes `H`, `F_r`, and `F_r-H` at all159 helper/native states and
all58 completed adjacency states in the five saved programs.

For example, the first gadget at `t=3` reverses the pair `(p1,b0)`. Its root
observer `(H,F_r,F_r-H)` is

```text
source:           (157,169,12)
helper birth:     (1575,1827,252)
first native:    (2769,5560,2791)
second native:   (1331,1441,110)
helper deletion: (213,249,36).
```

Thus neither the Hamilton count nor the defect is assumed invariant inside the
native program. Strongness is preserved by the separate proof in §5.

The coloured forest grammar can directly evaluate the eleven named edge-flip
macrosteps while retaining the root and seven suffix names and treating the
prefix anonymously. A Hamilton complement itself also changes internal prefix
arcs. It is **not** asserted to be one fixed-interface differential operator of
that grammar. The native compiler keeps the complete actual graph and guarded
path sidecar; a full-name forest could evaluate its individual edge changes if
the names had been retained from the start. No forgotten-name reactivation or
scalar defect reset is used.

## 7. Controls, reproducibility, and stopping boundary

The source is
[`native_six.py`](../../../../04-computation/multiport_grammar_20261103/native_six.py).
It contains the actual five programs at clocks1–5, a separate set-valued native
interpreter, literal endpoint-count recurrences, all intermediate matrices,
physical-root observers, and the compressed block-pattern controls. The three
bulk programs have identical compressed guard/cut/outside patterns. This finite
check supports the implementation; the structural guard lemma in §5 supplies
the unbounded quantifier.

The saved first five count pairs are

```text
31 -> 97, 73 -> 223, 157 -> 475, 325 -> 979, 661 -> 1987.
```

There are106 legally replayed native complements,53 temporary births/deletions,
159 strong helper/native states, and58 strong completed adjacency states.
No old name is removed. The first program uses18 complements; the other four
use22. The unpatched empty-prefix hostile is also recounted literally.

Reproduce with:

```text
python3 -B 04-computation/multiport_grammar_20261103/native_six.py
python3 -B -O 04-computation/multiport_grammar_20261103/native_six.py --output /tmp/native-six-optimized.json
```

Both modes pass **1,521 always-active gates** and reproduce the frozen
`native_six.json` receipt and491-byte stdout exactly. The
[independent audit](NATIVE_SIX_AUDIT.md) separately accepts the guard-block
and contraction arguments, replays every program and physical-root readout,
and checks unseen bulk sizes. Root promoted the status after that audit;
the producer source and receipt remain unchanged.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| `native_six.py` | 16,570 | `fea5f50b76d21b81b9d2312fbc6b9203ae3c465ea640c69e107acc67be824e96` |
| `native_six.json` | 144,516 | `5f239c3638cb88409a51b199fc95e2306f6da03e054f8acab93f4a3ff3e066aa` |
| `NATIVE_SIX.out` | 491 | `bc769537b6f86910e557beccfb76da91669fe0b7c950fd98d710e8bafdfe9185` |

This gives a constructive `44` map from the entire original `+4` family, with
the same old names and strong native intermediates. It does not give the
missing `24` map, arbitrary-depth two-branch closure, or global integer-spectrum
coverage. The separate [terminal preparation](PREPARED_FOUR.md) gives one
further `+4` from its actual bulk output, followed by an indefinite zero channel;
the next C2 digits2/4 are then blocked. Those return obligations remain
separate even though this formerly blocked branch now has actual guards.
