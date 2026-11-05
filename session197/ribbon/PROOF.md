# Bounded-backward-span tournaments: an independent frontier certificate

Session 197, 2026-11-05. This packet is developed independently from the definition below. It does not rely on the previously reported ribbon transfers or their strong-component criterion. The contextual Session 196 guide was read, but no mathematical claim here is imported from it.

## 1. Definition and scope

A span-`w` tournament is a tournament on `1,...,n` in which `i -> j` whenever `j-i>w`. The pairs at distance at most `w` can have either orientation. This is a property of the supplied order; an unlabelled tournament belongs to the class if it admits such an order. Write `H(T)` for its number of directed Hamiltonian paths, with distinct vertex sequences counted separately. All spectra below allow arbitrary order and do **not** impose strong connectivity.

A word `m_1 ... m_n` over `0,...,2^w-1` specifies a tournament. When vertex `j` is added, bit `i` of `m_j` says that `j` points to position `i` of the retained window `(j-w,...,j-1)`, with positions numbered `0,...,w-1`. Positions below `1` are padding and ignored. All remaining old vertices point to `j`. The word convention is exactly `graph` in `transfer.py`.

## 2. Finite frontier transfer

A selected set of path edges on the processed vertices is a disjoint union of directed paths, including singleton paths. Represent each component by `(start,end)`. Retain the last `w` vertex labels as `0,...,w-1` and replace every older endpoint label by `-1`. Keep a state only if at most one component start is `-1`.

This last restriction is necessary for every partial edge set that extends to a Hamiltonian path: a forgotten start has no selected incoming edge yet, and no future vertex can point to it, since their distance exceeds `w`. Hence it must be the eventual global start. Two forgotten starts can never be repaired.

The ending labels `-1` can appear on several components. They are not pooled: each is attached to its component's separately retained start. Thus no multiplicity is lost. Other vertices of a component need not be stored.

On adding vertex `v`, choose:

1. no incoming edge, or an edge from the end of one existing component to `v`;
2. no outgoing edge, or an edge from `v` to the start of one existing component.

The incoming choice is legal when its endpoint is forgotten or its mask bit is zero. The outgoing choice is legal only for a retained start whose mask bit is one. If both chosen components coincide, reject the choice because it closes a directed cycle. Otherwise merge the chosen components through `v`, or create its singleton if neither was chosen. Shift retained labels one position down; the new vertex obtains label `w-1`. Forget label zero and reject states with two forgotten starts. Different choices leading to the same state are counted with their multiplicity.

**Transfer theorem.** Starting from the empty state with weight one, the sum of the weights of one-component states after reading a nonempty word is exactly `H(T)`.

**Proof.** At each step these are precisely all possible choices of the incident selected edges of the new vertex in an acyclic graph of maximum indegree and outdegree one. Every such graph is a disjoint union of paths. Distinct sequences of edge choices yield distinct selected edge sets, so the recurrence neither omits multiplicity nor duplicates an edge set. A Hamiltonian path has a unique restriction to each prefix, and none of those restrictions is removed by the forgotten-start condition proved above. Conversely, each terminal one-component state represents a spanning directed path. Future compatibility depends only on the recorded endpoint labels: all forgotten ends point to all future vertices, no future vertex points to a forgotten start, and retained orientations are read from the word. Thus the aggregation is exact. QED.

There are finitely many states: at most `w` components with retained starts and one with forgotten start, and every endpoint belongs to the finite label set. Generating all reachable states gives 7, 27, 119 and 571 states for widths 1, 2, 3 and 4 respectively. These are counts of this explicit construction, not minimal-state claims.

## 3. All-order bounded-output closure

Fix `M`. Replace every coordinate by `min(M+1,x)` after each transfer. Nonnegative integer linear transfers commute with this coordinatewise saturation in the truncated semiring. In particular the saturated terminal sum determines `H(T)` exactly whenever `H(T)<=M`, and otherwise reports a number greater than `M`.

We may also discard every word prefix already having more than `M` Hamiltonian paths.

**Monotonicity lemma.** If `U` is an induced subtournament of `T`, then `H(U)<=H(T)`.

**Proof.** It suffices to insert a single vertex `v` into each Hamiltonian path `p_1...p_k` of `U`. If `v -> p_1`, put it first. Otherwise `p_1 -> v`. If `p_k -> v`, append it; in the remaining case take the first index at which the relation changes from `p_i -> v` to `v -> p_{i+1}`, and insert there. Thus each old path has an extension. Choose one deterministically. Deleting `v` recovers the old path, so the extensions are distinct. QED.

Run breadth-first search on saturated vectors, beginning with the empty-state unit vector. Extend by every mask, discard outputs above `M`, and deduplicate vectors. The search space is finite, bounded by `(M+2)^d` where `d` is the number of frontier states. Completion of this search is an **arbitrary-order** certificate, not a census up to a chosen vertex order: any longer word follows already visited vector transitions, and every discarded prefix is excluded by monotonicity. `explore.py` implements exactly this algorithm.

## 4. Verified finite-output conclusions

The dependency-free `verify.py` completes these searches:

- Width 2, `M=47`: 1,488 saturated vectors; the positive odd gaps are exactly `{7,21,47}`.
- Width 2, `M=1001`: 192,274 vectors; the positive odd gaps are exactly `{7,21,47}`.
- Width 3, `M=101`: 100,168 vectors; the positive odd gaps are exactly `{7,21}`.

No positive even value appears. The assertions about these bounded intervals come from the completed exact searches; no claim of global coverage is made.

In particular, **47 has minimum backward span exactly 3**. Width 2 excludes it at every order. At width 3 the word `0001321` gives a strongly connected order-7 witness. Its reversed pairs, written as `[smaller,larger]`, are

`[1,4], [2,5], [3,5], [4,6], [4,7]`.

All other edges point from the smaller to the larger vertex. Held–Karp subset dynamic programming independently gives 47, and a reachability computation verifies strong connectivity. This witness is not asserted to be order-minimal among all tournaments.

The negative conclusions also hold for the strong subclass. The positive coverage conclusions do not assert that every listed value has a strong witness. No interval-based strong-component test is used.

## 5. Exact all-order classification at width 1

Let `R_n` be the tournament with every adjacent pair reversed and all other edges forward. Put `h_n=H(R_n)` for `n>=1`, and `h_0=1` as an empty-product convention. Then

`h_0=h_1=h_2=1`, and `h_n=h_(n-1)+h_(n-2)+h_(n-3)` for `n>=3`.

**Proof.** In a Hamiltonian path of `R_n`, every descent must decrease the vertex label by one. Its maximal decreasing runs are therefore disjoint integer intervals, and they appear in their natural increasing interval order: a transition to a lower interval would be a forbidden larger descent, unless it continues the same run. A path is consequently a composition of `[n]` into consecutive intervals, each read in decreasing order. A transition between consecutive intervals of sizes `a,b` has forward length `a+b-1`, and is legal exactly when the two sizes are not both one. Thus `h_n` counts compositions of `n` with no adjacent parts equal to one. If `B=x^2/(1-x)` marks a part of size at least two, this composition generating function is

`(1+x)/(1-B(1+x)) = (1-x^2)/(1-x-x^2-x^3)`.

The recurrence follows. QED.

Every forward adjacent edge in a width-1 tournament separates an ordered join: all edges across that cut point forward. Splitting at these cuts gives blocks `R_n`. Hamiltonian paths in an ordered join visit its factors consecutively and in factor order, so their counts multiply. Conversely, ordered joins of the blocks remain width 1. Therefore the **complete width-1 spectrum** is the multiplicative monoid generated by

`1, 3, 5, 9, 17, 31, 57, 105, 193, 355, ...`.

The width-1 strong tournaments are exactly `R_n` for `n>=3`, together with the one-vertex tournament, up to the specified order. For `n>=3`, the descending adjacency path together with the forward edge `1 -> n` is a spanning directed cycle, proving strongness. `R_2` is not strong; a forward-cut decomposition into more than one block cannot be strong.

This classification is all-order and symbolic. It does not extend to width 2 by merely replacing the recurrence: backward length-two edges can bridge a forward adjacent cut.

## 6. Reproduction and limits

Run from the repository root:

`PYTHONPATH=session197/ribbon python session197/ribbon/verify.py`

Python 3 only; no third-party packages. Assertions must be enabled. The string witness encoding in `explore.py` is used only at widths 1–3; widths 4 and above require a delimiter or list encoding because mask values may have multiple digits. The general frontier transfer itself takes integer mask lists and has no such encoding limitation. The verifier compares the frontier count to independent Held–Karp counts on all 64 width-1 order-7 words, all 2,048 width-2 order-7 words, all 4,096 width-3 order-6 words, and 200 reproducibly randomized examples up to order 10. It reruns the complete bounded closures and rechecks every stored representative word with the transfer.

`VALIDATION.json` records the outcomes. `closure-2-47.json`, `closure-2-1001.json`, and `closure-3-101.json` give shortest-length representative words for the respective attained outputs. BFS witnesses are shortest among this word class because every extension adds one vertex.

Remaining question: do all positive odd integers outside `{7,21,47}` occur at width 2, or all outside `{7,21}` at width 3? These global assertions are not proved here. The bounded searches and the width-1 classification are not a substitute for an all-target construction.
