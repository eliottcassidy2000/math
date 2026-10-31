# Independent audit: cycle-realization and bounded packing

## Result first

**PASS, with one bibliographic correction.** The three mathematical results in
`RESULTS.md` are valid with their stated restricted-family
scope. In particular, **K14 cannot be the literal odd-cycle conflict graph of
any tournament**, even though 29 is an ordinary attainable Hamiltonian-path
count. The finite forbidden-subsequence classification is all-order, and the
entropy theorem concerns word presentations and isomorphism classes only.

The OCF reference in the source packet has been corrected to read **Irving–Omar, Corollary 19**, not
Corollary 20, for the linked v1. This was checked in both the
[v1 PDF, page 13](https://arxiv.org/pdf/2412.10572v1) and
[v1 HTML, section 3.4](https://arxiv.org/html/2412.10572v1#S3.SS4).
The formula and its attribution are correct; only its number needs repair.
The inherited `THM-002-ocf.md` repeats the same numbering error. The original audit left the source packet unchanged; this publication copy applies
the citation repair, recorded in `CITATION_CORRECTION.md`.

**No proof of global spectrum completeness, square-sum completeness, or a
density-zero theorem is supplied or implied. No publication was performed.**

## Finite-exact verification

All six entries in the source `SHA256SUMS` matched before replay. All three
source programs were rerun normally and with `python -O`; their six receipts
match the original corresponding receipts and each other. Published receipts
are in `receipts/`; filenames identify the program and optimization mode.

The new `independent_check.py` imports no executable source-packet code. It
constructs Boolean adjacency matrices, detects triangles by induced outdegree
one, counts Hamiltonian paths by literal backtracking without endpoint DP,
and recognizes subsequences using an iterator rather than prefix states. It
independently regenerated:

| Quantity | Verified result |
|---|---:|
| Triangle-core presentations, spine lengths 0 through 5 | 37,449 |
| Accepted presentations | 3,624 |
| Accepted by total order 3,4,5,6,7,8 | 1, 8, 64, 272, 869, 2,410 |
| Minimal forbidden length-3 / length-4 words | 240 / 35 |
| Distinct proper prefixes including empty | 84 |
| Maximal freely repeatable alphabets | the same nine binary sets |
| Independent explicit checks per execution | 38,207 |

The exact spectrum through 73 is independently recovered, including checking
that 7 and 21 do not accidentally occur in this restricted class:

```
1,3,5,9,11,13,15,17,19,23,25,27,31,33,37,39,
45,47,49,51,53,61,63,65,69,73.
```

Every accepted presentation had its Hamiltonian paths independently counted.
Literal cycles were separately enumerated by support and cyclic permutations
for 72 `(order,H)` representatives; the source replay checks literal-cycle
OCF on every accepted presentation. All nine ordinary witness arrays were
checked as actual bounded bitmask tournament adjacencies, and their exact key
set, paths, cycles, disjoint pairs, and OCF values were independently checked.

A hostile support-collapse control is the legal word `(1,4)` on five vertices:
it has H=13 and six literal odd cycles but only five distinct supports.
Replacing literal cycles by supports would falsely give 11. This independently
tests the distinction needed by the conflict-graph statement.

## Theorem 1: why the finite computation is all-order

1. A directed cycle's induced tournament contains a directed triangle. Replacing
   each member of a vertex-disjoint cycle family by such a triangle preserves
   disjointness. Thus the maximum packing sizes of cycles, odd cycles, and
   triangles coincide.
2. Tournament SCCs are linearly ordered. Every non-singleton SCC has a triangle;
   two such SCCs would give packing at least two. With packing at most one there
   is consequently at most one nontrivial SCC. Directed Hamiltonian paths must
   traverse the SCCs in order, so singleton SCC removal preserves H. The
   completely transitive case contributes H=1.
3. [Busch, Theorem 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf)
   gives a strong-tournament lower bound at least `5^((n-1)/3)`. The stronger
   residue-class constants in the paper are not needed. Since
   `5^8=390625 > 73^3=389017`, a strong representative with H at most 73 has
   order at most eight. There is no unstated asymptotic or census assumption.
4. Choose a triangle. Its complement has no triangle and is therefore
   transitive. Ordering that complement and recording its three core-edge bits
   gives exactly an eight-letter word. This is a surjection from presentations
   onto the relevant tournament classes, not an injective canonical labeling.
   Exhausting spine lengths zero through five covers all required strong
   representatives. Enumerating additional nonstrong presentations is harmless.

These points justify the claimed union over all tournament orders. They do not
assert that each listed value occurs at every individual order.

The nine values `29,35,41,43,55,57,59,67,71` have independently checked ordinary
witnesses but do not occur in the restricted spectrum. Their forced disjoint
triangle pair is therefore a genuine all-order necessary condition.

## Theorem 2: forbidden subsequences, phases, and entropy

Every triangle meets the core, since the spine is transitive. Two disjoint
triangles use at least three spine vertices (six total vertices minus the
three-vertex core) and at most four (each triangle uses a core vertex).
Their induced word is therefore bad at length three or four. Deletion to a
minimal word yields exactly the enumerated basis, and every such subsequence
retains the disjoint triangles when embedded in a longer word. This proves
both directions at arbitrary length; the length-five check is a cross-check,
not the logical source of the all-order claim.

All immediate deletions of every one of the 275 patterns were independently
verified good. Heredity then gives full proper-subsequence minimality.

For a freely repeatable alphabet, variation in two distinct core coordinates
allows a `1,0` pair of distinct spine occurrences for each coordinate. Concatenate
the chosen pairs, using separate occurrences: each pair together with its own
core vertex forms a triangle, and the two triangles are disjoint. Thus at most
one coordinate varies. Of the twelve cube edges, exactly three have the other
two fixed bits oriented so that one spine occurrence forms a triangle with
those two core vertices. Such a triangle can be made disjoint from a variable
coordinate triangle. The other nine edges are safe because all their triangles
contain the variable core vertex. This validates the stated binary alphabet
classification in both directions.

The prefix state is the set of realized proper prefixes. It only gains
coordinates; the empty prefix is present initially. Hence at most 83 transitions
are strict. At a fixed legal reachable state, any string of self-loop letters
is legal after the reaching prefix. Subsequence heredity implies its loop
alphabet is freely repeatable even when considered alone, so it has at most
two letters. Consequently there are at most **84 free phases plus 83 individual
strict-transition letters**, not 84 unrestricted contiguous binary blocks.

For specified strict positions and letters the automaton evolves
deterministically, with at most two possibilities at every other position. Thus
`L_m <= sum_j binom(m,j) 8^j 2^(m-j)` with `j <= min(83,m)`, exactly the stated
`2^m sum_j binom(m,j)4^j`. No independent choice of phase alphabets is being
assumed. The bound is valid but not claimed sharp.

## Theorem 3: general bounded-packing scope

For fixed positive k, a maximum packing of r at most k triangles supplies a
core of size 3r. Its complement is transitive; otherwise the packing could be
extended. There are finitely many possible labeled cores, independently of n.
For each such core, all cross edges are encoded by an alphabet of size
`q=2^(3r) <= 2^(3k)`.

Any forbidden packing of k+1 triangles uses at most 2k+2 spine vertices, because
each triangle meets the core. Keeping the entire fixed core and just those
spine vertices witnesses the violation. The legal language is therefore
characterized by finitely many forbidden subsequences of lengths at most
2k+2. Its nonempty proper-prefix coordinates number at most
`P_k=sum_(i=1)^(2k+1) (2^(3k))^i`. Shorter forbidden words cause no problem.

A recurrent alphabet varying in k+1 core coordinates allows k+1 separate `1,0`
spine pairs and hence k+1 disjoint triangles. It therefore varies in at most k
coordinates, so it has at most `2^k` letters. The same strict-position argument
gives `L_m <= 2^(km) sum_j binom(m,j)2^(2kj)` for `j <= min(P_k,m)`.

Summing over finitely many cores and r gives
`O(n^P_k 2^(kn))` **isomorphism classes**: each class has at least one such
presentation. No binomial factor selecting labeled core positions is required
for this unlabeled upper bound. The r=0 case is one transitive class per order.
The constants depend on k. This argument is not a bound on labeled tournaments,
not an injective representation theorem, and not an H-value density bound.

## K14 and the older poison-factor argument

The conflict graph on literal odd cycles is complete iff no two cycles are
vertex-disjoint. A realization of K14 would therefore have packing one and
OCF value `1+2*14=29`, contradicting the verified all-order restricted spectrum.
Meanwhile `I(K14,-1/3)=-11/3`, so `1+3x` is not a factor of its polynomial.
The six-vertex H29 witness has `I=1+10x+2x^2`; its two disjoint triangle pairs
are indispensable to attaining 29. All these implications are valid.

The historical output files corroborate the packet's reconciliation:
`k3_poison_complete.out` is empty; `k3_poison_fast.out` reports 1,099 graphs
through order five and then makes an unsupported universal factor claim.
K10, with polynomial `1+10x` and value 21 at x=2, refutes that universal claim.
The old H63 census records complete K31 conflict graphs, consistent with the
current independently recovered packing-one realization of 63. Historical
script sources are not present in the selected source snapshot;
the output statements were read directly, not represented as freshly replayed.

## Reproduction and integrity

From any directory:

```sh
python /path/to/package/independent_check.py
python -O /path/to/package/independent_check.py
```

The audit's normal and optimized JSON receipts are byte-identical. The three
source replays and their optimized counterparts are also included in
`receipts/`. `PROVENANCE.json` freezes the source identities and original audit evidence;
`receipts/verification.json` records portable replay equalities. `SHA256SUMS`
freezes this publication bundle. All validation uses explicit exceptions rather than
optimization-sensitive assertions.

Status: mathematical claims **VERIFIED / FINITE-EXACT plus audited all-order
arguments**; bibliographic repair **APPLIED IN PUBLICATION COPY**; unrestricted completeness
**OPEN / NOT ESTABLISHED HERE**. No action was denied or left uncertain.
