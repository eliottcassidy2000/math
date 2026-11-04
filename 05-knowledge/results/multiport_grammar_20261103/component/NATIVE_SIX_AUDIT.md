# Independent audit: the all-clock, all-strong native44 map

2026-11-03. **INDEPENDENT ANALYTIC / SOURCE / REPLAY AUDIT PASS.**
I read the full [NATIVE_SIX.md](NATIVE_SIX.md), its standalone source, and the
two cited inherited compiler proofs. Both independent normal/optimized runs
match the frozen receipt and stdout. No correction remains. Only this audit
sidecar was added; producer artifacts and previous frozen audits were preserved.

## 1. Physical statement and finite-clock boundary

The domain is exactly the original first-4 family, with fixed physical names
`p0,...,p(t-1),b0,...,b4,r`, count `21*2^t-11`, and root defect 12. The target
is on precisely those same old names, with count `63*2^t-29`. The bulk target
for t>=3 reblocks two old prefix vertices into the suffix; reblocking is an
identity of the literal graph, not a birth, deletion or physical relabelling.

The eleven bulk discrepancies are those independently audited in
[SIX_AUDIT.md](SIX_AUDIT.md). Their completed adjacency states stay strong
because the final suffix singleton remains marked. The final root defect is
46. Every old vertex and the same physical root survive.

The two finite patches are independently legitimate actual graphs on the same
names. At t=1 the five old b-names supply the listed five-vertex marked suffix,
giving `31->97` and defect 26. At t=2 the old p1 joins the five b-names in the
listed six-vertex suffix, giving `73->223` and defect 14. Their discrepancy
lists have lengths nine and eleven respectively. No scalar witness is treated
as a physical program without the actual named discrepancy compilation.

The empty-prefix hostile is retained: applying the bulk graph formula at t=2
gives actual count 219, not 223. The two patches therefore produce a
**piecewise all-clock algorithm**, not a single unchanged bulk suffix at the
two boundary clocks. The claimed resource counts follow: 18 complements at
t=1, 22 at every t>=2, one helper per discrepancy, no permanent births, no
old-name deletions, and peak order t+7.

## 2. Exact guard and helper compilation

The current proved supplier
[`recursive_seed_switching_20261103/transport/FOLLOWUP.md`](../../recursive_seed_switching_20261103/transport/FOLLOWUP.md)
constructs the pair-block guard by insertion. I recovered its actual theorem
and checked that the implementation agrees. Inserting a new vertex before
its first out-neighbour preserves a directed Hamilton word: every preceding
vertex dominates it, and it dominates the selected successor. For the selected
pair, all vertices preceding the first vertex beaten by either endpoint
dominate both endpoints. Choosing the pair's last endpoint to dominate that
successor makes every guard arc except the pair itself legal. Choosing the
pair orientation pays that final arc.

The two-move helper construction is the one in
[`adaptive_metagraph_reformulation_20261103/METAGRAPH.md`](../../adaptive_metagraph_reformulation_20261103/METAGRAPH.md),
§5. For a pair `u,v` traced consecutively in P:

* If v is not last, the initial helper beats only u. The first guard is P
  followed by the helper; the second inserts it between u and v. After both
  moves the helper beats exactly v and the old last vertex.
* If v is last, the initial helper loses only to v. The first guard begins
  with the helper; after both moves it loses exactly to u and the old first
  vertex.

These endpoint cuts are nontrivial for the stated domain. The first native
complement makes the two new links in the second guard directed. Restricted
to old pairs, the two protected-edge sets differ exactly at the chosen pair;
all other old pairs are changed either zero or two times. Deleting the helper
therefore reverses precisely that pair.

A native complement is involutive along its retained directed path. If the
pair word was legal in the reversed old tournament, swapping the helper
endpoints and reversing the two guarded moves is therefore legal, with the
formerly final cut used for the actual birth. The implementation retains the
inverse flag and the actual guard words; it does not pretend that the forward
birth cut is unchanged in this case. All helper names are fresh, and only that
helper is deleted.

The cited general compiler does not assert strongness. The new family-specific
argument below pays that additional predicate rather than importing it from
the inherited theorem.

## 3. Guard block invariance for every bulk size

Write k=t-2. At every completed bulk stage the first k vertices form TT_k,
each prefix vertex dominates every suffix vertex, and the root dominates every
prefix vertex. All selected discrepancy pairs are outside the prefix. Changes
to suffix/root arcs leave these uniform comparisons intact.

The insertion algorithm processes prefix vertices first. They append in their
transitive order. Every remaining suffix vertex loses to the entire prefix,
so it is inserted after the whole block. If the root is not selected for the
pair, it is processed last and inserts immediately before the prefix. If the
selected pair contains the root, its pair block instead inserts before the
prefix. A suffix-only pair never splits the prefix: neither endpoint beats
any prefix vertex. These choices depend only on the fixed outside arcs and
their uniform comparisons with the prefix.

At least five unselected suffix vertices remain after the prefix. Its first
and last word positions consequently have outside neighbours, and both ends
of the complete pair word are outside the prefix. Appending/prepending the
helper or inserting it into the selected outside pair leaves the block intact.
The explicit endpoint cuts above refer only to the selected pair and the old
first/last word vertices, all outside the prefix. This also pays the inverse
gadget's birth cut. Thus the compressed guards, pair orientation, helper case
and both endpoint cuts are independent of k>=1, and the helper's incidences
with the prefix are uniform.

This is a structural statement about the actual insertion loop, not an
inference from agreement at k=1,2,3. The frozen compressed-pattern checks
corroborate the implementation. Additional audit controls at k=4,7,17 also
agree and remain strong, but are not used to infer the unbounded quantifier.

## 4. All-intermediate strongness

After the first native complement, the induced prefix is exactly R_k: its
consecutive forward edges are retained and all other transitive edges are
reversed. For k>=3 the cycle

```
0 -> 1 -> ... -> k-1 -> 0
```

spans that prefix, proving it strong.

The nine outside vertices are the seven suffix vertices, root and helper.
Their induced graph is independent of k. An outside-to-prefix orientation is
uniform apart from possible first/last exceptions, since the prefix is one
consecutive guard block and only its entry and exit path edges escape the
complement. For k>=3 the first, last and interior classes are all nonempty.
The existence of an edge in either direction between an outside vertex and
the whole prefix is therefore independent of k.

Contracting R_k produces the same ten-vertex directed quotient for every
k>=3. Both directions between quotient vertices are allowed. Because the
contracted prefix is strong and other parts are singletons, the complete
graph is strong if and only if this quotient is strong. All eleven literal
k=3 quotients are checked and strongly connected. This proves first-native
strongness for every k>=3; the separately saved k=1,2 graphs supply the two
remaining cases where contracting the prefix would not be justified.

Before the first native move and after the second, the old restriction is one
of the already-proved strong completed adjacency states. The helper has an
old in-neighbour and an old out-neighbour by its nontrivial endpoint cut, so
adjoining it preserves strongness. Deleting it returns to the next strong old
state. This pays every native stage of the bulk algorithm. The t=1,2 patches
are finite literal programs; every one of their stages was independently
checked strong. Combining these facts gives the stated all-t theorem.

## 5. Independent implementation and observer replay

Both fresh runs of `native_six.py` reproduce **1,521 always-active gates**, the
144,516-byte receipt and 491-byte stdout exactly. The source's bit-matrix
compiler is checked by its separate labelled-set interpreter.

I additionally replayed the saved programs without importing the producer
compiler. Starting from each named source matrix, the audit reconstructed the
helper cut, checked that each guard is a permutation of the actual names and
is directed, complemented the whole graph, then restored the consecutive
guard edges. It compared both resulting matrices, deleted only the helper,
and checked that exactly the declared old pair was reversed. Independent
breadth-first reachability tested the actual saved matrices, not just their
stored strongness flags.

For every saved graph, a separate forward-append subset DP independently
computed h(S) and the number E_r(S) of paths ending at the fixed root r. It
then used

```
H=h(V),    F_r=sum_(S subset V) E_r(S)*h(V\S),    h(empty)=1.
```

This is a different recurrence direction from the producer's last-vertex
removal DP. The physical root is always at `names.index('r')`; the appended
helper is never substituted for it. All H, F_r and defect entries agree.

The full independent finite universe contains **106 legal native guards, 159
strong helper/native states and 58 strong completed adjacency states**. The
extra audit execution has 1,120 always-active checks, including the additional
unseen-prefix-size implementation controls. Its finite saved-program replay
does not call the producer compiler; only those separate larger-size controls
call the compiler under test. No count or defect invariance inside a native
gadget is assumed. The first t=3 gadget's displayed sequence
`(157,169,12)`, `(1575,1827,252)`, `(2769,5560,2791)`,
`(1331,1441,110)`, `(213,249,36)` is independently recovered.

The coloured endpoint-forest grammar handles the named arc-reversal macrosteps
with the declared active names. A whole Hamilton complement also changes
internal anonymous-prefix arcs; it is correctly not claimed to be one fixed
bounded-interface forest operator. The native program retains the full graph
and actual guard path. This distinction avoids hidden name resurrection or
collapsing the changing root defect to a phase.

## 6. Frozen pins, reproduction and scope

Independent producer-mode commands, from the worktree root:

```
python3 -B 04-computation/multiport_grammar_20261103/native_six.py --output /tmp/native_six_audit_normal.json
python3 -B -O 04-computation/multiport_grammar_20261103/native_six.py --output /tmp/native_six_audit_optimized.json
```

| Artifact | SHA256 |
|---|---|
| `NATIVE_SIX.md`, pre-promotion | `aaebe45afb4c75d46e503e4104a942776a6a7291e1be576301d01aaa6867d191` |
| `native_six.py` | `fea5f50b76d21b81b9d2312fbc6b9203ae3c465ea640c69e107acc67be824e96` |
| `native_six.json` | `5f239c3638cb88409a51b199fc95e2306f6da03e054f8acab93f4a3ff3e066aa` |
| Frozen and independent stdout | `bc769537b6f86910e557beccfb76da91669fe0b7c950fd98d710e8bafdfe9185` |

**PASS.** The accepted result is an explicit, old-name-preserving, all-strong
native +4 map on the entire original first-4 family, with two finite patches
and a uniform bulk mechanism. Its outputs are not proved to re-enter that
source family. It supplies no missing24 map, reusable two-branch physical
class or global integer-spectrum classification. Those statements are neither
assumed nor inferred. Parent owns promotion and Git integration.
