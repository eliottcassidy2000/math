# Marked completion plans consume proofs without promoting assumptions

**PROVED:** the selected-plan grounding, substitution, cut, and rank theorems
below, including strict first-hit certificate export. **FINITE-EXACT:** the
declared 64-plan universe and the named receipt controls. **CONDITIONAL:** an
open child port becomes a completed proof only when an authenticated supplied
proof or a finite grounded replacement is attached. No universal Collatz
completion or new all-height guard coverage is asserted.

## 1. Inheritance and the new interface

Finite authenticated common-future network connectivity is inherited from
[local/global grounding](paper_local_global_grounding_20261007.md), not a new
theorem of this package. Its root path, unit flow, and defect-cut equivalence
already separate local arithmetic from grounded evidence. The
[factorion completion model](factorion_completion_defects_20261007b.md)
distinguishes the least ROOT solution from an artificially completed map;
its actual nonroot cycles cannot be erased by an assumed completion edge.
[Poisson patching](collatz_poisson_patch_20261005.md) retains the boundary
forcing lost by naive local gluing. Finally,
[partitioned completion](collatz_partitioned_completion_20261004.md) compares
the final dependency with the original source and preserves actual guards.

The new operational interface is a **marked, source-labelled proof plan**.
An unverified terminal is a `Hole`, not an edge to ROOT. Replacing it consumes
actual certificates and preserves all other ports. The replacement may expose
a dependency cycle, which remains visible. Already grounded exports survive
every such replacement unchanged. This is an auditable way to use proposed
completion diagrams before all their terminal obligations have been solved.

The source is a finite collection of authenticated implications plus explicit
open labels. The target is either an actual first-hit ROOT word or an explicit
selected-plan defect. The map is finite dependency tracing followed by suffix
splicing. It preserves each original integer, ordered valuation words, and
the child identity at every interface. It intentionally does not preserve
the order in which independent holes were supplied. A compressed common-future
identity loses its proof dependencies; the retained plan is that sidecar.

The live concepts are ROOT seeds, actual forward prefixes, paid common-future
receipts, marked terminal ports, dependency rank, and signed route clocks.
The cheapest hostile is the actual `3 -> 5 -> 3` proof cycle in section 4.

## 2. Typed ports and actual transport

Let
\[
 U(n)=(3n+1)/2^{v_2(3n+1)}
\]
on positive odd integers. An actual valuation word is verified at its exact
source and may stop at 1 but may not execute a letter after reaching 1.
All numeric packet fields are exact integers; booleans and floating aliases
are rejected. The streaming verifier retains only the current integer, not
the entire orbit. Its running time and bit cost depend on the supplied word
and intermediate heights.

A selected plan has exactly one node at every listed source:

* `Root(n,w)` contains a checked first-hit word from `n` to 1.
* `Hole(n,label)` asks for a proof of this exact `n`, with no implicit edge.
* `Transport(n,h,a,b,z)` contains checked paths `U_a(n)=U_b(h)=z>1`.

Each transport's child port is explicitly listed. ROOT itself is the literal
`Root(1,())`. A forward transport has `b=()` and `h=z`; a paid transport also
requires `h<n`. The general common-future type does not claim that inequality.
If a supplied paid receipt already ends at 1, the adapter returns the direct
`Root(n,a)` rather than hiding its proved terminal behind a dependency.

Given `Root(h,w)`, determinism and the first-hit boundary force `w` to begin
with `b`. Hence
\[
 a\,w[|b|:]
\tag{1}
\]
is the actual first-hit ROOT word for `n`. The consumer checks the child
identity and this prefix before export. It does not accept a word for another
integer that merely has the same rank, residue, or apparent endpoint.

Two consecutive transports can be compressed while retaining an open terminal.
For paths `n -(a,b)-> h -(c,d)-> q`, the two words at `h` are comparable
prefixes. If `|b|<=|c|`, use `(a c[|b|:],d)`; otherwise use
`(a,d b[|c|:])`. This is the exact word version of the inherited clock rule.
A vacuous composite at the same source is an identity implication, never a
ROOT proof. `conditional_link` returns both the compressed implication and
the still-open `Hole`.

## 3. Grounding, ranks, and persistent defects

For a finite selected plan, tracing dependencies from a supplied source has
exactly three outcomes:

1. **GROUNDED:** a checked `Root` is reached;
2. **OPEN:** a labelled `Hole` is reached;
3. **CYCLE:** a previously visited port is reached.

It halts after at most the number of listed ports. This is completeness for
the **selected plan**, not for every mathematical proof or every alternative
rule bank. An unselected valid alternative can change the result.

### Rank theorem

The whole finite selected plan is grounded if and only if it has no holes
and admits a natural-number rank strictly decreasing from each transport's
source to its child. The forward direction assigns each node its distance
in selected dependencies to a `Root`. The reverse direction follows strict
descent and the fact that a finite terminal other than a `Root` is forbidden.
Induction using (1) then exports actual first-hit words. Ranks at supplied
ROOT nodes need not be zero in a verifying witness; the least proof-height
rank returned by `ranks` assigns zero there.

This rank is a proof-dependency height, not the number of actual odd Collatz
steps. On a grounded path ending at a supplied word `w`, the exact latter is
\[
 \tau(n)=|w|+\sum_{\text{selected edges}}(|a|-|b|).
\tag{2}
\]
Individual terms may be negative. `rooted_clock` requires a grounded path;
a formally consistent open or cyclic clock is not a ROOT deadline.

For partial plans, the complement of the least grounded set is a checkable
defect cut: it contains no `Root`, and every selected transport from a member
stays in the cut. Its endpoints may be holes or cycles. Such a cut proves
non-grounding in this plan; it does not prove nonconvergence of its integers.
Conversely a closed cut with no ROOT seed cannot contain a grounded node.

### Same-port substitution and monotonicity

`refine` replaces exactly one declared hole with a finite set of typed nodes
containing that same source. Existing other ports must be retained identically;
conflicting redefinitions are rejected. New nodes may refer to old ones, so
the resulting plan is traced again rather than assumed complete.

Every previously grounded path avoids the replaced hole. Therefore its nodes
and supplied ROOT word remain unchanged, and its exported word stays exactly
the same. The grounded set can only grow, although an old open path may become
a cycle. When the replacement is itself grounded without passing through the
hole again, all its ancestors become grounded by (1). This is a concrete
finite certificate substitution theorem, not a promise that every hole has
such a replacement.

`supply_root` also permits an independently checked ROOT word at an existing
cyclic or transport port. It selects that literal anchor in a new immutable
plan; the original plan, including the previous receipt, remains available.
This can break an exposed cycle. It never overwrites a port merely on the
strength of an assumed home predicate. Previously grounded sources remain
grounded, and uniqueness of the actual first-hit word preserves their exports.

At the schema level, an independently proved rule returning, for **every**
positive odd `n>1`, either an actual ROOT word or an authenticated dependency
`h<n` would give completion by strong induction. The universal production and
guard premise is essential. A successful finite table does not supply it.
More general schemas may use any independently proved well-founded rank.

## 4. The locally paid cyclic hostile and its repair

Start with the actual forward observation `U_1(3)=5` and an open ROOT port
at 5. One might attempt to discharge it by the locally paid receipt
\[
 U_{()}(5)=U_{(1)}(3)=5,\qquad 3<5.
\]
This arithmetic is valid, but it changes the dependencies to `3 -> 5 -> 3`.
There is no grounded ROOT seed in that component. Its signed clock increments
are `+1,-1`, so clock cancellation does not repair the absent anchor.
No strictly decreasing dependency rank exists. Both integers do in fact
reach ROOT; the failure is the proposed proof's circularity.

Supplying the actual word `(4)` for 5 repairs the original open plan. The
consumer exports `(1,4)` for 3. Supplying a word for a different source,
silently overwriting another port, or padding ROOT with the formal letter 2
is rejected. The same test would detect an artificially completed factorion
cycle: the changed edge must remain marked until authenticated.
The same checked anchor at 5 also repairs the cyclic plan itself through
`supply_root`; its immutable input continues to record the earlier CYCLE.

### Pure backward implications cannot create a ROOT seed

For `n>1`, a pure predecessor rule has the form `U_b(h)=n`, with `b` nonempty,
and treats `home(h)` as its premise. Its child cannot be 1: a strict first-hit
word cannot execute even its first letter from 1. Thus a finite bank consisting
solely of such predecessor implications cannot ground a nonroot source from
the singleton seed 1. Independent nonroot ROOT seeds or genuine forward
terminal information are required. This generalizes the same boundary in the
native `G1/G5/G17` bank; it says nothing against their usefulness for routing
to an already supplied seed. In particular odd targets divisible by 3 have
no nonempty odd Collatz predecessor at all.

## 5. Actual native-head completion and finite controls

The existing [two-anchor head decoder](collatz_twoanchor_head_decoder_20261007b.md)
exports an actual receipt for source 63829, child 2365, and join 281.
Installing this receipt with `Hole(2365,...)` remains OPEN. The separately
supplied child word

```text
3,1,1,3,3,2,1,3,1,1,3,4,1,3,1,2,3,4
```

is checked at 2365. Filling that exact port exports the same 15-letter
first-hit source word as the inherited head compiler. This tests the new
consumer; it is not a new native guard family or new source coverage.

### A newly supplied J2 child, with no orbit search in the consumer

The companion
[two-twos grounding packet](collatz_two_twos_grounding_20261007c.md)
supplies a source-specific, independently authenticated word for
`M_6125=2^6125-1`. It has 28,762 odd steps and total valuation 51,712.
The inherited two-anchor receipt at `M_6129` has the strictly smaller child
`M_6125`. The marked consumer first installs that receipt with an OPEN port
at the child. It then consumes `root_word(6125)` and exports the actual
`M_6129` word, of rank 28,762 and cost 51,716. A separate companion splice
agrees letter for letter, and both are replayed strictly with constant state
memory. The original immutable plan still reports OPEN.

These are two specific enormous integers, not an exponent progression proved
home by extrapolation. The companion owns discovery and the frozen data;
this program calls its verifier and receipt constructor, never `discover`.
The general receipt family remains conditional at every unsupplied child.

The independent finite selected-plan universe is all 64 choices on ports
`{3,5,7}`. Each port chooses a hole, its supplied ROOT word, or one of two
checked common-future implications through 5; literal ROOT is also included.
Across the 256 port queries, the outcomes are 142 GROUNDED, 78 OPEN, and
36 CYCLE. Every grounded export is replayed, every conditional open transport
is independently discharged using its supplied leaf control, and every
grounded plan's least rank and every defect cut are checked. Hole fillings
also verify monotonicity and byte-for-byte persistence of earlier words.
Malformed source types, wrong child identities, bad ROOT boundaries, rank
aliases, unlisted ports, and conflicting replacements are hostile controls.

Reproduction from the repository root:

```text
python -B -X utf8 04-computation/experiments/collatz_marked_completion_20261007c.py
python -B -O -X utf8 04-computation/experiments/collatz_marked_completion_20261007c.py
```

The program makes no production orbit search and contains no function that
converts an assumption of home-ness into a certificate. Its only anchors are
supplied words checked at their immutable source integers.
