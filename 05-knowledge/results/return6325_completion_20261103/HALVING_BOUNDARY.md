# Pair orientations have a sharp attachment-depth budget

**PROVED, scoped all-size statements; same-agent proof review.** Independent
finite controls are in `halving_response.py` and its JSON. This refines the
two-boundary-slots mechanism of
[BOUNDARY_BUDGET B1/B4](../prime_lane_ports_five_groups_20261102/BOUNDARY_BUDGET.md)
using signed cancellation. The older reversal-cancellation idea is recorded
in [THM-068, position-character-decomposition, part ii](../../../01-canon/theorems/THM-068-position-character-decomposition.md).
Neither precursor is presented as a new result here.

**Concurrent integration:** incoming [TM3](../digit_memory_all_bases_20261103/TOURNAMENT.md)
independently proves the general equal-profile-module degree bound, of which
HB1 is a special case. [INCOMING](INCOMING.md) credits the overlap and extends
the shared mechanism to specified new-new edge differences. The full rank
classification and recurrences are in [HB5/HB6](RANK_DEPTH.md).

## HB1. Births bound the degree of hidden pair interactions

Let a tournament have N disjoint two-vertex modules {a_i,b_i}. A module
has the same relation from either member to each vertex outside that pair.
Independently reverse their internal arcs, obtaining T_x, x in {0,1}^N.
Attach k named vertices by fixed cuts, with fixed arcs between new vertices,
independent of x. All other old arcs stay fixed. Then the Hamilton-path
count of the extension, as a multilinear Boolean function of x, has degree
at most 2k. The same bound holds for any fixed ordered path-cover count with
endpoint requirements naming only vertices outside the variable pairs.

**Proof.** Take a mixed difference on D, with |D|>2k, fixing other bits.
Expand the count as a sum over paths and orientations. If a path does not
use a selected pair's internal arc, summing that orientation cancels it.
Every surviving selected pair is therefore a consecutive two-vertex block.
At most 2k such blocks can be adjacent to new vertices. Choose the first
selected block with no new neighbour, in a fixed ordering of the pairs.
Swap its two vertices in the path and reverse its orientation bit. Relations
to neighbouring old vertices are unchanged by the module property; any
path endpoint remains permitted because the endpoint guards are external.
This is a sign-reversing involution. The choice of untouched block is
unchanged, so there is no unmatched term. For an ordered path cover perform
the swap inside its component; empty bins cause no change. QED.

This is a degree bound for this operation menu. Complements changing old
external arcs can destroy the module hypothesis. Naming an endpoint inside
a variable pair also escapes the stated guard scope. Neither operation is
silently included.

## HB2. An edge between new vertices spends two slots

The difference obtained by reversing one specified new-new edge has degree
at most 2k-2 in the old pair bits. Only paths using that edge contribute to
its difference. Its two incidences occupy two of the 2k slots that could
touch old pair blocks. Apply the same involution to any larger selected set.
The argument also applies to ordered covers with the guards of HB1.

In particular, reversing the edge between two new readers cannot change
a degree-three interaction. The actual C3[TT2,TT2,TT2] experiment has mixed
response -4 for either orientation of that new edge. This is a finite
instance of HB2, not a numerical coincidence.

## HB3. Positive laws with growing exact memory

Put uniform probability on even-parity x, or on odd-parity x. For N>2k,
all expected responses allowed by HB1 coincide for these two positive laws.
Their difference is exactly the top Walsh coefficient, which HB1 kills.
These are distributions on actual tournaments. No negative source weights
or nonrealizable port tensors are needed for the indistinguishability.

This is an ensemble/joint-law statement. It is not a lower bound on the
number of integer-valued queries needed to identify one unknown state;
a single integer can encode several bits. OBSERVER_COMPLETION gives a
literal example where state identification precedes linear-law recovery.

## HB4. Strong sources attain every budget

For each k>=1 set N=2k and use the strong source

    T_x = C3[TT_(4k) with its 2k disjoint adjacent pairs oriented by x, 1,1].

The two singleton vertices are fixed external port names. All T_x are
isomorphic and have H=2^(4k)+1. Halving the specified pairs inside the large
module gives C3[TT_(2k),1,1], with the same named anchors and unordered pair
fibres, independently of x. This is modulewise halving; the anchors are
retained, not contracted.

Group the 4k internal vertices into k successive four-vertex intervals.
In each interval add one reader, beaten by positions 0 and 2 and beating
positions 1 and 3. Keep the enlarged intervals in ordinal order. Put every
reader in the large module of the same C3 wrapper. This fixes every cut
independently of the hidden bits and preserves strongness.

For one five-vertex interval the four responses, indexed 00,01,10,11, are

|Pair bits|H|F2, two ordered bins, either empty allowed|
|---|---:|---:|
|00|3|48|
|01|5|64|
|10|5|64|
|11|9|84|

These tiny entries can be counted over the 5! permutations and four cuts;
`halving_response.py` does so independently of subset convolution. Their
alternating sums are 2 and 4. Ordered-bin counts multiply across ordinal
sums, by restriction and unique concatenation in each bin. For any S,

    H(C3[S,1,1]) = 3H(S)+Q2(S) = H(S)+F2(S),

where Q2 has two nonempty bins. Indeed a Hamilton path visits S in either
one or two runs; there are three one-run orders and one ordered two-run
skeleton. Consequently the full 2k-bit alternating sum after the k readers
is exactly

    2^k + 4^k > 0.

HB3 excludes every attachment with fewer than k new vertices. The displayed
attachment succeeds with k. Thus the minimum attachment depth separating
these two positive parity laws is **exactly k**, for every k, with strong
sources and targets and both exterior ports retained. Direct subset counts
independently check every state for k=1,2 (orders seven and twelve).

The inherited matching-cap moment theorem
[BC7](../boundary_carry_returns_20261103/PROOFS.md) bounds the number of
exponential modes for a fixed clone query. HB4 bounds a different resource:
how many actual new vertices can test a high-order correlation among old
module orientations. A short scalar recurrence and a growing joint-law
carrier can therefore coexist.
