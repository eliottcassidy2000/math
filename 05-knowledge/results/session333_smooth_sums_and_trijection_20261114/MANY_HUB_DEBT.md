# Many hubs: the full-overlap factorial and the sharp next scale

Consider a tournament consisting of a transitive spine plus a fixed
number k≥1 of outside hubs. A full-overlap frame is a Hamilton path
order of the hubs whose k+1 spine runs permit β spine vertices to
be freely assigned to every run, after fixing the necessary extreme
vertices. For fixed k, as β→∞ the minimum ordinary Hamilton count
over such tournaments is

    k! (k+1)^β + Θ_k((sqrt(k(k+1)))^β).

The exact two-hub family and stronger finite checks are in
[TWO_HUB_FRAME_DEBT.md](TWO_HUB_FRAME_DEBT.md).

## Lower bound: spacetime order does not kill the other copies

For each hub j, let a_j be the spine predecessor of that hub and
b_j its spine successor. Each is locked by the frame. The β fully
eligible sites lie after every b_j and before every a_j; in
particular all 2k anchors are distinct. Keeping only these anchors,
the β sites, and the hubs gives an induced core.

For **every** permutation of the hubs, preserve each hub's two
locked joins a_j→hub_j→b_j. Every newly adjacent pair of joins
is compatible because all b_j lie before all a_j along the spine.
The sum of per-run eligibility for an ordinary non-anchor site i is

    Σ_j ([i<a_j]+[i>b_j]) − (k−1),

using [b<i<a]=[i>b]+[i<a]−1 when b<a. Thus even the multiplicity
for each site is independent of the ordering of the hubs. In the
induced core there are k! distinct marked frames of count (k+1)^β
each, regardless of the directed edges among the hubs. Their
contribution is k!(k+1)^β.

Now pick any β-site v, with l earlier and r later sites. For any
hub j, either v beats the hub, and can replace its right/incoming
anchor a_j, or the hub beats v, and it can replace its left/outgoing
anchor b_j. These replacements keep all b_j before all a_j. In
the first case the earlier sites retain k+1 choices and the later
sites retain k; in the second case the two groups exchange roles.
The displaced old anchor has at least one choice in an outer run.
There are k distinct single-anchor replacements for this v, and
each can be arranged in all k! hub orders. Different replacements
and different v give disjoint frames (same marked order, different
fixed extrema), all separate from the original k! frames. Therefore

    H−k!(k+1)^β ≥ k! k Σ_(j=1..β)
                   min((k+1)^(j−1) k^(β−j),
                       k^(j−1) (k+1)^(β−j)).

The central summand alone is Ω_k((sqrt(k(k+1)))^β). Hamilton-path
insertion extends all these distinct paths injectively back to any
larger host.

## Upper bound: physical blocked hosts

Choose spine order

    (b_1,...,b_k), (p early sites), (q late sites), (a_1,...,a_k),

with p+q=β and |p−q|≤1. Let b_j be beaten by exactly hub j,
and a_j beat exactly hub j. Every early site beats all hubs,
and every hub beats all late sites. Orient hub-to-hub arcs
arbitrarily. The tournament is strong: each hub's b_j→...→a_j
backward-return cycle contains the bulk, and the intervals cover
the full spine and overlap. The k! prescribed frames are legal.

For the remaining Hamilton words, classify using the order of
the hubs and the extremes of their spine runs, with missing
endpoints allowed in a run that is empty. For fixed k there are
only O_k((β+1)^(2k)) such frames. In all count estimates, choices
of the 2k original anchors not locked by a particular frame cost
at most a fixed factor (k+1)^(2k).

If any run is empty or an interior run is a locked singleton, or
if a frame uses a right/incoming anchor among the original left
anchors or a left/outgoing anchor among the original right, then
at most k choices remain at every bulk site. These frames cost
O_k((β+1)^(2k)k^β) in total. The same is true if some right anchor
is in the early block and some left anchor is in the late block:
every site is either after the former or before the latter and
therefore misses some run.

The other noncanonical case, up to reversing sides, has at least
one right anchor in the early block and none of those left-side
deviations. All left anchors are then the prescribed b_j. If the
earliest right anchor in the early block lies t places before its
end, there are O_k((t+k+1)^k) possible right-anchor lists: all
their early entries belong to that tail of length t+1. All bulk
sites after the first replaced anchor have at most k options.
Thus the total of this case is bounded by a constant times

    (k+1)^p k^q Σ_(t≥0) (t+k+1)^k (k/(k+1))^t,

and the geometric-polynomial sum converges for fixed k. The
mirror case contributes O_k(k^p (k+1)^q). The low-choice polynomial
term above is smaller than these because p and q differ by at
most one.

The k! original frames are the only ones allowing all bulk sites
the full k+1 choices. Their contribution is exactly k!(k+1)^β.
Taking the two sided bounds proves the claimed sharp scale. In
particular augmenting the number of hubs does not make a full
overlap coefficient independently selectable: it increases the
unavoidable multiplicity from 2 to k!.
