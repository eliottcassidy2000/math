# Two guards: literal tournament seams versus Moore-deck distance

The [matching-fibre and one-apex inverse](FIBRE-TRANSFER.md) separates true outputs from multiple descriptions. This pass keeps four transfer coordinates fixed, compares moves in the two underlying spaces, and shows why the physical tournament move must live outside the graph payload.

## 1. The six-site guard is a pair, not a run

Let \(Y\) be the middle block in the [graph-address fork](../moore_merged_count_addresses_20261124/README.md#3-a-graph-preserving-physical-2-fork), and let \(O\subset Y\) be its apex outcontacts. For ordered partitions of \(Y\) into two directed paths, allowing empty lanes, define:

- \(Q\): all covers;
- \(R\): covers whose right lane is empty or begins in \(O\);
- \(J\): covers whose left lane is nonempty and ends outside \(O\).

The four guards are \(\sigma=(H(Y),Q,R,J)\); the reading is \(H(Y+a)\). For a binary prefix \(w\) ending in 0, with value \(B\) and one-apex reading \(h(w)\), the extended reading from the existing fork formula is
\[
 2^K(QB+R)+H(Y)(h(w)-1)+JB+H(Y+a).
\]
Thus a literal internal reversal with fixed marking \(O\), fixed \(\sigma\), and increment 2 in the reading works for **every** prefix and guard width.

The exact subset-DP census fixes the displayed Hamilton spine \(0\to1\to\cdots\to n-1\) of \(Y\), considers every optional arc marking and every fixed \(O\), and tests all reversals of optional internal arcs. Arc pairs are counted once in the increasing direction. It gives:

| \(n\) | Marked guard-preserving +2 edges | Longest directed chain in heights |
|---:|---:|---:|
| 2–5 | 0 | none |
| 6 | 48 | 2 vertices |
| 7 | 96 | 2 vertices |

All 48 six-site edges read \(111\to113\), with guard vectors \((29,320,158,162)\) or \((29,320,191,129)\). These are counts of **spine-marked labelled** moves, not of unmarked isomorphism classes. All seven-site examples are padding, with no new edge after trimming:

- A new dominating source, also beating the apex, leaves both readings and \(H(Y)\) fixed and sends \((Q,R,J)\) to \((2Q,R,2J+H(Y))\). To see the last entry, place the source in either lane; it creates one extra eligible type when it occupies a previously empty left lane.
- A new dominated sink, also beaten by the apex, sends \((Q,R,J)\) to \((2Q,2R,J)\) and leaves both readings fixed.

Hence the 48 source-padded moves and 48 sink-padded moves exhaust the 96 observed moves. The images are disjoint: an arc in both would remain an eligible five-site edge after removing both extreme sites, because each guard transformation is injective when \(H(Y)\) is fixed. This rules out simply making this particular common-spine guard one site larger to obtain three consecutive heights. It does not rule out different interfaces, apex-contact edits, or a reversal that does not preserve the chosen spine.

Setting \(K=1\) in the existing six-site fork gives
\[
 H_-=802B+29(h(w)-1)+427,\qquad H_+=H_-+2.
\]
This gives \(\Theta(X)\) descriptions below \(X\) in a fixed leading prefix. It is **not** an injective list of integers: \(B=659056,666074\) both give 546598033. The replay verifies this collision. At width 23, with prefix 1010 and last bit 0, 262144 words yield 262012 lower integers. This finite count is not a density theorem. The longer guard was for uniquely decoding the payload from the height, not for making the physical seam. A new selector could use two outputs to recover the payload even when neither is individually injective; the different linear \(B\)-slopes of two fixed guard choices make their paired height vector injective. It would still need an exact target inverse for a chosen scalar height.

## 2. Frozen probes give a minimum Moore-deck edit distance

Use the [rooted matching deck](../moore_merged_count_addresses_20261124/README.md#2-exact-count-addressed-encoding) of a Moore graph of degree \(k\), girth 5 and diameter 2. At a fixed labelled root the \(k\) branches have width \(k-1\); every pair of branches has one perfect matching. Normalize all matchings involving branch 0 to the identity. Compare **two valid decks** with these same root/heads/labels/normalization. Draw a support graph on the branch indices, marking a branch pair when its matching differs. Let \(t\) be the number of nonisolated branch indices. The other \(k-t\) branches are frozen, including branch 0.

**Probe bound.** At any satellite \(z\) with a changed neighbor, at least \(k-t+1\) of its old neighbors are replaced. Indeed choose one new neighbor \(v\) in branch \(j\). Each frozen branch \(\ell\) has a unique vertex \(w_\ell\) adjacent to \(v\), equally in both decks. In the new graph \(z,v,w_\ell\) is a two-step path; \(z,w_\ell\) are not adjacent by girth 5 and this nonadjacency is unchanged by frozen incidence. Their unique old common neighbor \(u_\ell\) is a deleted neighbor of \(z\); otherwise the new graph would have two common neighbors. Moreover the \(u_\ell\) are all distinct: two repetitions would form the old square \(v,w_\ell,u_\ell,w_{\ell'},v\). None of them can be the deleted old partner of \(z\) in branch \(j\): that partner and \(v\) already have the branch head as their common neighbor and cannot also share \(w_\ell\). Hence the \(k-t\) probes demand \(k-t\) different deletions besides the old partner in \(j\).

Every active branch therefore has support-degree at least \(k-t+1\). Since it has at most \(t-1\) support neighbors, \(2t\ge k+2\), and there are at least \(t(k-t+1)/2\) modified matchings. In particular for \(k=57\):
\[
 30\le t\le56,\qquad \#\text{modified matchings}\ge56.
\]
At \(t=30\) the lower bound is 420; each of the 30 branches would have at least 28 changed incident matchings. Both facts are independent of automorphism hypotheses. Freezing all matching rows of any 28 designated branches (including 0) guarantees **at most one** valid normalized completion. Every changed permutation removes at least two edges and adds two, so there are at least 224 literal reversed arc positions between distinct valid native decks **in the same displayed site order**, if a common complement spine is fixed. It says nothing about whether even one such deck exists. It is a separation result, not a minimum-distance construction.

The Petersen normalized deck is unique by this bound (\(k=3\) permits no nonzero \(t\)). Hoffman–Singleton is a finite control: for one fixed root and hub, regrouping its \(6!\) non-hub branch orders and \(6!\) common labels gives 720 coordinate descriptions at distance 0, 10800 at 8, 72000 at 12, 216000 at 14, and 218880 at 15. The minimum positive distance there is 8 (the general probe bound gives 6). These figures are descriptions of a **fixed graph** under marked coordinate changes, not counts of distinct graphs. The replay constructs the graph by the standard two pentagon-families formula \(y=x+ij\bmod5\) (pentagons with step 1 and pentagrams with step 2 internally).

## Connection

A true graph-payload edit can require many simultaneous contacts; the guarded tournament seam changes a separate arc whose four exported coordinates stay fixed. More auxiliary width is not automatically more scalar freedom: the seven-site census proves it does nothing in the first extension of this interface. For the next pass, test a different port signature or graph-safe multi-contact transfer rather than search one larger padding. Conversely, for a Moore search, freeze exact probe rows and seek a global solution or obstruction in one shot instead of treating close payload edits as likely valid.

Replay with `node guarded-fibres.js 7` (the full seven-site census takes longer than the six-site check). Run `node moore-deck-distance.js` for the normalized deck sample and full Hoffman–Singleton marked relabel census.
