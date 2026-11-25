# Cut cancellation: partial matchings and exact two-step native seams

First overnight Chan–Pak ↔ degree-57 pass. It improves the uniform native lower bound and gives an exact small-case test of what that large baseline can and cannot say about a one-arc difference.

## 1. Count partial matchings, not only a fixed maximum matching

All degree-57 conclusions are conditional on the existence of the Moore graph. Keep the native chart and its marked mate from [the preceding pass](DEGREE57-FACTORS.md), and split any complement spine into balanced halves of size 1625. For opposite sets of equal size \(\ell\), the shifted Moore spectrum gives
\[
 e(X,Y)\ \ge\frac{115\ell^2-15\ell(3250-\ell)}{6500}
 =\frac{\ell(\ell-375)}{50}.
\]
Greedily **choose edges** until making a crossing matching of size 1250. At the successive choices the unused sides have sizes \(\ell=1625,1624,\dots,376\). Thus there are at least
\[
 P=\prod_{\ell=376}^{1625}\left\lceil\frac{\ell(\ell-375)}{50}\right\rceil
\]
ordered sequences of choices, and at least \(P/1250!\) (unordered) matchings. For each matching, use the ordered backward blocks, a fixed Rédei prefix and suffix on the unmatched sites, and the dummy-anchor switching lemma. The forbidden crossing list still has at most 52929 entries, with degree at most 59. Among the \(1250!\) ways to order the blocks, at least a fraction \((1012/1013)^{52929}\) are valid. Now multiply. The factorial cancels, rather than appearing as a loss:
\[
 H(T),H(T^\ast)\ge P(1012/1013)^{52929}>2^{16256}.
\]
The injection is exact even though the unmatched sets vary. The first and last 375 sites in each expanded path are the unmatched prefix and suffix; the middle 2500 split at fixed positions into the ordered backward pairs. Thus the path recovers the matching and its block order. Count the mate in its converse. The upper bounds \(2^{28860}\) and \(2^{28863}\) remain unchanged. The uniform exclusion range for \(E_p\pm1\) improves to \(2\le p\le8128\), since \(E_p+1<2^{2p-1}\). The stronger selected-spine floor \(2^{18706}\) is a different, existential guarantee and still stands.

### A huge class of paths cancels in the exact single-arc difference

Take a within-half arc for which at least one endpoint has a crossing graph neighbor. Fix that backward pair first, then apply the same partial-matching argument, starting the product at the top size 1624. For each resulting matching the chosen endpoint is in the middle, not either of the unmatched paths. No within-half arc occurs in the middle or in its joins, and neither unmatched path can contain the chosen arc. Thus the same expanded paths survive even if that arc is reversed. Their number exceeds \(2^{16241}\) in either native chart (and its corresponding single-arc edit). This is a common baseline, not a bound on the difference. In particular, it does not rule out \(\pm2\) between enormous heights. At most 47 vertices on either side can have no crossing neighbor, by the balanced matching result, so the argument applies to most internal adjacent contacts.

## 2. Read the exact difference by side changes

For any tournament with a fixed balanced cut, define \(F(z)=\sum_P z^{j(P)}\), with \(P\) a directed Hamilton path and \(j(P)\) the number of steps remaining on the same side of the cut. A within-half edit cannot change \(F(0)\), and its exact Chan–Pak step is \(\Delta F(1)\). There is a general signed version: for a within-half arc \(u\to v\), \(\Delta F=F_{\rm reversed}-F=z(Q_{vu}-Q_{uv})\). Here \(Q_{uv}\) counts paths of the marked contraction of the forced block \(u\to v\): entering uses the arrow into \(u\), leaving the arrow out of \(v\), and the single block occupies their original side. \(Q_{vu}\) uses the opposite endpoint conventions. Both polynomials grade by same-side steps of the contracted path. This follows by cancellation for paths not using the edited pair, and collapsing the one forced step in each surviving direction. The two contracted objects need not be tournaments.

For the Petersen graph, the replay takes all 796 backward masks supplied by complement spines on ten sites. It tests within-half adjacent swaps that preserve legality of the displayed spine. Each such swap reverses exactly one tournament arc and decodes the same Petersen graph. The mask edges below count undirected marked swaps up to graph automorphism (not isomorphism classes of unmarked tournaments). There are no loops because swapping an adjacent pair and preserving the same mask would make the transposition a graph automorphism, impossible for Petersen: the pair is nonadjacent and has only one common neighbor, so their neighborhoods differ.

| Cut edges | Crossing perfect matchings | \(F(0)\) | Spine masks | Within-half swaps | Absolute height gap 2 |
|---:|---:|---:|---:|---:|---:|
| 5 | 1 | 240 | 8 | 0 | 0 |
| 7 | 0 | 42 | 94 | 86 | 0 |
| 9 | 1 | 124 | 150 | 220 | 0 |
| 11 | 2 | 156 | 544 | 1024 | 10 |

All ten gap-2 marked swaps lie in the same cut type; there are five distinct scalar pairs:
\[
10325\leftrightarrow10327,\quad
11305\leftrightarrow11307,\quad
11641\leftrightarrow11643,\quad
13289\leftrightarrow13291,\quad
13787\leftrightarrow13789.
\]
The 1330 moves also include 18 with height gap zero. The fully alternating reading is not simply a function of the number of backward perfect matchings: an alternating path starting on the left needs only four backward pairs, so even a cut with no perfect matching contributes. In particular, a simple expectation that larger alternating baseline implies more small steps does not fit the observed table.

Two examples show **different mechanisms** at the same shared \(F(0)=156\). In increasing height direction, one seam \(11305\to11307\) has
\[
 \Delta F=11z^2-12z^3+2z^4-5z^5+4z^7+2z^8.
\]
Another, \(13289\to13291\), has
\[
 \Delta F=-48z-261z^2-242z^3-151z^4+127z^5+353z^6+147z^7+77z^8.
\]
Thus a literal +2 can be either quiet at the first nonconstant layer or the residue of hundreds of changes in several opposing layers. The fully alternating layer supplies a protected baseline; exact control must come from the signed contraction polynomial. The cut data give a small-model filter, not evidence that the four table rows exhaust anything at degree 57.

## Next transfer tests

- Look for a structural explanation of why the 11-edge cut alone admits these Petersen seams; first compare the two forced-block polynomials on its automorphism orbits. Test the same statistic on a small censused section of Hoffman–Singleton complement orders or exact phase charts with known \(\pm2\) edits.
- A sufficient rule for scalable exact +2 is \(Q_{vu}-Q_{uv}=2+(z-1)R(z)\) with the cancellation coming from a reversible pairing local to the readout, not just from value substitution. Search for such a rule in a guarded rail, then test which of its cut constraints are simultaneously compatible with Moore square-free incidence.
- Improve the uniform floor further by counting many extendable partial matchings or by changing the number of leftover sites. The current product only uses edges forced by spectral mixing at **every** greedy step.

Reproduce with `python3 partial-matchings.py` and `c++ -O3 -std=c++17 petersen-layers.cpp -o petersen-layers && ./petersen-layers`.
