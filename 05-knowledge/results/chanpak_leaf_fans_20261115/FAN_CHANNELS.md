# Exact fan invariance, and what happens to the unused ports

**Model.** Use the chronological tournament shorthand: adding a leaf at q means a new vertex x beats q and loses to every other old vertex. In the reversed-edge graph relative to chronological order this adds the undirected edge xq. It does not mean adjoining a leaf to an ordinary graph and counting dominating sets or independent sets. This is the unrestricted-parent version of the leaf geometry in the [older chronological-caterpillar work](../../../session199/hidden-structure/PROOF.md).

Write H=H(T). For a port q use S_q and d_q for total and odd-minus-even weighted starts; every nonempty path starting at q is weighted by the Hamilton count of its complement. Put K_q=H(T−q). The signed weight d_q is exactly the number of Hamilton paths ending at q. One way to see the latter is the squarefree path-polynomial identity in the [rank-one complement calculation](../chanpak_heisenberg_target_20261115/CUT_HEISENBERG.md): if H(X) is the subset Hamilton series, P_q and E_q its start and end path series, then H(X)P_q(−X)=−E_q(X); take the coefficient of the full vertex monomial. In particular 0≤d_q≤K_q, since deleting a terminal q injects into Hamilton paths of T−q.

## First leaf, exact activation, and geometric continuation

The [backward-edge gain](../chanpak_whip_seam_20261115/BACK_CROSS.md) gives H' = H+S_q−d_q. Define the defect with the sign convention
```
δ_q=H'−2H=S_q−d_q−H.
```
With q kept as the port, the new state is
```
S'_q=2S_q−K_q,   d'_q=K_q,    K'_q=K_q,
δ'_q=δ_q+2(d_q−K_q).
```
A nonempty path from q either omits x or appends x at the end of a path of old order at least two. If it omits x then x is a suffix sink in the complement (since q was used); if it appends x the weight is unchanged and its parity reverses. The singleton q is the sole old path that cannot append x. This proves both weighted-start updates. Deleting q leaves x a sink and preserves K.

So the *first* attach can lower the defect. It remains invariant precisely when the port already had d_q=K_q, equivalently all Hamilton paths of T−q end at a vertex beating q. After that first attach it always holds for every subsequent attachment at that port. For L=S_q−K_q and η=δ_q+2(d_q−K_q), the exact fan count after k≥1 leaves is
```
H_k = 2^k L−η.
```
In particular a stable unit fan with δ=±1 has H_k=2^k(H+δ)−δ. The original minus-one alone is insufficient: the literal chronological rows (256,509,505,497,481,449,129,257,64) have H=5 and (S,d,K,δ)=(6,2,3,−1) at q=0. The first count is 9, the defect becomes −3, and the second count is 15, not 17. A genuinely stable +1 site is q=0 in rows (6,12,8,1), giving 5,11,23,47,95,… exactly.

## The newborn never supplies a new +1 during a growing step

For the newly attached x, the even weighted starts are the old odd starts at q (prepend x before q). Thus
```
δ'_x = 2d_q−H.
```
If q is not a sink, an old Hamilton path ending at q has an injective partner not ending at q: remove q and reinsert it immediately before its first outneighbor along the remaining path. If that outneighbor is not first, the preceding vertex beats q by definition of first. Removing q recovers the original remaining path, so the map is injective. As H is odd, d_q≤(H−1)/2. Hence δ'_x≤−1 during every nonneutral append. If q is a sink then d_q=H and the append is neutral; the newborn is +1 only in the special H=1 case. In particular moving to the newborn is not a general way to refill the positive branch of a binary tree.

## Every individual fan leaf has an exact defect

Let the k new fan vertices be t_1,…,t_k in birth order, and focus on t_i. Put a=i−1, b=k−i, retain L,K=K_q,d=d_q,η from the original T. Its exact next-leaf defect is
```
δ_(t_i) =
  L·2^a (2·3^b−1−2^(b+1))
  +(K−η)2^b +(K+η)[b=0]+(d−K)[a=0]+η.
```
Here every older vertex other than q beats the fan, and each fan vertex beats q; the whole fan is transitive in birth order. Split a weighted even path starting at t_i into:
1. **Only fan vertices.** Choose an odd subset of the b later vertices (size l), and its complement contributes H_(k−1−l). Substituting the fan formula, including the correction H_0−(L−η)=d−K, gives 2E=L·2^a(3^b−1)−η(2^b−[b=0])+2(d−K)[a=0 and b odd].
2. **Old vertices only q.** The path enters q last, so choose an even number of later vertices. The contribution to 2E is K(2^b+[b=0]).
3. **An old q-start segment of order ≥2.** Earlier fan vertices can be unused or after the old segment; later ones can be unused, before q, or after the old segment. Both used runs have forced internal order. Half of the 2^a3^b assignments have the required parity, corrected by the sign sum, which vanishes for a>0 and is (−1)^b for a=0. Summing old path weights gives 2E=L·2^a3^b+(d−K)(−1)^b[a=0].

Add these contributions, use (−1)^b+2[b odd]=1, and subtract H_k. The 3^b versus 2^b channels are the effect of allowing a used fan vertex to participate in either run; pure fan prefixes cannot see that extra choice.

**Consequence.** On a stable +1 fan, L=H+1 and d=K. Every fan vertex except the newest has b≥1, with 2·3^b−1−2^(b+1)≥1, hence δ≥3. The newest has defect 2K−H_(k−1), tending to minus infinity. Moreover a stable +1 root must have at least two old outneighbors: if it only beats p, write A=T−q. Then its start data satisfy d_q=K−d_p(A); stability forces d_p(A)=0. The leaf gain at p gives H(T)=K+S_p(A), and therefore δ_q=−K. A sink is negative as well. The [cross-port law](CROSS_PORTS.md) now applies to every other old vertex p: at least one old q-outneighbor is distinct from p, so its positive 3^k coefficient forces δ_p to +infinity. **Sufficiently far along any fixed stable +1 fan, the only unit port is q itself.** This excludes rearming the other sign while remaining on that same fan indefinitely, even if any *fan* vertex is available to be picked. It is not an exclusion of a target-dependent move to another core before the fan grows.

Replay: `python replay_fan_channels.py` checks 2250 individual fan-port defects with independent literal Hamilton subset DP, as well as the original −1 counterexample and the positive ladder.
