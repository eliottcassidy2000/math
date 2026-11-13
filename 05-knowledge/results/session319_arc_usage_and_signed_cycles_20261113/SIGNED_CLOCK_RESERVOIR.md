# Equal-clock signed Collatz reservoirs: exact carry, hull, and bounded inverse

**Proved for a fixed finite bank of equal-clock valuation blocks.** No new claim of completeness of positive or negative Collatz cycles is made.

## Source lineage

- [Pillai and the cycle gate](../collatz_mod6_20260921_pillai_convergents_cycle_gates.md): the inherited exact cycle equation, integer obstruction and the clocks \(2/1,1/1,3/2,11/7\) for the known signed cycles. For \(3n+1\), the known positive odd cycle is \(1\), while the known negative odd cycles are \(-1\), \((-5,-7)\), and the seven-cycle through \(-17\). “Known” is not exhaustive. Sign reflection conjugates \(3n+1\) on the negative side to \(3n-1\) on the positive side.
- [The negative-cycle shadow](../collatz_negative_cycle_shadow_20261004.md): inherited negative \((1,2)\) shadow obstructs a bounded positive common-future bank with explicit intercept guard. Crucially a rational point can feed an integer endpoint, so integrality must be retained.
- [Signed resonance and pair mass](../collatz_signed_resonance_pair_mass_20261013b.md): signed phases and exact counts are separate from asymptotic or finite-state pictures.
- [Radix block transport](../collatz_radix_families_20261018p.md): retain exact carry and low-digit branch, not only leading size. This motivates the two endpoints and digit congruences below.

## 1. Independent switches, but a common clock

For odd \(n\), let \(U(n)=\operatorname{oddpart}(3n+1)\) on either sign. A word of exact valuations \(a_1,\ldots,a_r\) has total \(s=\sum_i a_i\), formal map
\[
F_w(n)=\frac{3^r n+B(w)}{2^s},
\qquad
B(w)=\sum_{i=1}^r 3^{r-i}2^{a_1+\cdots+a_{i-1}}.
\]
Swapping two adjacent letters \(a,b\) starting at \(i\) changes only the term \(i+1\) in the sum. Thus, for the original minus the swapped word, the exact carry difference is
\[
3^{r-i-1}2^{a_1+\cdots+a_{i-1}}(2^a-2^b).
\]
Disjoint pairs give independent additive binary carry switches: every choice keeps the clock and exactly the same exponent sum before each later pair. The switch bank is a literal subset-sum system, not just an analogy. It does **not** say that all resulting carry targets satisfy the native return gate \(2^s-3^r\mid B\), or that arbitrary switches preserve a supplied source's exact valuation cylinder.

## 2. Equal-clock concatenation gives an exact weighted hull

Fix a finite bank of words all of the same clock \((r,s)\). Put \(p=3^r,\ q=2^s\). Since \(r,s\ge1\), \(p,q>1\) are coprime and unequal. Write \(z_i=B(w_i)/(q-p)\) for the single-block rational fixed points. For a concatenation \(w_{i_0},\ldots,w_{i_{k-1}}\), compose the maps. The unique formal fixed point is
\[
z=\frac{\sum_{j=0}^{k-1}p^{k-1-j}q^jB(w_{i_j})}{q^k-p^k}
 =\frac{\sum_{j=0}^{k-1}p^{k-1-j}q^j z_{i_j}}
        {\sum_{j=0}^{k-1}p^{k-1-j}q^j}.
\]
The last identity uses \((q^k-p^k)/(q-p)=\sum_j p^{k-1-j}q^j\). All weights are positive. In particular the start and every block endpoint of a closed itinerary, by cyclically rotating the blocks, lie in the fixed compact interval \([\min_i z_i,\max_i z_i]\). The sign of \(q-p\) determines which sign side the formal fixed point is on, since all carries are positive. Thus equal-clock rearrangements cannot change a positive candidate to a negative candidate.

This yields a finite native inverse for a bank of arbitrary finite size and unbounded concatenation length. Use the odd integer points of that interval as vertices. At each vertex retain a word only if it matches the **actual** next \(r\) valuations and its endpoint remains in the interval. Every closed bank trajectory appears as a cycle in this finite directed graph, and every such cycle is an actual signed Collatz cycle (expanded into its original odd steps). There is at most one outgoing edge per integer, because the next \(r\) valuations are unique. Consequently it is enough to test simple cycles of length at most the number of odd integers in the interval, regardless of the allowed concatenation length. This does not bound the odd vertices internal to one block by the same hull, and it gives no reduction to a finite bank for arbitrary Collatz trajectories.

## 3. Two endpoint words: an exact obstruction stronger than interval density

Suppose a two-word bank has integer single-block fixed points \(a<b\), and put \(d=b-a\). For a mixed word, a putative integer start \(n=a+t\) requires \(0<t<d\). In the preceding formula, encode the choice of \(b\) by \(\epsilon_j\). Then
\[
d\sum_j\epsilon_jp^{k-1-j}q^j
=t\sum_jp^{k-1-j}q^j .
\]
Reducing modulo \(q\) forces \(t\equiv d\epsilon_0\pmod q\). Reducing modulo \(p\) forces \(t\equiv d\epsilon_{k-1}\pmod p\). Hence if \(d\le\max(p,q)\) no interior integer is possible, at **any** length. This obstruction precedes the native valuation test and does not assume the endpoint integers are consecutive.

There is also a length-sensitive sieve when the spread is larger. Write \(t/d=u/v\) in lowest terms. Then \(v\) must divide \(S_k=\sum_jp^{k-1-j}q^j\). Because \(S_k\) is a unit modulo both \(p\) and \(q\), any factor 2 or 3 in \(v\) permanently excludes the target. If a prime \(\ell\mid v\) does not divide \(p-q\), reduction gives \((q/p)^k\equiv1\pmod \ell\), so the corresponding order divides the number of blocks. These are necessary, not sufficient, conditions.

For example the negative \((1,2)\) and \((2,1)\) blocks have \((p,q)=(9,8)\), carries 5 and 7, and single-block fixed points \(-5,-7\). Any mixture closes formally strictly between those two values, but cannot close on an integer by the congruence theorem. Both pure words recover the two actual phases of the same cycle. Likewise at the seven-cycle clock \((p,q)=(3^7,2^{11})\), the seven actual phase points are within \([-91,-17]\). More directly than interval density, reduction of a fixed point modulo \(p=2187\) forces it to agree modulo \(p\) with the final block's single-block point, since all seven points are integers. The interval has width 74, so the fixed point equals that final point. The exact next seven valuations at a cycle phase select its single correct rotation and return to that phase. Thus *even allowing all seven rotations as separate letters*, no mixed closed block word exists. This last argument applies to any number of integer endpoints in an interval of width less than \(\max(p,q)\): choose the first or last block using the appropriate modulus; when these endpoints themselves are actual single-block fixed points, the resulting integer fixed point selects its unique actual block at all subsequent phases.

## 4. Transfer and its limit

The applicable analogy is a reservoir of local switches with stable clock, not equality between the Collatz and tournament state spaces. The Collatz carry exposes a failure mode that an exact tournament proof must rule out separately: independent changes can give rich formal ranges and preserve a coarse return type, while a native integrality or readiness gate blocks every mixed return. For a proposed tournament reservoir, first write the conserved clock and the exact signed response of each independent switch. Then show that the combined exact response lies in the allowed target fiber **and** its final interface can be read natively. The equal-clock lemma cannot be used to conclude that either property follows from the other.

## 5. A synchronized shadow can preserve a gap while equality remains impossible

This is the closest exact analogue to the reservoir switch, but it has a different output and does not transfer a tournament recursion. Let equal-clock blocks have actual integer fixed phases \(a,b\), with \(a\ne b\). For every \(k\ge1\), a shared displacement \(L\) divisible by \(2q^k\) gives actual repeated words from both starts \(a+L,b+L\). To check the guard, two odd sources congruent modulo \(2^{s k+1}\) have identical valuation words through total cost \(s k\) if one source realizes that word: after every exact step their difference is multiplied by \(3/2^{a_i}\), and enough factors of two remain to preserve the next valuation. The two exact outputs are
\[
a+(p/q)^k L,\qquad b+(p/q)^k L.
\]
Thus the leading exponential part is synchronized, and the gap \(b-a\) survives at every depth. Taking \(L\) large positive gives this comparison entirely on positive integers even when the reference phases are negative. For the \(-5,-7\) phases, take \(L=2\cdot8^k M\), \(M\ge1\). The positive starts are \(2\cdot8^kM-5,\ 2\cdot8^kM-7\), the positive outputs are \(2\cdot9^kM-5,\ 2\cdot9^kM-7\), and the output gap is exactly 2.

By contrast, equality after repeating the two blocks \(k\) times, starting at any integer sources \(n_a,n_b\), requires
\[
n_b-n_a=(b-a)\bigl(1-(q/p)^k\bigr).
\]
Coprimality implies \(p^k\mid b-a\). In particular there is **no exact endpoint equality** at any positive depth for \(-5,-7\); and any fixed two distinct integer anchors have a finite upper bound on the depth of exact equality. This criterion is necessary even for formal affine maps on integer sources, so native valuation constraints cannot rescue a failed divisibility test. In the positive \(3n+1\) fixed phase \(1\), there is no second known positive integer phase of the same primitive clock. The comparison is therefore a positive-source shadow of known negative structure, not an assertion of a positive nontrivial cycle.

The tournament charts are better in one specific sense: they realize identical absolute \(H\) on a selected dyadic ray with distinct reservoirs, where the two Collatz phase charts above can only preserve a gap. They are worse in another precise sense: they presently lack the Collatz word's explicit native cylinder at the output. The nonzero tournament \(Q_2\) discrepancy and the Collatz residual \(b-a\) both show that agreement of a large-scale term does not identify the continuation data. The Collatz congruence test is a model for which state coordinate to keep, not an identification of \(Q_2\) with a Collatz phase.
