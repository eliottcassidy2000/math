# Changing the marked host is not enough when both cores come from the binary-apex alphabet

**Statements.** The two independent core profiles need not have equal ratios. For every fixed finite anchor menu, if both substitutable cores are ordinal sums of one-apex tournaments, the number of distinct outputs through X is at most
$$
X^{log_3 2}exp(O(sqrt{log X})).
$$
Allow a *changing* host whose marked outlet is an apex over a transitive spine. If the middle core W is cyclic and both cores are ordinal sums of arbitrary one-apex tournaments, at least 2% of odd targets are absent in lower relative density. If U is cyclic as well the lower fraction missed exceeds 73%. This does not restrict hosts with a nontransitive marked deletion, cores outside that alphabet, or the separate lane where W is transitive.

## 1. Where the power 3 enters a two-slot count

Let h and Q respectively count Hamilton paths and two ordered covering paths, with empty slots allowed, and put P=Q−h. For a nonempty tournament, Q≥2h by putting a full path in one slot and leaving the other empty. More strongly Q≥(n+1)h by cutting each full path. Ordinal sums multiply both h and Q: split each of the two paths at the ordinal boundary.

A strong one-apex tournament of order n≥3 has a transitive spine and a marked apex whose outgoing binary word starts 1 and ends 0. Every cycle contains the apex. From the [odd-cycle packing formula](https://arxiv.org/html/2412.10572v3#S3.SS4.Thmcorollary20), an odd cycle of length l contributes 2 to h and 2^(n+2−l) to Q. In the latter count, a slot must contain the entire cycle but each remaining vertex independently picks its slot. An eligible pair of spine positions i<j, marked 1 then 0, has r=j−i−1 available interior vertices. Any even-sized subset of these gives exactly one odd cycle; the spine order determines its internal order. Its aggregate contributions to h and Q are
```
2 if r=0, otherwise 2^r;     2^(n−2−r) (3^r+1), respectively.
```
The empty packing contributes 1 to h and 2^n to Q. In particular the first and last spine vertices alone imply
$$
Q \ge 2^n+2\cdot 3^{n-3}+2.
$$
Thus binary freedom in the one-apex word meets ternary growth in its two-slot readout: an unselected interior position is free to occupy either complement slot, or it is selected for the apex cycle. Parity picks the even subfamily and gives the +1 correction.

## 2. Any fixed menu, independent ratios

Let D be the set of distinct simultaneous (h,Q) profiles of strong cyclic one-apex tournaments. The number of profiles in D obtainable at order n is at most 2^(n−3). Give the singleton profile (1,2) its own factor. Every ordinal sum of one-apex blocks has a profile encoded by a commutative multiset from D and singleton copies, by multiplicativity; duplicates in the final product can only lower the count.

Set alpha=log_3 2. For real s>alpha near alpha, the multiset series with weight Q^(-s) is bounded by
$$
Z(s)=(1-2^{-s})^{-1}\prod_{(h,Q)\in D}(1-Q^{-s})^{-1}
       =\exp(O((s-\alpha)^{-1})).
$$
For convergence and the estimate, sum the logs of the factors, then use `sum_(n≥3) 2^(n−3)(2·3^(n−3))^(−s)`, a geometric series of ratio 2·3^(−s). Profile duplicates across orders are harmless upper bounds.

For the [physical two-core host formula](../chanpak_quadratic_anchor_atlas_20261115/README.md), let the marked host coordinates be A,K,B, with B≥K≥1. For arbitrary U,W in the stated alphabet its output satisfies
$$
N=(A+K)h_Uh_W+P_U(BQ_W-Kh_W)\ge (B/4)Q_UQ_W.
$$
All terms are nonnegative and we used P_U≥Q_U/2, h_W≤Q_W/2. For each host, Rankin's inequality bounds the number of possible pairs by (4X/B)^s Z(s)^2. Use s=alpha+1/sqrt(log X) and sum over the fixed host menu. In particular this reaches a density-zero fraction of primes and prime cofactors of 7p, without any synchronized assumption.

## 3. Changing the binary host: exact target decoder and an absolute capacity gap

Choose any host R with R−q transitive. Any leading zeros in the spine give an initial ordinal block dominating the rest of the glued tournament, so they can be removed without changing the count. After their removal q reaches all of R; the glued tournament is strong. Mark b_i=1 when q beats the ith spine vertex (of m). The first bit must be 1 unless m=0: the first spine vertex is only reachable by a direct arc from q. Set t=sum b_i 2^(m−i) for one-based i. Then A=h(R), K=1 and B=1+t: the first slot can contain q alone, or have first spine vertex i with b_i=1, whereupon each of the m−i later spine vertices has two choices. The complement is transitive. Because the first bit is 1, each t>0 encodes exactly one word; t=0 is the singleton host.

By converse, the one-apex count A equals that for the reverse-complement word. In the [canonical suffix decoder](../chanpak_collatz_exact_addresses_20261114/BINARY_SUFFIX_DECODER.md), that word has low-bit address t, so A−1=p(t) with the all-order recurrence
```
p(0)=0,    p(2u+1)=p(u),    p(2u)=p(u)+u+(u mod 2).
```
For each fixed profile pair put a=h_Uh_W, b=P_UQ_W, C=P_U(Q_W−h_W)+2a. The absolute inverse is not a quadratic or norm equation:
```
N = C + b*t + a*p(t).
```
The cited suffix split applies verbatim with h=a,S=b and decodes t from a prescribed target in O_(a,b)(1+sqrt(N)) candidate tests. This produces a literal host, not just a congruence. It does not promise a solution. By the same harmonic law, the exact number of presentations for a fixed pair is (X/a) ln(1+a/b)+O_(a,b)(sqrt X).

### A quotient that halves the host capacity

A terminal run of 1s in the host word is a transitive ordinal suffix of R. Move this suffix to the front of W. This changes neither vertices nor arcs: q and all earlier spine vertices dominate the suffix, and the suffix dominates old W and U. It preserves the profile alphabet and the cyclicity of W. Consequently we may and do count only canonical hosts with empty word or final bit 0, that is, even t. This is an exact gluing redundancy, not an identification of merely equal outputs.

For a uniform capacity upper bound use only N≥C+bt and C≥b/2. If X<b/2 there are zero presentations; otherwise their number with even t is ≤X/(2b)+1≤3X/b. For each fixed pair the limsup divided by X is at most 1/(2b). To safely sum over infinitely many pairs, bound the tail by 3 times the sum of 1/b over that tail; the next calculation proves absolute convergence.

Put Z_c=prod_D (1−1/Q)^(-1). A cyclic ordinal profile uses at least one cyclic strong block. Thus the sum of 1/Q over the distinct profiles (h,Q) of such cores is at most 2(Z_c−1), including the singleton-factor geometric series. Their orders are at least 3, so P≥3Q/4. The sum of 1/P over transitive U is E=sum_(n≥1) 1/(2^n−1). Therefore the capacity upper bound per X for arbitrary U and cyclic W is at most
```
(Z_c−1) (E + (8/3)(Z_c−1)).
```
If U is also cyclic, omit E. Counts depend only on the profiles (h,Q) and t, so duplicates can only lower these bounds.

Here is a small exact certificate that Z_c<611/500. The rational product over all distinct strong cyclic one-apex profiles through order 15 is <3041/2500 by the two slot-pair formulas in §1 (4,115 profiles; see replay). Beyond that, reverse-complement preserves both h and Q, and its fixed words number at most 2^((n−3)/2), with equality at odd n only. For n≥16 the number of word orbits is hence at most (65/128)2^(n−3), because fixed words are at most 1/64 of the word count. The log tail is at most
```
(65/128) (81/16)(2/3)^16 < 4/1000,
```
from Q−1≥2·3^(n−3) and log(1+1/(Q−1))≤1/(Q−1). Since exp(z)≤1/(1−z), Z_c <(3041/2500)(1000/996)<611/500. Also E<1607/1000 by summing through n=12 and bounding the tail by 2/8191 (the tail terms are at most (8192/8191)2^(-n) for n≥13). These rational bounds give limsup absolute capacity <0.49 X, or <0.135 X if both cores are cyclic. As there are X/2+O(1) odd targets, the corresponding lower omitted fractions are at least 2% and 73%, respectively. This is not a claim about the full spectrum.

## Remaining opening

The separate lane W transitive, cores whose cyclic structure is not an ordinal sum of one-apex blocks, and anchors whose marked deletion is not transitive evade the variable-host bound. A host search should not freeze all of those escape coordinates simultaneously. The target decoder above permits host changes with fixed cores; the anchored linear progression in the two-core note permits changing a core at fixed host. Neither inverse supplies the missing structural coverage theorem.

Replay: `python replay_host_capacity.py`; independent subset/endpoint DP checks every printed profile through order 8 and 32 literal variable-host glues; see [receipt](replay_host_capacity.txt). The all-order estimates use the proofs, not propagation of checked ranges.
