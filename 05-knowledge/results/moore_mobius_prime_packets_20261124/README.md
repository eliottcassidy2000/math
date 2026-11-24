# Moore cycle-lift heights: Möbius inversion, actual prime packets, and a local budget

Only degrees 2, 3 and 7 are replayed from existing Moore graphs. Everything about a degree-57 Moore input is conditional. The arithmetic deck of permutation partitions is finite-exact but much larger than the set of jointly realizable graph readings.

## Two cycle charts with exact Chan–Pak overlap

At a root, each pair of leaf branches of a Moore graph is joined by a perfect matching. Pick three branches and traverse their matchings cyclically. The return permutation on $d-1$ leaves has no fixed points (a fixed point would be a triangle), and a cycle of length $\ell$ lifts to a graph cycle of length $3\ell$. Put each cycle in traversal order in a tournament block, with all arrows forward except the first–last arrow. A block on $m$ vertices has $1+2^{m-2}$ Hamilton paths: there is one avoiding the reversed arrow, and every other path assigns each internal vertex independently to the increasing segment before or after that arrow. The ordinal join gives
$$
 H_3=\prod_\ell f_3(\ell)^{c_\ell},\qquad f_L(\ell)=2^{L\ell-2}+1,
$$
where $c_\ell$ is the return permutation's cycle census. This is a literal height of the tournament *chart*, not the Hamilton count of an orientation of a Moore graph.

For a graph pentagon, take its five sets of exterior neighbors. Consecutive sets have no contacts, and each opposite pair is perfectly matched (the unique length-two route to that central vertex). The resulting pentagram lift has $5\ell$-cycles with $\sum\ell c_\ell=d-2$. Its chart has $H_5=\prod_\ell f_5(\ell)^{c_\ell}$.

| Moore degree | Three branches | Pentagon exterior |
|---:|---|---|
| 2 | not applicable | $1$ (empty product), one pentagon |
| 3 | $17$, all 10 rooted triples | $9$, all 12 pentagons |
| 7 | $17^3=4,913$, all 1,750 rooted triples | $9^5=59,049$, all 1,260 pentagons |
| 57 (conditional arithmetic minima) | $17^{28}$, all transpositions | $9^{55}$, separate exterior pentagons |

The four small nontrivial heights already arise in the repo's [older one-apex bank](../../../01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md): $g(16)=17,\ g(8)=9,\ g(5024)=4913,\ g(83922)=59049$. The two large minima also have automatic ordinal-join witnesses, without any graph. Their graph significance is the cycle structure they encode, not the existence of a separate tournament of that height.

The minimum for $H_3$ follows because $f_3(\ell)/17^{\ell/2}$ increases strictly for $\ell\ge2$. The least nonminimal even part is $4$, with multiplier $1025/17^2$; odd parts come in pairs, and two 3's cost only $129^2/17^3$, strictly less. Thus the next partition value is $129^2 17^{25}$. For $H_5$, replacing $\ell$ singleton parts with a part of length $\ell\ge2$ raises height, since $2^{5\ell-2}>9^\ell$. Neither minimum is asserted for any degree-57 root or pentagon.

## Two Möbius transforms before counting any residues

The old [closed-walk/Witt calculation](../../../01-canon/theorems/THM-502-the-closed-walk-census-ladder-and-coefficient-law.md) recovers primitive periodic walks from traces; a general tournament walk can still revisit vertices. Here the operator is a **permutation** $P$: its primitive walks are precisely the matching cycles. Exact fixed-point counts satisfy
$$
 F(m)=\operatorname{tr}P^m=\sum_{\ell\mid m}\ell c_\ell,
 \qquad
 \ell c_\ell=\sum_{a\mid\ell}\mu(a)F(\ell/a).
$$
This is the old [zeta/moment mechanism](../../../01-canon/theorems/THM-509-the-two-layer-baby-hodge-dichotomy.md) on a carrier with no hidden overlap corrections.

There is a second inversion directly on factors. Put $E=L\ell-2$ and let $b(E)=c_\ell$ on the supported exponents, zero otherwise. Then
$$
 2^E+1=\prod_{\substack{a\mid E\\E/a\text{ odd}}}\Phi_{2a}(2),
 \qquad
 W(a)=\sum_{j\ge1,\ j\text{ odd}} b(aj),
 \qquad
 b(E)=\sum_{j\ge1,\ j\text{ odd}}\mu(j)W(Ej).
$$
Here $W(a)$ counts formal cyclotomic packets. One must **not** identify a packet with a distinct rational prime. If a prime has exact order $2E$ for base 2 it divides $2^{E'}+1$ iff $E'/E$ is an odd integer; on those aliases the valuation is $v_p(2^E+1)+v_p(E'/E)$ by LTE. This is a multiplicity rule, not a parity rule. For each of the 55 lengths in each chart the replay verifies that there is a prime in its block factor absent from all smaller factors: iteratively strip all gcds with smaller factors; a nonunit remainder certifies existence. Thus descending on lengths reconstructs the exact partition from $H_L$, even though several lengths can share a prime. In the 56-sheet three-branch chart all **75,547** arithmetic partitions give distinct heights.

Actual small packets, before taking any quadratic character:

| Three-branch cycle length | Block factor |
|---:|---:|
| 2 | $17$ |
| 3 | $3\cdot43$ |
| 4 | $5^2\cdot41$ |
| 5 | $3\cdot2731$ |
| 6 | $65,537$ |

The degree-57 top bands illustrate what would be lost by compressing to the sign.
- $v_3(H_3)$ counts **all** odd-length cycles with multiplicity: $L=3$ gives $3\ell-2\not\equiv0\bmod3$ and LTE gives valuation exactly one when $\ell$ is odd.
- $v_{17}(H_3)$ counts cycles with $\ell\equiv2\bmod8$, each with valuation one on this width. The first alias of length 2 is length 10; if a marker appears 25 times, there is not room for any such alias. Thus the exact high bands are $28$: all 2's; $27$: impossible; $26$: a 4 and 26 twos; $25$: a 6 or two 3's together with 25 twos. The latter two have $v_3=0$ or $2$, respectively. Notably the **second-smallest height** is in the $25$-band, not the $26$-band.
- Modulo 7 each factor is 3, so the **full** residue gives cycle count modulo 6. Taking its Legendre symbol only recovers monodromy sign. A useful quick consistency check: products over the four triples in any four branches must be squares modulo 7 because each matching sign is counted twice. The test does not constrain the nonabelian part of the matching data.

There is even a factoring-free recovery of the number $c$ of cycles from a valid degree-57 height. For width $n=56$, $2^{3n-2c}<H_3\le 2^{3n-2c}(17/16)^c$. Intersect this one-sided window with $H_3\equiv3^c\bmod7$. Since $c\le28$, the windows for counts congruent modulo 6 are disjoint. This shortcut recovers how many, not which lengths; the prime and packet steps recover the latter.

## Local multiplicities must fit one common graph

Each length-two cycle in a three-branch return gives an induced Petersen on the root, the three parents, and the six cycle vertices. A length-one exterior component likewise gives the central pentagon plus an exterior pentagon and five spokes, another induced Petersen. Conversely in any such Petersen, relative to a root or a pentagon, there is one of the corresponding cycles.

At fixed root $r$ and fixed unordered pair of spokes $i,j$, consider the pentagons through those spokes. In an induced Petersen containing these spokes, exactly two pentagons use both; each arises from one of the two matched leaf pairs between the two branches. Hence the following **exact local budget** uses the decoded multiplicities, not their signs:
$$
 \sum_{C\ni(r;i,j)} c_1(C)=2\sum_{k\ne i,j} c_2(r;i,j,k).
$$
Summing over spoke pairs gives $\sum_{C\ni r}c_1(C)=6\sum_{\tau\text{ at }r}c_2(\tau)$. At degree 57, subtracting each side from its own maximum gives the same local and rootwise identities for deficits (55 minus the exterior singleton count, 28 minus the triple transposition count). In particular the all-minimum triple and all-minimum pentagon conditions are equivalent globally. The local budget is stronger than that equivalence: it tests partly specified decks.

## What was borrowed from the RH proof

[Alpöge’s announcement](https://alpo.ge/remark.pdf) singles out **creating multiplicity** in the spirit of Stepanov and then controlling a smoothed prime sum. The [verification and classical reading guide](https://www.math.uwaterloo.ca/~mrubinst/alpoge_claude_grh/) distinguish the machine-checked RH/Dirichlet-GRH statements from the still-in-progress human exposition. The concrete transfer here is to look for **prime-power divisibility with exact multiplicities**, then force the independent local readings under a common budget. For instance a $17^{25}$-divisibility gate gives a short length budget; the actual 3-adic valuation selects the branch that has two 3-cycles; the spoke-pair equation checks the multiplicity against independent pentagon charts. This improves on borrowing only the sign or a vague square-root heuristic. No RH estimate is used or claimed to control primes dividing these special finite heights.

For a second instance of the local-to-global transfer, see [the array completion and spectral prediction](PENTAGON.md). Replay: **python3 replay.py**, which checks both chart types, the local identity, the small array constraints, the prime packet existence certificates, and the entire width-56 arithmetic deck.

Follow-up: [the full residual spectrum for every exterior type, the inverse from walks, and an exact completion tolerance](RESIDUAL.md).
