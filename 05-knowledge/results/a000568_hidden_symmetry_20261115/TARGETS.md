# What the symmetry census can rule out at an exact odd target

## A complete order-five Burnside slice for $3\nmid h<30{,}035$

If $H(T)=h<30{,}035$ and $3\nmid h$, an order-five automorphism is possible **only** in the following strata, and each occurs with exactly the indicated number of strong cores. Every core has ten vertices and ordinary automorphism group $C_5$.

| $h$ | 8,245 | 15,395 | 15,505 | 15,565 | 15,745 |
|---:|---:|---:|---:|---:|---:|
| strong cores | 1 | 2 | 1 | 1 | 1 |

Every other vertex is a singleton SCC. Hence at order $n$ multiply the core inventory by $\max(n-9,0)$ (place the core anywhere in the ordered word of singletons). These are exact class counts for this nonidentity slice, not the class counts for the whole $H=h$ fiber. In particular every realization of $5^a$ for $1\le a\le6$ is rigid. That extends the $5^3$ part of the inherited [prime-cube theorem](https://github.com/eliottcassidy2000/math/blob/main/session197/orbit-wildcard/PRIME_CUBES.md); it does not prove those realizations exist.

**Proof.** Choose an order-five automorphism. Each moved orbit of five vertices is regular, with $H=15,c_2=90$. If all its interactions with the other moving orbits are uniform, it is a module, since fixed vertices also see it uniformly. Its whole nonempty ordered-cover vector is divisible by three: the first two coordinates are $15,90$ and $r!\mid c_r$ handles $r\ge3$. Thus for $3\nmid h$ there must be a mixed pair. Index both by translations in $\mathbf F_5$. A mixed pair has 4 skew masks for each internal tournament and 30 nonconstant difference masks, for 480 tables. Independent endpoint and startpoint subset recurrences yield the complete histogram:
```
8175:80  8245:40  8265:40  15395:80  15465:120
15505:40  15565:40  15745:40
```
A third moving orbit would already force $H\ge8175\cdot15$ by the disjoint path-merge inequality (the earlier orbit-spectrum packet). A fixed vertex distinguishes the two orbits only by seeing their entire blocks oppositely. Across 960 such augmented tables both DP recurrences find the minimum $30{,}035$. Thus below that cutoff the mixed pair's union $U$ is a module, and it is strong.

If the quotient containing the contracted $U$ is strong and nonsingleton, there is a directed triangle containing its marked vertex. Its lift gives $H\ge3H(U)+c_2(U)\ge12H(U)$, as each of the nine cuts of every spanning path injects into ordered nonempty two-covers (concatenation in the original order recovers the path and cut). If $U$ itself is an SCC, any other nontrivial SCC would multiply the count by at least five, since the total is prime to three. In either case any nontrivial enlargement exceeds the cutoff. The only remaining growth is singleton padding.

An order-25 subgroup of $S_{10}$ is Sylow and contains single supported 5-cycles. Such a subgroup would expose independent 5-rotations in the two blocks, incompatible with a mixed mask. Therefore the diagonal $C_5$ in each table is Sylow in the automorphism group. By Sylow conjugacy an isomorphism between two fixed tables can be adjusted in the destination to normalize this chosen action. Its normalizer consists of independent translations, a common nonzero multiplier, and possible exchange of the two blocks. These operations give the class table above in the companion replay. Finally, no factor three is allowed by $3\nmid h$; the only potential additional odd prime in $S_{10}$ is seven. Since no factor 25 is possible and $v_7(10!)=1$, that would give a group of order 35, which is cyclic, and a permutation of order 35 cannot fit in ten vertices. So the group is exactly $C_5$.

## Divisor 7: coupling, and a physical change from $7q$ to a multiple of 63

For a prime cyclic automorphism with moving orbits indexed by $\mathbf F_p$, join two orbits in an interaction graph when their cross-mask is mixed. A blockwise translation vector $(b_i)$ is an automorphism exactly when $b_i=b_j$ along each interaction edge: a nonempty proper subset of $\mathbf F_p$ has no nonzero translation period. Fixed-vertex incidences and uniform masks impose no additional conditions. This subgroup is $C_p^{\#\text{components}}$. Consequently at $H=7q$ with $q\ge11$ prime the interaction graph must be connected. The residual Burnside search in [MODULES](MODULES.md) can reject disconnected tables before counting their paths.

One intact rotating block is different. The cyclic interval and Paley 7-circulants have, respectively, first coordinate $175,189$ and common tail
$$
(c_2,\dots,c_7)=(1638,6678,15120,20160,15120,5040).
$$
These inherited profiles are recorded in [uniform prime transport](../uniform_prime_transport_20261103/THEORY.md). Dividing by seven, the tail entries have gcd $18$. In any unmarked host quotient $Q$, the coefficient of $c_1$ is exactly $H(Q)$. Hence a replacement from interval to Paley changes the actual count by $14H(Q)$ while leaving every multi-run contribution unchanged. The interval case has $H/7\equiv 25H(Q)\bmod18$, so a prescribed $7q$ needs $H(Q)\equiv13q\bmod18$ and $H(Q)\le q/25$. The Paley case has $H/7\equiv27H(Q)\bmod18$ and is therefore always divisible by **63** in an unmarked host. The singleton and triangle hosts give $189,2205$, with gcd 63, so this is a sharp contextual divisor. For instance the [inherited first symmetric $7q$ witness](https://github.com/eliottcassidy2000/math/blob/main/session196/symmetry/PROOF.md) gives the paired outputs $4151=7\cdot593$ and $4221=63\cdot67$. The [inherited diagonal example](https://github.com/eliottcassidy2000/math/blob/main/session197/orbit-wildcard/PROOF.md), recovered up to cross-mask relabeling in the replay, has two coupled Paley blocks but is prime to three. Thus the structure producing 63 is independent module action, not just the presence of copies of Paley 7.

These filters identify which symmetries a target can support. The remaining Chan–Pak obligation is to produce a rigid candidate, or an admissible coupled candidate, at each absolute target.
