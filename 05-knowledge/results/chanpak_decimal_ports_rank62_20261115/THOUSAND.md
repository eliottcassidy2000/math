# What the decimal rungs preserve, and what their old arithmetic did not

Write $L_k=1000^k+1$ for $k\ge1$, with the transitive Hamilton count $1$ as a separate seed. It is not $L_0$, which would be 2. The exact arithmetic relation does not itself give a graph operation.

## The historical decimal triple

[Session 82](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session82_decimal_aliquot_square_sum_20261028/DECIMAL_TRIPLE.md) used
$$
B=t+1,\quad A=t^2+1,\quad C=t(t+1),\qquad A+B=C+2,\quad B+C=B^2 .
$$
At $t=1000$ its third label was $C=1{,}001{,}000$, not 1. The square-sum edge $B+C=B^2$ is genuine in that old graph. It is not an edge between decimal rungs: those are all 1 mod 4, so sums of two of them are 2 mod 4. The label $C$ is even, so it is not an ordinary Hamilton count in a tournament. The shared idea worth carrying forward is the old [inverse-aliquot fibre repair](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session82_decimal_aliquot_square_sum_20261028/FIBER_PORTS.md): two individually possible steps do not compose unless the same *named* port appears in both.

## Same physical core, different addresses

Use the directed triangle with outgoing rows $(4,1,2)$ and only vertex 0 marked, in the two-tape construction of [THM-5192](https://github.com/eliottcassidy2000/math/blob/main/01-canon/theorems/THM-5192-carrier-forms-the-insertion-identity-monotone-unit-cuts-and-odd-cycle-crowding-at-7-and-21.md). The count is
$$
12XY+6X+3Y+3[p(X)+p(Y)]+5+2[X\text{ odd}]+[Y\text{ odd}].
$$
Here $p(0)=0,\ p(2u+1)=p(u),\ p(2u)=p(u)+u+(u\bmod2)$. The literal endpoint guards make each output strong.

| Target | $X$ | $Y$ | Strong order |
|---:|---:|---:|---:|
| 1,001 | 4 | 18 | 12 |
| 1,000,001 | 1 | 61,481 | 21 |
| 1,000,000,001 | 804 | 103,594 | 31 |

The last realization is the one in the [preceding target-compiler session](../chanpak_cone_collision_placement_20261115/README.md), not a new range claim. The old [23-vertex certificate for 1,000,001](https://github.com/eliottcassidy2000/math/blob/main/session199/physical-adjoint/TARGET_1000001.json) is a different realization. The useful new comparison is that all three positive rungs use the exact *same* strong carrier. The addresses themselves give no radix recurrence.

There are also distinct ordinal (reducible) realizations. $1{,}000{,}001=101\cdot9901$; the factor 101 has order-seven rows $(98,48,107,67,13,88,18)$, and 9901 is a literal order-four-core carrier with out-neighbor rows $(8,1,3,6)$, mark mask 1 and addresses $(178,1)$, of order 14. For the billion rung, $1{,}000{,}000{,}001=1001\cdot999001$, and the latter factor has a 21-vertex host with the same four-vertex rows, mark mask 8, and addresses $(1500,20)$. Ordinal products therefore give order 21 and 33 witnesses, respectively. The identity is a component product, not proof of a compatible marked interface.

## Arithmetic clocks, not component births

For even $t$ and $g=\gcd(u,v)$,
$$
\gcd(t^u+1,t^v+1)=
\begin{cases}t^g+1,&u/g,\ v/g\text{ both odd},\\1,&\text{otherwise}.\end{cases}
$$
Modulo an odd common divisor, $a=t^g$ is a unit and $a^{u/g}=a^{v/g}=-1$. Bézout makes $a=\pm1$; an even quotient forces the common divisor to divide 2, and both odd force $a=-1$. The reverse divisibility in the odd case is immediate. Therefore adjacent positive rungs are coprime, so a nontrivial earlier strong component cannot be kept as an ordinal factor of the next rung. The subfamily $L_{2^j}$ is pairwise coprime; a finite fixed palette of nontrivial SCC counts cannot supply it by ordinal products.

| Prime(s) | Order of 1000 modulo the prime | Visible in exponent $k$ | Valuation there |
|---|---:|---|---|
| 7, 11, 13 | 2 | $k$ odd | $1+v_p(k)$ |
| 101, 9901 | 4 | $k\equiv2\bmod4$ | $1+v_p(k/2)$ |
| 19, 52579 | 6 | $k\equiv3\bmod6$ | $1+v_p(k/3)$ |

These follow from direct multiplicative orders and ordinary odd-prime exponent lifting. The first and second rows never coexist; the first and third coexist when $k\equiv3\bmod6$. This says which factors are arithmetically available, not which are exposed as tournament components.

The inherited [fixed-port spectral result](https://github.com/eliottcassidy2000/math/blob/main/session197/fixed-mixed/PROOF.md) localizes any periodic fixed-port approach: its candidate bases are $d^r e^\ell$ with $d,e\le k+1$ for $k$ fixed ports. Hence three ports cannot supply radix 1000; at four ports the factorization requires bases 2 and 5 and three births in each direction. A one-sided fixed-port period requires at least nine ports because $d^r=1000$ forces $d=10$ or $1000$. This is a necessary test only. The stronger [unmarked growing-module obstruction](DECIMAL_PORTS.md) explains why even a larger fixed acyclic substitution need not realize the coefficient one.

## The billion rung's aliquot orbit

The older orbits met at 1. The next positive rung again enters at 1, not by a shared earlier source:
```
1000000001 → 413350399 → 107622305 → 39974047 → 929673 → 439227
→ 212853 → 70955 → 18037 → 1079 → 97 → 1 → 0
```
Its only positive square-sum edge is $929673\to439227$, with sum $1170^2$; their gcd is 9, and the normalized edge still sums to $390^2$. For comparison the 1001 orbit has $343\to57$ with sum $20^2$, while the 1,000,001 orbit has none. The new edge has $929673=9\cdot53\cdot1949$; for prime $p\notin\{3,53\}$ the more general sum is $\sigma(9\cdot53\cdot p)=702(p+1)$, square exactly if $p=78z^2-1$, with $z=5$ here. These are separate arithmetic graphs, not a recursion of tournament realizers.

Run `python verify_thousand.py` for the literal marked-core arithmetic, independent direct subset controls, factor clocks, and aliquot trace.
