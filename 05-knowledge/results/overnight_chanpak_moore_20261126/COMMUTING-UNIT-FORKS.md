# Two vertex-disjoint reversals give exact commuting unit forks

2026-11-26. This passes beyond the [fixed binary left-arm obstruction](ISOLATED-SPINE-REVERSALS.md): one nonbinary six-vertex core suffices. A unit means a change of two in the Hamiltonian-path count. The construction is exact and all four tournaments are strongly connected. It does not cover every prescribed odd integer.

## 1. A physical square with heights $N,N-2,N+2,N$

Use left root 1 on vertices 0 through 5, with outgoing masks $(6,60,56,49,33,1)$. Reverse $3\to0$. The new masks are $(14,60,56,48,33,1)$. The marked profiles $(H,E,R)$ are $(23,8,24)$ and $(23,10,24)$. Here $E$ counts spanning paths starting at the root and $R$ counts ordered two-path covers whose first path ends at the root, allowing the second path to be empty. The replay counts all three quantities by complementary-subset convolution and checks connectivity.

The right block is the converse binary core with address $c$ and right root $r$. Its dual profile is $(F(c-1),[c\text{ odd}],c)$. All cross arcs go right to left except the marked root bridge left-to-right, so the [seam](TRANSFERRED-BOUNDARIES.md) has count

    H F(c−1) + R c − E[c odd].

For odd $c$, the marked left edit lowers the count by exactly two. For the right edit, choose $i\ge1$ satisfying $i\equiv1\pmod{11}$, let $\ell=2^i$, and choose **any even suffix** $0\le r_0<\ell$. Write $Z$ for its $i$-bit reverse-complement and $\epsilon=Z\bmod 2$. Put

    β = (24ℓ−2)/23 + Z + ε,
    q = β or β−1,
    c = 2ℓ q + r₀ + 1.

The order of two modulo 23 is 11. Thus $\beta$ is integral and even: the numerator $24\ell-2$ is even and the divisor odd, and $Z+\epsilon$ is even. The selected bit $i$ of $c-1$ is zero, and $q+[q\text{ odd}]=\beta$. Turning on that bit reverses the right-root spoke, a physical edge disjoint from the marked left arc. The [suffix formula](README.md) gives

    ΔF = Z+ε−β,
    ΔN = 24ℓ+23(Z+ε−β)=2.

The new right address remains odd, so the seam's marked endpoint multiplier stays one. Consequently the four orientations have counts $(N,N-2,N+2,N)$: their mixed difference is zero. Both left orientations are strongly connected. The right root can reach every site via its highest 1-bit and its decreasing spine. The flipped bit is below that highest site, and every right vertex has a cross arc back to the left. Thus all four orientations are strongly connected.

With $i=1,r_0=0$, the replay counts **$(557,555,559,557)$** at $c=13$ on 11 vertices and **$(791,789,793,791)$** at $c=17$ on 12 vertices by direct full-tournament subset DP. The all-order construction is not inferred from these checks.

## 2. Many heights and the remaining inverse

On each valid scale $\ell$, there are $\ell$ pairs $(r_0,q)$, as the suffix ranges over even values and $q$ over the two choices. Both $q$ and $Z$ are $O(\ell)$ and $q=\Omega(\ell)$, so all heights on that scale are $\Theta(\ell^2)$. The [square-root window](README.md) applies with the blocks interchanged by global converse: for a fixed right-edit midpoint the value $Z$ lies in an interval of at most 26 consecutive integers at this scale, after the constant $H$ and $R$ are fixed. For each such $Z$ there are at most two prefixes. Hence at least $\ell/52$ heights on that scale are distinct. Summing or selecting across the geometric scales $i\equiv1\pmod{11}$ proves that this explicit family produces **$\Theta(\sqrt X)$ distinct central heights up to $X$**. The upper bound here is only for this particular family, not for the spectrum.

The arithmetic guard on the address-9 family is genuinely broken. At $i=23,r_0=0,q=\beta$ the central height is 13,516,917,716,131,391 and the lower midpoint is 1 modulo 9. Further, at $i=12$ the 4,096 parameter choices have distinct heights; they hit all residues modulo each of 81, 63, 99, 231 and 135 (separate modulus tests, not arbitrary simultaneous powers). Neither fact hits the exponent-61 perfect midpoint itself. The earlier exact inverse already excluded all rooted cores through order seven at exponent 61 for this binary-spoke lane, including these order-six cores. The next target is to change the core profile or allow a second nonbinary side.

## 3. Two more physical cores

For any isolated left change $\Delta E=2$ with $(H,R)$ fixed, require $R\ell\equiv\sigma\pmod H$ for $\sigma\in\{-2,2\}$. Then $\beta=(R\ell-\sigma)/H+Z+\epsilon$ gives a right shift $\sigma$. Set the starting left orientation so its own edit has shift $-\sigma$. The same even-suffix construction works when $q>0$.

| Left masks | Root | Left arc endpoints | $(H,R,E_{\rm low})$ | Right indices |
|:--|--:|:--|:--|:--|
| $6,60,56,49,33,1$ | 1 | $0,3$ | $(23,24,8)$ | $1\bmod11$: $+2$ |
| $12,53,56,50,33,1$ | 1 | $2,4$ | $(25,28,8)$ | $4\bmod20$: $-2$; $14\bmod20$: $+2$ |
| $12,45,56,48,35,1$ | 1 | $3,5$ | $(29,31,9)$ | $14\bmod28$: $-2$; $0\bmod28$: $+2$ |

Each core's two orientations are checked literally; the replay also checks instance recurrences through bit index 90. This construction explains why the universal binary-spine obstruction is not an obstruction to independent arc behavior: a nonbinary profile can preserve both ordinary channels while varying the endpoint channel by two, and its coprime coefficients let the other arm solve the signed unit equation.

## 4. A direct cover update for the unique-terminal grammar

Given a nonempty rooted core $(C,u)$ on the left and a converse binary block $A$ on the right, make the same seam and mark its right root $r$. Write $V_C$ for the two-cover count with second path starting at $u$ and $Q_C(2)$ for all ordered two-slot covers with empties allowed. At the new root the terminal count is unchanged, and the new two-cover count is

    e_G=e_C,    V_G=c Q_C(2)−V_C.

For $e_G$, a spanning path ending at $r$ has the unique ordered prefix $A-r$ before a spanning path of $C$ ending at $u$, then takes $u\to r$. For $V_G$, neither of the two cover paths can use $u\to r$: the second path already starts at $r$, and the first avoids $r$. Both therefore travel from an $A$-run into a $C$-run. Naively this gives $V_AQ_C(2)=cQ_C(2)$. The sole forbidden join is $r\to u$. It requires the second $A$-run to be the singleton $r$, making the first $A$-run the unique path of transitive $A-r$, and the second $C$-run to start at $u$. This subtracts exactly $V_C$. The replay verifies the formula on 192 literal orientations. It supplies an explicit two-cover update for the [Fibonacci grammar](UNIQUE-TERMINAL-GRAMMAR.md).

Replay with `node commuting_unit_forks.js`.
