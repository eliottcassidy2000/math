# The six locally universal rigid blocks hit just one square shifted down by two

**PROVED.** Let \(\mathcal B=\langle5,17,37,77,129,751\rangle\) be the multiplicative monoid of heights obtained by ordinal joins of the six rigid tournaments, allowing a singleton and the empty product. For positive integers \(s\),
\[
  s^2-2\in\mathcal B\quad\Longleftrightarrow\quad
  s=113,\qquad 113^2-2=12\,767=17\cdot751 .
\]
The [six-block construction](../chanpak_tiling_viewpoint_20261123/PRECEDENCE-FIBER.md) supplies a 15-vertex tournament at that exact height. The two possible orders of its strong blocks give different ordinary classes; both have the same transitive pair-precedence tournament up to relabeling. This classifies the entire six-block family, with unbounded total order. Outside this family, other shifted squares occur, starting with \(5^2-2=23\).

The contrast with modular selection has the same quantifiers at every square. For every fixed odd \(s\ge3\) and every \(a,b,c\ge0\), some word of the six blocks gives \(H\equiv s^2-2\pmod{2^a3^b11^c}\). Of these targets, only the displayed one is an equality. A shift by two is the coordinate change: a missing height 7 has shifted label 9, and missing height 21 has shifted label 23. The result concerns the whole quadratic family of labels.

## Four blocks disappear before any Pell search

Every odd prime dividing \(s^2-2\) admits 2 as a quadratic residue. A direct reduction modulo \(5,37,11,3\), respectively, rules out factors \(5,37,77,129\): these primes are \(3\) or \(5\bmod8\). Thus a candidate is \(17^a751^b\). Since \(s\) is odd and \(17\equiv1,751\equiv7\bmod8\), \(b\) is odd. Write
\[
 (D,y)=
 \begin{cases}
 (751,17^{a/2}751^{(b-1)/2})& a\text{ even},\\
 (12767,17^{(a-1)/2}751^{(b-1)/2})& a\text{ odd}.
 \end{cases}
\]
The exact condition is \(s^2-Dy^2=2\), with \(y\) supported on \(17,751\).

## The odd-power description of the Pell solutions

For either of these \(D\), let \((x_0,y_0)\) be its positive norm-two solution with least \(y_0\). Set \(\alpha=x_0+y_0\sqrt D\) and \(\epsilon=\alpha^2/2=(x_0^2-1)+x_0y_0\sqrt D\), an integral norm-one unit. Every positive solution has
\[
 s+y\sqrt D=\alpha^n/2^{(n-1)/2},\qquad n\text{ positive and odd}.
\]
Indeed \(\epsilon>1\) and \(\epsilon=\alpha/(2/\alpha)\), so multiplication by an integral power of \(\epsilon\) takes any positive norm-two element to \((2/\alpha,\alpha]\). For its coordinates \((x',y')\), use \(2\sqrt D\,y'=z-2/z\) in that interval. This gives \(-y_0< y'\le y_0\). The norm equation excludes zero. Minimality of \(y_0\) forces the reduced element to be \(\alpha\). Undoing the reduction gives the formula and positivity of \(y\) forces \(n\ge1\). In particular every such \(y\) is a multiple of the odd part of \(y_0\), by expansion of the odd power.

For \(D=751\) the smallest solution is
\[
 (x_0,y_0)=(2700614460969,98546821297),\qquad 271\mid y_0 .
\]
Here is a short finite certificate of minimality. In lowest terms (the norm forces coprimality) each positive solution satisfies \(0<x/y-\sqrt{751}<1/(2y^2)\), and is consequently a continued-fraction convergent by Legendre's criterion. The norms of the successive convergents through the displayed one are
\[
 -22,21,-6,37,-15,33,-14,9,-39,10,-31,21,-30,13,
 -27,26,-15,5,-35,18,-7,25,-3,42,-11,2 .
\]
Their denominators increase. The replay regenerates the list by the integer square-root algorithm. Since 271 is coprime to both allowed primes, no solution survives this branch. [Legendre’s convergent criterion, Proposition 1.6](https://www.imo.universite-paris-saclay.fr/~arnaud.durand/files/topics_mnt.pdf)

For \(D=12767\) the least solution is \((113,1)\). Write \(n=2k+1\) and \(y_n\) for the coefficient of \(\sqrt D\) in the normalized odd power. At each \(p=17,751\),
\[
 v_p(y_n)=v_p(n).
\]
Indeed in the binomial expansion of \(2^k y_n\) the leading term \(n113^{n-1}\) has this valuation. For each later term, indexed by \(1\le j\le k\), use
\[
 v_p\!\left(\tbinom{n}{2j+1}113^{n-2j-1}D^j\right)
 \ge v_p(n)-v_p(2j+1)+j>v_p(n).
\]
The strict inequality holds since \(p^j>2j+1\) for \(p\ge5\). If \(y_n\) has only the two allowed prime factors, this proves \(y_n\mid n\). For \(k\ge1\), its first positive binomial term alone gives \(y_n\ge n(113^2/2)^k>n\). Thus \(n=1\), forcing \(y=1\), \(s=113\), and \(a=b=1\).

The argument gives a reusable test for other banks. First discard an entire block whenever one of its prime divisors forbids the square residue of the shift. For the shift two, each remaining exponent-parity pattern gives one norm-two equation. A least-solution denominator containing a forbidden odd prime excludes the whole pattern. If every allowed odd prime divides \(D\) and is at least five, valuation growth at those primes is bounded by \(v_p(n)+v_p(y_0)\); positive growth of the real solution then eliminates all later odd powers whenever \(x_0^2>2\). The primes with global modular selection and the primes dividing the denominator coordinate play distinct roles.

Run `python3 -B bank.py` in this directory for the finite certificate and the exact seed construction.

## Adjoining every allowable small height leaves only two squares

**PROVED.** Let \(\mathcal M\) be generated by all attained heights in \([1,25]\). The known low spectrum gives precisely the positive odd numbers there except 7 and 21. For positive integers \(s\),
\[
 \begin{array}{c|c}
 \text{allowed generators} & \{s: s^2-2\text{ occurs as a product}\}\\ \hline
 \mathcal M & \{5\}\\
 \mathcal B & \{113\}\\
 \langle\mathcal M,\mathcal B\rangle & \{5,113\}.
 \end{array}
\]
Both entrances are literal: \(5^2-2=23\) and \(113^2-2=17\cdot751\). Thus mixing the two inventories in any order or multiplicity creates no further shifted odd square. The [ordinary height exclusions, THM-1370](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) concern the integers 7 and 21. The displayed intersection restricts only these ordinal constructions, whose path counts multiply.

Among the new positive heights through 25, the only nonunit candidates passing the prime-divisor test are 17 and 23. There are now three possible prime factors \(17,23,751\). As \(23,751\equiv7\pmod8\), exactly one of the last two exponents is odd. Removing a square leaves precisely four possible \(D\), namely \(23,391,751,12767\), and now \(y\) may have any of the three prime factors. For \(D=391\), the least norm-two solution is \((2709,137)\), whose denominator 137 excludes the branch. The convergent norms through that solution are \(-30,9,-15,13,-19,18,-15,2\). For \(D=751\), the previously found divisor 271 still excludes every solution.

For either remaining \(D\), \(y_0=1\). At an allowed prime \(p\nmid D\), work in \(\mathbf F_p[t]/(t^2-D)\), allowing the split case. Since \(p\) is odd, \(t\) is a unit and the reduced element \(s+y_nt\) equals its conjugate exactly when \(p\mid y_n\). The ratio of the two conjugates of the normalized odd power equals \(\epsilon^n\). Thus \(p\mid y_n\) exactly when the multiplicative order of \(\epsilon\) in this finite ring divides \(n\). The three needed orders follow by repeated multiplication.

| \(D\) | \(p\) | Order of \(\epsilon\bmod p\) | An exclusion whenever \(p\mid y_n\) |
|---:|---:|---:|---|
| 23 | 17 | 9 | \(3\mid n\), so \(7\mid y_3=49\) divides \(y_n\) |
| 23 | 751 | 125 | \(5\mid n\), so \(y_5=2351\), coprime to all three primes, divides \(y_n\) |
| 12767 | 23 | 11 | \(857\mid y_{11}\) divides \(y_n\), and 857 is coprime to the three primes |

The divisibility in this table holds for odd \(r\mid n\): expand the normalized power at \(n\) as the odd power \(n/r\) of the solution at \(r\) divided by a power of 2. Each \(y_r\) here is odd, so \(y_r\mid y_n\). Consequently a surviving solution uses only primes dividing \(D\). The ramified valuation inequality above gives \(y_n\mid n\); the first positive binomial term then forces \(n=1\) and \(s=x_0\). This gives exactly the displayed roots. The two smaller monoids inherit the upper bound. The first has no factor 751 and the second has no factor 23, which settles the separate classifications.

## Minimum physical flips at both surviving targets

**FINITE-EXACT.** Each surviving shifted square and its square are also the two Hamilton counts of an actual single-arc flip. The rows below use zero-based outgoing masks for the lower endpoint. Both flips have the minimum possible number of vertices even for representing just the lower count.

| \(s\) | Lower height | Upper height | Minimum order | Lower rows | Pair flipped |
|---:|---:|---:|---:|---|---|
| 5 | 23 | 25 | 6 | 2, 16, 11, 35, 13, 23 | 0,4 |
| 113 | 12,767 | 12,769 | 10 | 594, 360, 539, 209, 898, 349, 916, 807, 525, 42 | 4,9 |

The direct subset recurrence in the replay gives both counts before and after changing the indicated bit pair. Maximum heights at orders five and nine are 15 and 3357 respectively; a common singleton pads every smaller order. At 12,767 the 15-vertex ordinal witness and the lower endpoint of this 10-vertex seam give different representatives of the same integer. The two-component ordinal witness itself cannot be the lower endpoint of a unit flip between heights at least five, by the strong active-core criterion. [Maximum-height census](../path_spectrum_two_angles_metagraph_20261106/README.md) · [Active-core theorem](../chanpak_arc_fiber_grammar_20261122/README.md)


At the minimum orders, the complete marked census reveals how many representatives support the actual move. A marked edge records the lower tournament and its editable unordered pair up to simultaneous relabeling.

| Flip | All lower classes | Eligible lower classes | All upper classes | Reached upper classes | Marked edges | Edges after converse merging |
|---|---:|---:|---:|---:|---:|---:|
| 23 to 25 | 4 | 4 | 2 | 2 | 4 | 2 |
| 12,767 to 12,769 | 502 | 24 | 489 | 22 | 24 | 12 |

Every participating endpoint here is distinct from its converse. Each eligible lower class has exactly one marked edge up to isomorphism; each is rigid because its height has no prime divisor at most its order. All four ordinary order-six lower classes work at the smaller shifted square, whereas fewer than one in twenty of the order-ten lower classes work at the larger one.

Run `python3 -B bank.py --census` for the complete marked counts, using nauty's `gentourng` and a C++17 compiler. To cover any marked flip, delete one of its endpoints and enumerate every way of restoring it to the parent. The replay uses all 12 ordinary parents of order five and all 191,536 of order nine, testing every contact mask and every incident flip. It links the [all-apex height oracle](../chanpak_arc_fiber_grammar_20261122/birth_cut_cover.cpp) and directly recounts both sides of each retained edge. For the quotient, canonically refine vertex colors by outdegrees to color cells, then minimize over all permutations left within ties; repeat with the editable endpoints initially distinguished as an unordered set. Converse merging takes the smaller of the two marked canonical keys. The census sizes of the ambient integer fibers are from [the order-by-height atlas](../path_spectrum_two_angles_metagraph_20261106/data); sum the ordinary-class multiplicities grouped by cycle level.

## Which missing small generator changes the square rows?

**PROVED.** If all odd labels \(1,3,\ldots,25\) are admitted as *formal* multiplicative generators, the exact set of roots is
\[
                 \{s:s^2-2\text{ is a product}\}=\{3,5,11,45,235\}.
\]
Restoring the single forbidden generator 7 to the actual small bank gives this whole five-element list. Restoring 21 alone gives only 5. The generator 21 has the factor 3, and \(s^2-2\) is never divisible by 3. In the coordinate of shifted odd squares the two missing generators therefore have sharply different global effects.

The prime-support test now leaves exactly \(7,17,23\); the exponents of 7 and 23 have opposite parity. This gives \(D=7,119,23,391\). The last case still has the forbidden divisor 137 in every possible \(y\). In the other three cases the least norm-two solution has \(y_0=1\), with \(x_0=3,11,5\) respectively. We use the same odd-power index \(n\) as above. At an unramified allowed prime \(p\), let \(r\) be the order of \(\epsilon\bmod p\). When \(r\mid n\), expanding the normalized \((n/r)\)-th power of \(x_r+y_r\sqrt D\) gives
\[
                    v_p(y_n)=v_p(y_r)+v_p(n/r).
\]
Indeed \(p\mid y_r\), \(p\nmid x_r\); each higher binomial term has at least \(2j-v_p(2j+1)>0\) additional factors of \(p\) over the leading valuation. If \(r\nmid n\), then \(v_p(y_n)=0\). Together with the ramified valuation formula this gives the following explicit bounds whenever \(y_n\) has only allowable prime factors.

| \(D\) | Unramified orders for allowed primes | Bound | Positive odd indices left by \(y_n\ge n(x_0^2/2)^{(n-1)/2}\) |
|---:|---|---:|---|
| 7 | \(r_{17}=3,\ v_{17}(y_3)=1;\ r_{23}=12\) | \(y_n\le17n\) | 1, 3 |
| 119 | \(r_{23}=11,\ v_{23}(y_{11})=1\) | \(y_n\le23n\) | 1 |
| 23 | \(r_7=3,\ v_7(y_3)=2;\ r_{17}=9,\ v_{17}(y_9)=1\) | \(y_n\le833n\) | 1, 3, 5 |

For example \(833=7^2\cdot17\); the factors \(p^{v_p(n)}\) multiply to at most \(n\). All order and valuation checks use integer arithmetic in the replay. At \(D=23,n=5\) the forbidden divisor 2351 removes that index. The surviving roots and exact factorizations are
\[
 3^2-2=7,\quad 5^2-2=23,\quad 11^2-2=7\cdot17,\quad
 45^2-2=7\cdot17^2,\quad 235^2-2=23\cdot7^4.
\]
The missing generator 7 blocks these four products in the small ordinal grammar. Only the integer 7 itself is excluded among their ordinary Hamilton counts.

In fact each of the other three lost products admits an exact physical single-arc flip to the square at its first possible order. This recovers the missing multiple of 7 without using a height-7 factor.

| \(s\) | Actual flip | Minimum order | Lower rows | Pair flipped |
|---:|---:|---:|---|---|
| 11 | 119 to 121 | 8 | 246, 124, 120, 241, 224, 192, 128, 6 | 4,7 |
| 45 | 2023 to 2025 | 9 | 476, 113, 66, 310, 420, 389, 56, 78, 198 | 6,8 |
| 235 | 55,223 to 55,225 | 11 | 1766, 500, 952, 691, 1601, 1872, 1164, 304, 601, 194, 910 | 0,10 |

The exact census omits 119 at order seven. The maximum at eight is 661 and the maximum at ten is 15,745, excluding the other two at smaller orders. Each pair above is checked by the same independent endpoint recurrence. [Order-by-height spectrum](../path_spectrum_two_angles_metagraph_20261106/README.md)
