# A single arc invisible through the first return

**PROVED.** Given any two nonempty tournaments \(A,B\), put two singletons \(u,v\) around them in the order
\[
 u\longrightarrow A\longrightarrow v\longrightarrow B\longrightarrow u,
 \qquad A\longrightarrow B .
\]
Let \(G_0\) complete the missing edge by \(u\to v\), and \(G_1\) by \(v\to u\). Both completions are strong, since each original vertex lies on a directed cycle following the four blocks and can reach all other blocks. This completion turns differences between the ordered path-cover counts of the modules into a literal marked-arc difference. For \(a_j=Q_j(A), b_j=Q_j(B)\), set
\[
 W_j=a_{j+1}b_j-a_jb_{j+1}, \qquad
 Q_j(G_1)-Q_j(G_0)=j(W_j+W_{j-1}),\quad j\ge1 .
\]
Here \(Q_j\) permits empty slots, \(Q_1=H\), and \(W_0=0\). In particular the signed Hamilton-count change is \(a_2b_1-a_1b_2\). The two scalar counts differ by two precisely when the integer vectors \((a_2/2,a_1)\) and \((b_2/2,b_1)\) have determinant of absolute value one. The numerators \(a_2,b_2\) are even by exchanging the two slots.

### The determinant and its sign

Only covers using the flipped arc contribute to the difference. Choose its slot in \(j\) ways and contract its two consecutive endpoints to one ordered atom. Write \(x,y\) for nonempty path runs in \(A,B\). In the new order \(v,u\), the portion preceding the atom is empty or finishes in \(A\), while the portion following it is empty or begins in \(A\). This gives \((1+x)(1+x+xy)\). In the old order \(u,v\) use \(B\) instead, giving \((1+y+xy)(1+y)\). Each other slot admits the four patterns \(1,x,y,xy\). The signed compressed polynomial divided by \(j\) is therefore
\[
 (x-y)\bigl(1+(1+x)(1+y)\bigr)((1+x)(1+y))^{j-1}.
\]
A monomial \(x^i y^k\) counts ordered covers of \(A,B\) by \(i,k\) nonempty runs. Consequently \((1+x)^s(1+y)^t\) evaluates to \(a_s b_t\). Substitution gives the stated sum of two adjacent determinants.

**Equality criterion.** For \(d\ge1\), the two completions agree at slots \(1,\ldots,d\) exactly when the module ratios \(a_j/b_j\) are constant for \(j=1,\ldots,d+1\). Indeed \(W_0=0\), so vanishing in successive slots is equivalent to \(W_1=\cdots=W_d=0\), and all \(b_j>0\). In particular all slots agree exactly when the two monic module polynomials coincide: an indefinitely constant ratio gives proportional polynomials, and their leading terms then force the ratio to be one. This identifies exactly how many consecutive observations can fail to see the single arc.

### A 32-vertex pair resolved only at depth two

Take the two order-15 modules that already agreed through three slots and differed at the fourth. Each row below specifies an ordinal sum of the two factors listed, with outgoing masks in vertex order. [Prior three-slot collision](../chanpak_induced_gap_three_slot_descent_20261118/README.md)

| Module | First factor | Second factor |
|---|---|---|
| \(A\) | \(8,9,3,100,15,23,55\) | \(144,5,161,7,78,27,47,122\) |
| \(B\) | \(32,97,19,7,11,28,61\) | \(112,97,67,7,14,28,56,127\) |

Their common reads through grade three are \(17325,5105184,317017071\), and \(a_4=8535609600,\ b_4=8540125440\). The determinant makes \(G_0,G_1\) agree at slots 1 and 2, with third-slot difference \(3a_3(a_4-b_4)\). For the triangle return \(R\), the known transfer gives \(H(RG)=Q_1+Q_2\) and \(H(R^2G)=-Q_1+5Q_2+2Q_3\). Thus the same physical flip is neutral on orders 32 and 34 and changes the count on order 36 by exactly \(-8,589,590,219,427,840\). [Prior return transfer and triangular inversion](../chanpak_bott_cayley_periods_20261115/CAYLEY.md)

| Number of returns | Count on the \(u\to v\) side | Count on the \(v\to u\) side |
|---:|---:|---:|
| 0 | 51,861,365,876,187 | 51,861,365,876,187 |
| 1 | 587,178,721,172,950,515 | 587,178,721,172,950,515 |
| 2 | 1,644,123,578,453,985,193,299 | 1,644,114,988,863,765,765,459 |

For the absolute values, enumerate each compressed slot in the four-block source, recording which singletons occur and how many nonempty \(A,B\) runs appear. Combine slots with disjoint singleton masks and retain the union containing both singletons. The number of fillings by \(i\) nonempty runs of \(A\) is \(\sum_{j=0}^i(-1)^{i-j}\binom ij a_j\), and similarly for \(B\). The checker compares this separate absolute calculation against literal small hosts and the signed determinant; it also checks the least one-return witnesses at orders 6 and 7. The order-32 example is an existence claim. [Exact replay](return_pair.py)

### A unit flip on a square with two arbitrary modules

**PROVED.** This four-block construction has an absolute gap of two if and only if the two module orders are 1 and 2, in either order. The counts in that case are 9 and 11.

At \(j=1\) the formula reduces to \(Q_2(A)H(B)-H(A)Q_2(B)\). In terms of the integer vectors \((H(A),Q_2(A)/2)\) and \((H(B),Q_2(B)/2)\), a unit half-gap is exactly the unimodular determinant condition for two Farey neighbors. For every module of order at least two, \(Q_2\) is divisible by four. Indeed the complementary subsets form an even number of pairs, each contributing twice an odd product. If both modules have order at least two, their second coordinates are even, which rules out unimodularity in all sizes. If both are singletons the determinant vanishes. With one singleton and an \(n\)-vertex module \(A\), the absolute gap is \(Q_2(A)-2H(A)\), the count of ordered covers with both slots nonempty. Every nonempty proper subset has at least one path, as does its complement. Thus the gap is at least \(2^n-2\), exceeding two once \(n\ge3\). For \(n=2\) the transitive pair gives equality. For instance the lower row masks \((14,28,24,16,1)\) and the flip of \((0,3)\) have counts 9 and 11.

In particular, take **any** exact perfect straddle at an odd exponent at least five, after stripping its common singleton ends. Both types of two-step path between the switched vertices must occur by the nested-contact theorem. If all other vertices belong to these two types, at least one arc in **each** direction must occur between the types. Otherwise they are exactly the two modules of the square just classified (after exchanging the endpoint names if necessary). Thus in the absence of a common in-neighbor or out-neighbor, the straddle must break the uniform orientation between the two flanks. [Common square and attachments](SQUARE-EARS.md)
