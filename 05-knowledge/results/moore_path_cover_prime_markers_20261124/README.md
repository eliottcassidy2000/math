# Prime factors of ordered path covers: Moore degrees 2, 3, 7 and 57

There are three different kinds of primes here: factors of graph matrices, factors forced by a symmetry of a specific tournament, and factors that appear when its path-cover polynomial is evaluated. They need not agree. All degree-57 encodings below are **conditional on a Moore graph**; the polynomial and cell results themselves are unconditional.

## More counts for the existing charts

Write $Z_k(T)$ for the number of ordered $k$-path covers, with empty paths allowed. Using the [Grinberg–Stanley formula](https://arxiv.org/abs/2307.05569) color an entire directed odd cycle with one color and otherwise color vertices independently. Thus, with $n=|T|$,
$$
 Z_k(T)=\sum_{\mathcal C}2^{|\mathcal C|} k^{n-|V(\mathcal C)|+|\mathcal C|}
       =k^e F(k^2).
$$
Here $\mathcal C$ ranges over disjoint directed odd-cycle packs, and $e$ is the smallest nonzero exponent. This is the older repo's [odd-cycle conflict enumerator](../../../06-writeups/research-document-for-mathematicians.md) with an additional grading by vertices covered. At $k=1$ that grading disappears.

| Chart | $Z_1$ factored | $Z_2$ factored | $Z_3$ factored | $\gcd_{k\ge1}Z_k$ |
|---|---|---|---|---:|
| C5 cyclic | $3\cdot5$ | $2^3\cdot3\cdot5$ | $3\cdot5^2\cdot7$ | 15 |
| C5 paired 01432 | $7\cdot13\cdot71$ | $2^5\cdot7\cdot601$ | $3^2\cdot37\cdot4153$ | 1 |
| Petersen sparse spine | $13\cdot401$ | $2^3\cdot31\cdot421$ | $3^2\cdot7\cdot37\cdot463$ | 1 |
| Petersen compact | $3\cdot5\cdot19\cdot29$ | $2^5\cdot3\cdot5\cdot17\cdot19$ | $3^2\cdot5\cdot7\cdot13\cdot367$ | 15 |
| Petersen hidden-star | $3\cdot7\cdot4003\cdot9923$ | $2^4\cdot3\cdot533017853$ | $3^2\cdot115811\cdot380231$ | 3 |

[profile.py](profile.py) enumerates each odd-cycle count by subset/endpoint DP, assembles the closed polynomial, and prints the exact values and residual prime factors (using GNU factor) through $k=6$. The displayed fixed gcds also follow from finite differences: values at $0,\ldots,n$ determine the gcd of all integer evaluations, since the forward-difference coefficients in the binomial basis are integers. No full $Z_k$ for the compact 50-vertex Hoffman–Singleton chart is claimed here. There is a guaranteed factor 5 for that chart at **every** $k$: global residue translation simultaneously in all its five-site cohorts is a cyclic automorphism, and every nonidentity automorphism acts without fixed points on the ordered list of directed path words. This also proves the forced 5 for the cyclic C5 and compact Petersen charts. Neither sparse Petersen nor the paired C5 chart has that forced common factor.

A useful sieve separates forced and incidental factors:
- Parity: $Z_k\equiv k^n\bmod2$, so 2 can never divide every evaluation. Also $Z_k\equiv k^n Z_1\bmod3$ for all $k$, because all exponents have the same parity. Thus 3 divides every evaluation exactly when it divides $Z_1$.
- For odd $p$, $p$ divides every evaluation iff $x^{(p-1)/2}-1$ divides $F(x)$ over $\mathbb F_p$: the nonzero square values of $k^2$ are precisely the roots of that polynomial, and the $k=0$ case contributes zero to $Z_k$. Thus $p\le n-e+1$ for a fixed odd prime. To be forced *after* stripping $k^e$, the extra condition $F(0)\equiv0\bmod p$ is needed. For example 5 divides the full compact Petersen readings at every $k$, but not every normalized reading.
- Outside of fixed factors, a nonconstant $F(k^2)$ has unbounded prime factors even when $k$ is confined to any arithmetic progression. Indeed apply the elementary Schur argument to the resulting integer polynomial $G(j)$ after picking $G(0)\ne0$: if a fixed finite set contains all its primes, evaluate at $j$ a large multiple of $G(0)$ times their product. Then $G(j)/G(0)$ is 1 modulo every prime in the set but has absolute value greater than one. The fact that $k$ can grow here is not sufficient for the earlier fixed-$k$ unit problem.

## A height-preserving Moore code with actual prime markers

In the [original faithful edit chart](../moore57_tournament_observers_20261124/moore-tournament-representations.md), two cells of height 9 encoded each bit. **They have the same entire profile**, not just the same height:
$$
 Z_k=k^5+6k^3+2k .
$$
Therefore even knowing every scalar $Z_k$ gives no information about the literal bit word in that deck. Its constant-height encoding is a stronger example of the earlier [Golay example](../chanpak_h23_golay_chains_20261124/README.md) where the ordinary count collapses many different contact placements.

A replacement pair exists on six vertices. Take the common spine $0,1,\ldots,5$ and reverse the nonspine pairs corresponding to masks 40 and 42 in lexicographic pair order, least significant bit first. These strong tournaments differ in just one free arc, have $H=37$, and have profiles
$$
 P_0(k)=k^2(k^4+12k^2+24),\qquad P_1(k)=k^2(k^4+14k^2+22).
$$
The exhaustive fixed-spine search in [cells.py](cells.py) finds no strong same-height, different-profile, one-free-arc pairs through order five and 122 labeled pairs at order six. The 9-cells above do not contradict this because their profiles agree. At $k=2$ the new cells give $352=2^5\cdot11$ and $376=2^3\cdot47$. At $k=3$ they give $1917=3^3\cdot71$ and $2061=3^2\cdot229$. This is a literal way for the numbers $2,3,11$ from the earlier modular-selection discussions to meet: evaluate at 2 or 3; 11 is a factor of **one particular bit**, not a universal geometric prime.

In fact the residuals $P_0(k)/k^2,P_1(k)/k^2$ have no common odd prime other than 37. Their difference is $2(k^2-1)$; any common odd prime forces $k^2\equiv1$, then divides 37. It occurs exactly for $k\equiv\pm1\bmod37$, with common valuation one (otherwise the residual difference would force $k^2\equiv1\bmod37^2$). The common 2-part is exactly 2 if $k$ is even, and 1 otherwise; substitute $k^2=4u$ to see the two residual valuations 3 and 1 in the even case. At the marker primes 11 and 47, $k=2$ is a simple root of the appropriate residual with derivative 80 or 88, respectively. Hensel lifting plus CRT lets one prescribe any two positive exact exponents at these separate primes by choosing $k$, without either prime dividing the other cell's residual. This is **a readout CRT**, unlike [the strip CRT](../chanpak_tournament_sl2_20261124/STRIP-CRT.md), which programs residues of actual heights. Neither supplies a general absolute target inverse.

For a Moore graph use the same literal closed-neighborhood code as before. The bit string at $v$ is $1_{N[v]}$ with $d+1$ ones among $m=d^2+1$ sites. Replace bits with the six-site cells and take ordered joins. All heights are $37^m$; specifically $v_{11}(Z_2)=m-d-1$, $v_{47}(Z_2)=d+1$, and $Z_3$ has the same exponents on 71 and 229. To recover distance rather than just degree, build the same join on the XOR of two codewords: then $v_{47}(Z_2)$ is exactly their Hamming distance. Adjacent and nonadjacent vertices have neighborhood intersection 2 or 1, respectively.

| Moore degree | Exponents on 11 and 47, single row | $v_{47}$ on paired XOR, adjacent / nonadjacent |
|---:|---:|---:|
| 2 | 2, 3 | 2 / 4 |
| 3 | 6, 4 | 4 / 6 |
| 7 | 42, 8 | 12 / 14 |
| 57 (conditional) | 3192, 58 | 112 / 114 |

The numbers 11, 47, 71 and 229 are chosen by the tournament cells. All rows have equal weight and equal unmarked counts even in the new code. The XOR is a prescribed physical pairing operation, not information recoverable from the scalar of either row by itself. More generally, introducing 47 here does not get around the [fixed-prime obstruction](../moore_chanpak_fixed_prime_rays_20261124/): a fixed deck still supplies only finitely many primes at fixed $k$.

## The graph-derived primes are different

Let $A$ be the Moore graph adjacency and $X=A+I$. On the centered space, $A^2+A=(d-1)I$, and $(2A+I)^2=(4d-3)I$. On the constants, the eigenvalues of $A,X$ are $d,d+1$. The centered eigenvalues and their multiplicities follow from the quadratic and trace zero. This gives:

| Degree | Centered eigenvalues of $A$ | $\lvert\det A\rvert$ | $\lvert\det X\rvert$ | Discriminant $4d-3$ |
|---:|---|---|---|---:|
| 2 | $(-1\pm\sqrt5)/2$, twice each | $2$ | $3$ | 5 |
| 3 | $1^{[5]},(-2)^{[4]}$ | $3\cdot2^4$ | $2^7$ | $3^2$ |
| 7 | $2^{[28]},(-3)^{[21]}$ | $7\cdot2^{28}3^{21}$ | $2^{24}3^{28}$ | $5^2$ |
| 57 | $7^{[1729]},(-8)^{[1520]}$ | $3\cdot19\cdot2^{4560}7^{1729}$ | $2^{5188}7^{1520}29$ | $3^2 5^2$ |

So at degree 57, **19** comes from the constant mode of $A$, **29** from that of $A+I$, and **7** in the latter from a centered eigenspace. **3 and 5** instead collide the centered quadratic roots while $A+I$ stays invertible modulo either. This separates the [prior 5-primary Moore Smith calculation](../moore57_switching_flags_20261124/LOCAL-SMITH.md) from both literal bit primes and tournament counts. For $X=A+I$ at degree 57 the ranks modulo $2,7,29$ are respectively 1521, 1730 and 3249. For 7 and 29 split off the constants (3250 is a unit) and use the two distinct centered eigenvalues of $X$; only the $-7$ eigenspace or the constant eigenvalue 58 vanishes, respectively. Modulo 2 the characteristic polynomial of $X$ has a 1520-dimensional 1-primary space and a 1730-dimensional 0-primary space. Since $X(X+I)=J$ and $XJ=0$, $X$ has rank one on the 0-primary space, and is invertible on the other. These are matrix primes, not forced divisors of $Z_k$ for a hypothetical host.

## A way to make new primes at **fixed** path-cover width

Replace one vertex of a directed triangle by a transitive $a$-chain, leaving the other two singleton. Every directed odd cycle contains both singletons and an odd-size subset of the chain, traversed in chain order; no two such cycles are disjoint. Therefore the full profile is
$$
 Z_k=k^{a+2}+k\bigl((k+1)^a-(k-1)^a\bigr).
$$
At $a=3$ this is *exactly* the old height-nine cell polynomial. Its first three readings grow as follows:

| Chain $a$, not Moore degree | $Z_1$ | $Z_2$ | $Z_3$ |
|---:|---|---|---|
| 2 | 5 | $2^5$ | $3^2\cdot13$ |
| 3 | $3^2$ | $2^2\cdot3\cdot7$ | $3\cdot137$ |
| 7 | $3\cdot43$ | $2^2\cdot3\cdot11\cdot37$ | $3\cdot22817$ |
| 9 | $3^3\cdot19$ | $2^2\cdot3\cdot7\cdot17\cdot29$ | $3\cdot47\cdot6823$ |

There is an exact bridge to the [old 57 atlas](../number_bridges_57_27_16_20261028/ATLAS.md): the height at chain length 9 divided by the height at length 3 is $(2^9+1)/(2^3+1)=\Phi_{18}(2)=57$. The order of 2 modulo 19 is 18, so 19 divides the height precisely for lengths $9\bmod18$, whereas 3 already occurs at length 1. This makes 19 a new factor of the cyclotomic step, independently of its appearance in the Moore constant eigenvalue. Along lengths $a=2^j$, the height is exactly the Fermat family in the [Pépin/Collatz session](../s58_fermat_pepin_collatz_20261026/README.md). Pairwise coprimality gives new factors whether or not there are more Fermat primes.

More generally, **fix any positive evaluation $k$ and any finite prime alphabet $S$**. Only finitely many lengths $a\ge1$ have $Z_k$ supported on $S$. To prove it, divide the displayed three-term equality by the output $Z_k$ and enlarge $S$ by the primes of the nonzero bases $k-1,k,k+1$; omit the zero term for $k=1$. This is a nondegenerate S-unit equation: neither positive term is equal to the absolute negative term because $a\ge1$, and terms of one sign cannot cancel. [Evertse–Schlickewei–Schmidt](https://annals.math.princeton.edu/2002/155-3/p04) gives finitely many normalized tuples. The ratio between two terms of distinct positive bases determines $a$. Hence *simultaneously* at $k=1,2,3$, every sufficiently large length has a factor outside any prescribed finite alphabet in each of the three readings. This provides a concrete source of primes that the fixed-deck substitution lacks. It does not show that two such variable fillers balance to a unit, or that they hit prescribed Chan–Pak heights.
