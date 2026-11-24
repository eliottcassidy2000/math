# Simultaneous selection inside one strongly connected tournament

**PROVED.** For each \(M=2^a3^b11^c\) with nonnegative exponents there is an \(n_0(M)\) such that **at each prescribed order \(n\ge n_0(M)\)** every odd-representable residue modulo \(M\) occurs as the Hamilton-path count of a strongly connected, rigid tournament whose odd-precedence relation is transitive. All constructions use the same six-vertex outer tournament and the six known height seeds. The [ordered-join selector](../chanpak_tiling_viewpoint_20261123/PRECEDENCE-FIBER.md) and the [independent strong construction using modulus-dependent caps](../chanpak_tiling_viewpoint_20261123/STRONG-CRT.md) give the inherited local theorem. Here the outer tournament has six vertices for every modulus, and an extra 3-adic coordinate lets one keep that same outer tournament. The integer output and an effective finite search for a prescribed absolute target within this family are given below.

## What survives substitution

Let \(D(T)\) direct \(i\to j\) if an odd number of Hamilton words put \(i\) before \(j\). The [pair-precedence theorem](../chanpak_tiling_viewpoint_20261123/ORDER-PARITY.md) obtains this from \(T\) by toggling all pairs in the vertex set of every simple directed cycle. Replace each vertex \(i\) of an outer tournament \(T\) by a nonempty module \(B_i\), orienting cross arcs as in \(T\).

**Substitution law.** Between two singleton blocks, take the direction from \(D(T)\). Between distinct blocks of which at least one is nonsingleton, take the original direction from \(T\). Inside a block of size at least three take \(D(B_i)\). Inside a two-vertex block, toggle its original arc precisely if an odd number of simple directed cycles of \(T\) contain \(i\). These four rules give the full labelled precedence tournament.

Here is the parity calculation. Fix an exterior cyclic word with at least one exterior vertex and indexed gaps for \(r\ge1\) nonempty visits to a module of size \(k\). Root at the least exterior label used to count each cycle once. Require a specified set of \(t\) module vertices to appear. The other vertices are optional. For prescribed subsets assigned to the runs there are an odd number of internal path fillings, by Rédei on each induced tournament. By inclusion-exclusion, the parity of the number of fillings is that of
\[
 \sum_{j=0}^r(-1)^{r-j}{r\choose j}j^t(j+1)^{k-t},
\]
with \(0^0=1\). If \(0<t<k\), every term is even. For either \(t=0\) or \(t=k\), the sum has parity \(2^{r-1}\), so exactly the single visit survives. For an exterior prescribed pair, applying the \(t=0\) rule to the other modules contracts them back to single vertices. One may root at a prescribed exterior vertex throughout. For a pair in different blocks with one nonsingleton use \(0<t=1<k\) in that block. For a pair inside one block retain all internal cycles, and use \(t=2\) for cycles leaving the block. When \(k=2\) its surviving single visit can root the contracted cycle while disposing of the other modules. These are exactly the four clauses.

In particular inflating *every* vertex into three ordered vertices produces a fixed point of \(D\) containing the original tournament as quotient. To make a strong nonsingleton outer tournament become transitive in precedence after inflation requires at least two singleton blocks; otherwise all its cross directions survive. Two suffice: for a cyclic triple with one block \(B\) of size at least three and two singletons, the arc between the singletons reverses and the result in precedence is singleton, \(D(B)\), singleton in order. This explains the structural role of the singleton boundaries.

## The six-vertex wrapper

Use outgoing masks \((2,4,9,35,15,23)\) and replace vertex 0 by \(B\), with \(|B|\ge3\). The outer source is strong. In the precedence tournament on the five *exterior* vertices their order is \(5,2,4,3,1\), computed in the full six-vertex source. Vertex 0 loses to the first four and beats the last. Consequently the expanded tournament has precedence order \(5,2,4,3,D(B),1\) whenever \(D(B)\) is transitive. It is strong: a directed outer cycle through the port permits travel between any two module members, and all other outer reachabilities lift. It is rigid because any automorphism preserves the unique precedence order.

Let \(h=H(B)\) and let \(q=Q_2(B)\) count covers by two *ordered* paths, allowing empty slots. A full word can visit the expanded vertex at most twice. Every return goes from its unique outneighbor 1 to an inneighbor, and 1 can occur only once. The six-vertex skeleton has 17 valid full words with the port visited once and six with it visited twice. Enumerate the latter by reading 0 twice nonconsecutively and every other letter once. Thus the exact output is
\[
                     17h+6(q-2h)=5h+6q.                    \tag{1}
\]
In an ordinal join, both \(h\) and \(q\) multiply: concatenate the component paths independently in each occupied slot, including empty component paths. The seed data, verified directly by the [replay](selector.py), are

| Seed \(i\) | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|
| \(d_i=|B_i|\) | 4 | 6 | 6 | 8 | 9 | 9 |
| \(h_i\) | 5 | 17 | 37 | 77 | 129 | 751 |
| \(q_i\) | 32 | 224 | 352 | 1904 | 4884 | 17636 |

For a word of these blocks with multiplicities \(e_i\) and \(t\) singleton letters the order of the outer tournament is \(n=5+t+\sum d_ie_i\). Write \(h=\prod h_i^{e_i}\). Formula (1) becomes
\[
 H_{\rm out}=5h+6\,2^t\prod q_i^{e_i}
 = h\left(5+6\,2^{n-5}\prod\lambda_i^{e_i}\right),
 \qquad \lambda_i=q_i/(2^{d_i}h_i).
\]
The first expression is integral; the last will be read at 3.

## The extra ternary coordinate

Each \(\lambda_i\) is a 3-adic unit congruent to 1 modulo three, cancelling the single factor 3 in both the numerator and denominator for the 129 seed. For any \(a\ge3,b\ge2,c\ge1\), the four unit seeds \(5,17,37,751\) generate the full group
\[
 (\mathbb Z/M)^\times\ \times\ (1+3\mathbb Z_3)/(1+3^b\mathbb Z_3)
\]
The second factor records principal 3-adic units modulo \(3^b\). Each seed is sent to \((h_i,\lambda_i)\). The first projection already generates by the [six-bank proof](../chanpak_tiling_viewpoint_20261123/PRECEDENCE-FIBER.md). The additional coordinate is a 3-group. The 3-primary part has two cyclic factors; to generate it, it suffices to span its quotient by cubes since a finite abelian 3-group equal to its own cubes is trivial. In coordinates \((h_i^2,\lambda_i)\bmod9\), the seeds 5 and 751 give \((7,4)\) and \((7,1)\), respectively. They span both copies of the order-three group \(\{1,4,7\}\). The other Sylow projections are unchanged. A subgroup surjecting on each primary part of a finite abelian group is the whole group. Smaller moduli follow by projection.

Now prescribe \(n\) and the residue \(r\). For the inner height request the CRT target
\[
 h\equiv r/5\pmod{2^a11^c},\qquad
 h\equiv r/(5+6\,2^{n-5})\pmod{3^b}.
\]
Both divisions use units. Prescribe the simultaneous second coordinate \(\prod\lambda_i^{e_i}=1\bmod3^b\). First include \(c\) copies of the 37-seed: its second-slot count 352 is divisible by 11 and its height is a unit at all three primes. Include enough 129 and 77 blocks to supply the truncated 3- and 11-valuations of the desired \(h\). The unit group acts transitively on every fixed valuation stratum of a prime-power ring; choose the correcting unit at each prime and combine by CRT. The preceding product-group theorem supplies that unit and the required correction of the ternary coordinate *independently*, by nonnegative powers of the other four seeds. All these choices depend on finitely many states at fixed \(M\), including the phase of \(2^{n-5}\) at 3. Choose a word in advance per state. At sufficiently large \(n\), padding inside the wrapper gives \(t=n-5-\sum d_ie_i\ge a\) and \(|B|\ge3\). The second term in (1) vanishes modulo \(2^a11^c\), while the fixed ternary coordinate supplies its prescribed value modulo \(3^b\). This proves the theorem at the *same* order for all residues, with both parities of that order allowed.

## Order bounds and the sharp residue wheel

For \(a\ge3,b\ge2,c\ge1\) the four correcting generators act on a group of order \(g=\varphi(M)3^{b-1}\). A shortest directed word to a reachable state has at most \(g-1\) letters: a repeated state gives a deletable segment. Each letter costs at most nine vertices. The preliminary blocks have total size at most \(6c+8c+9b\). Thus
\[
 n_0(M)\le 5+\max(3,a)+14c+9b+9(g-1).
\]
One can replace the last term by the maximum shortest *vertex cost* in the same group. At \(M=792=8\cdot9\cdot11\) the 720 states have maximum such cost 101, giving \(n_0(792)\le141\). These are uniform sufficient orders.

Modulo 24 the formula simplifies further. For \(|B|\ge2\), pair each subset of \(B\) with its complement in the definition of \(Q_2\). There are \(2^{|B|-1}\) pairs, and each contributes twice an odd integer. Hence \(4\mid q\), and (1) reads \(H_{\rm out}\equiv5h\pmod{24}\). The ordinary transitive-precedence fiber supplies every odd residue on eleven vertices, and adding singleton components supplies the same residues at every higher inner order. As 5 is a unit, the fixed shell converts those witnesses to strong rigid ones at **every order at least 16**. Using just the six fixed seed types plus singletons inside this shell has its own higher threshold of 27, obtained by shortest vertex-cost products modulo 24.

The [three-parent construction](STRONG11.md) realizes the twelve residues at eleven. The [twelve fixed marked tournaments](EXTENSIONS.md) supply every order from 11 by changing just one chain length. Smaller orders miss residue 11 even in the entire transitive-precedence fiber. Consequently eleven is also the first order from which the **strong** fiber covers every odd residue at *every subsequent order*. This sharp lower bound is specific to the precedence restriction.

## Exact integer boundary

For a given positive odd \(N\), enumerate exponent vectors with \(5\prod h_i^{e_i}+6\prod q_i^{e_i}\le N\). Such vectors lie in a finite box since each \(q_i\ge32\). For each vector the remaining quotient \((N-5\prod h_i^{e_i})/(6\prod q_i^{e_i})\) must be an integral power of two. Its exponent is the only possible padding. Keeping the module size condition makes this an exact integer inverse for the stated wrapper.

**FINITE-EXACT.** At either neighbor of the exponent-61 perfect number, 22,337 exponent vectors survive the size inequality and none returns that neighbor. Thus even strong connectivity, rigidity, one prescribed order at the modular stage and arbitrarily deep simultaneous local selection can coexist with a concrete failure of integer membership in the same family. The wider square and general Chan–Pak constructions remain available beyond it.

## Many ordinary isomorphism classes at a single exact height

**PROVED.** For fixed \(M\), at every sufficiently large order \(n\) and for each admissible residue \(r\), there is a *single integer* \(h_{n,r}\equiv r\pmod M\) represented at that order by at least
\[
 \frac{[z^L]\,(1-z-z^4-2z^6-z^9)^{-1}}{(L+1)^4}
\]
different strong rigid ordinary isomorphism classes with transitive precedence. Here \(L=n-N_M\), with one constant \(N_M\) for the modulus. These are equalities of integer heights, strengthening the previous [exponential class bound per residue](../chanpak_tiling_viewpoint_20261123/STRONG-CRT.md). The minimum possible order for a given height is a separate condition.

Use the inner alphabet consisting of a singleton and the four unit-height seeds 5, 17, 37 and 751. The coefficient counts its ordered words of total size \(L\). There are at most \((L+1)^4\) choices of the four nonsingleton multiplicities; together with \(L\) they determine the singleton multiplicity too. Choose a largest group of words having identical multiplicities. Each word in it has *exactly the same* \((h,q)\), since both reads multiply. Its height is a unit at 2, 3 and 11, and its normalized second read is a principal unit at 3. The correction argument above therefore works after this prefix: absorb its first and second reads into the two independently prescribed unit targets. The required corrections and final padding fit in a fixed budget \(N_M\), uniformly in the prefix and in the phase of \(n\). In the extended exponent regime one can use the sufficient order bound in the preceding section as this budget.

Use the identical correction word and trailing padding for the whole group. The wrapper formula gives exactly the same integer height for all of them. Its unique precedence order identifies the entire prefix vertex interval. On that interval the ordered list of strongly connected components recovers every letter, including the two differently weighted six-vertex components. Hence distinct words give distinct ordinary isomorphism classes even after wrapping. A small example on 15 vertices uses the two words of seed heights 5 and 17 in opposite orders. Both give height \(5\cdot85+6\cdot7168=43433\).

The coefficient grows at least as a positive constant times \(\rho^L\), where \(1=\rho^{-1}+\rho^{-4}+2\rho^{-6}+\rho^{-9}\), \(\rho\approx1.53850\). Indeed the coefficient recurrence and this equality induct from the first nine positive coefficients. Also \(Q_2(B)\ge2^{|B|}\), since every subset and its complement each admit at least one Hamilton word. Thus these high-multiplicity integer heights themselves tend to infinity with the prescribed order.
