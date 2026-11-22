# What a Hamiltonian count remembers about thirteen independent flips

**PROVED.** Write \(E=2^{30}(2^{31}-1)\). The 40 odd integers from \(E-39\) to \(E+39\) freely generate a commutative monoid under multiplication inside the ordinary tournament spectrum. Its inner 26 generators come in thirteen single-arc pairs. Taking repeated ordinal sums of their strong witnesses gives arbitrarily large cubes in the labelled arc-reversal graph. On each of these cubes the Hamiltonian count determines exactly how many copies of each kind of arc have been reversed. [Interval and physical pairs](README.md)

## A private divisor from a narrow interval

Let distinct integers \(n_1,\ldots,n_s>1\) lie in an interval of integral width \(W\ge1\). From each \(n_i\) remove all prime powers supported at primes at most \(W\), leaving \(r_i\). If every \(r_i>1\), then the \(n_i\) are multiplicatively independent. Indeed any prime dividing \(r_i\) exceeds the nonzero differences \(|n_i-n_j|\) for \(j\ne i\), so \(r_i\) is coprime to every other generator and to \(n_i/r_i\). In a product \(P=\prod n_i^{k_i}\), the largest \(a\ge0\) such that \(r_i^a\mid P\) equals \(k_i\). This gives a decoder using only division by the residual integers. The condition is sufficient and makes no primality requirement on a residual.

For our 40 generators use \(W=78\). The smallest remaining divisor after the removal is 29,861,325,211, at offset \(-3\). This is a finite division test on 40 integers. Each generator already has a tournament, and ordinal sums multiply their Hamiltonian-path counts. Thus every nonnegative exponent vector yields its own distinct scalar count. The zero vector is realized by a singleton. [Divisions and band certificates](scalar_band.py)

## The exact quotient by repeated choices

For \(-6\le j\le6\) choose a strong tournament pair \(T_{j,0},T_{j,1}\) differing on one arc, with counts \(a_j=E+4j-1\) and \(b_j=E+4j+1\). Fix any multiplicities \(m_j\ge0\) of positive total \(M\). In an ordinal list of blocks, use \(m_j\) positions of type \(j\). Choose either orientation separately at every position, keeping each block's vertex labels. The resulting \(2^M\) tournaments form an induced, isometric \(M\)-cube in the labelled arc-reversal graph. Two outputs differ on precisely the chosen internal arcs corresponding to their differing bits. The Hamming distance therefore equals the number of reversals required. All outputs are distinct even up to isomorphism, since the ordered strongly connected components recover the blocks and their different path counts recover the bits.

If \(k_j\) positions of type \(j\) use the upper orientation, then
\[
      H=\prod_{j=-6}^{6}a_j^{m_j-k_j}b_j^{k_j}.
\]
The private-divisor lemma recovers each \(k_j\) from this one integer. Exactly \(\prod_j (m_j+1)\) heights occur. The fiber at the vector \((k_j)\) consists of exactly \(\prod_j {m_j\choose k_j}\) tournaments in the cube. It is the orbit of the permutations of the choice positions within each type. These permutations act on the cube of choices. The ordinal ordering of the components still distinguishes the corresponding tournaments up to isomorphism.

Collapsing equal heights *in this cube* gives the Cartesian product of the thirteen chains \(0,\ldots,m_j\), since a single coordinate reversal changes just one \(k_j\) by one. For one copy of each type, the 8,192 cube vertices have distinct scalar heights. For \(m\) copies of each type, the same description gives \((m+1)^{13}\) heights from \(2^{13m}\) tournaments. The increment in height along a flip of type \(j\) equals twice the product of the counts of all the other blocks. This gives a precise scale change for a microscopic gap of two. The construction supplies no bound on gaps between unrelated spectrum values. For more than one block these ordinal sums have several strong components.
