# One balanced target detects both open Chan–Pak requirements

For a tournament let \(H\) be its number of directed Hamilton paths. **PROVED.** Chan–Pak is equivalent to the following single family of positive tests. For every odd prime \(p\ne3,7\), some tournament with exactly three nonsingleton strongly connected components has height \(21p^2\). This is an open existence requirement. At every such \(p\) the target belongs to the same proper-divisor family \(F=U+S\).

The inherited definition uses \(D(N)=\{d\mid N:1<d<N\}\), with \(F=|D(N)|\), \(U\) its number of primes and \(S\) its number of squarefree integers. The prime divisors are also squarefree. For \(N\ge2\) equality holds exactly at \(p,p^3,p^2qr\), with three distinct prime factors in the mixed case. The inherited odd-square bracket theorem singles out the primes \(2,3,11\) by a different test, keeping a prime and its double between the same consecutive odd squares. In particular the forbidden tournament count 7 is a divisor-balanced prime, while 21 has \((F,U,S)=(2,2,2)\). The link below uses the three *choices of the squared prime on fixed support*, where multiplication of heights has a literal tournament interpretation. [Divisor classification](../trijection_reset_20261020/arithmetic.md) · [Odd-square brackets](../procgen_numerology_20260923_odd_square_brackets.md) · [THM-1370, the two excluded heights](../../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) · [Earlier Berggren selector on balanced divisors](../chanpak_berggren_divisor_twins_20261123/README.md)

The divisor-balance set also contains **no pair** \(N,\ell N\) for any prime \(\ell\). Its possible total numbers of prime factors, counted with multiplicity, are \(1,3,4\). A single prime multiplication could therefore only take the cube case to the mixed case, and the exponent three would persist. In particular the central move \(N\mapsto3N\) preserves the tournament spectrum by an ordinal triangle while always leaving divisor balance when started there.

## The three choices on the support \(\{3,7,p\}\)

A tournament is the unique ordered ordinal sum of its strongly connected components, and their heights multiply. Each nonsingleton strong component has odd height greater than one. Let \(B_d(z)\) count ordinary strong isomorphism classes of height \(d>1\) by order. In particular \(B_3=z^3\), \(B_7=B_{21}=0\), and \(B_{49}\ne0\). For \(i\in\{3,7,p\}\) write \(W_i(z)\) for the polynomial counting singleton-free ordered words of **three** strong components of product \(3\cdot7\cdot p\cdot i\). Then

| Squared factor | Height | Forced multiset of component heights | Exact order polynomial |
|---|---:|---|---|
| \(3\) | \(63p\) | \(3,3,7p\) | \(W_3=3z^6 B_{7p}\) |
| \(7\) | \(147p\) | \(3,49,p\) | \(W_7=6z^3 B_{49}B_p\) |
| \(p\) | \(21p^2\) | \(3,p,7p\) | \(W_p=6z^3 B_pB_{7p}\) |

Four prime tokens distributed among three factors greater than one permit exactly one merger. A lone 7 and the pair \((3,7)\) are forbidden, which forces the contents in each row. Four components are impossible. All listed component orders are freely realizable by ordinal sums, and their positions survive ordinary isomorphism. The formula also gives
\[
                  W_3(z)W_7(z)=3 z^6 B_{49}(z)W_p(z).
\]
Thus the first two choices detect \(7p\) and \(p\) individually and the third detects their conjunction. Singleton components at arbitrary positions multiply each grade-three series by \((1-z)^{-4}\). The eventual third forward difference of the total ordinary class counts reads \(W_i(1)\), by the [fixed-count series](../chanpak_support_birth_grammar_20261122/FIXED-COUNT-SERIES.md).

A witness at height \(p\) or \(7p\) has only one nonsingleton component, since 7 is absent and \(p\) is prime. Consequently Chan–Pak implies the third row is nonempty for every eligible prime. Conversely, those rows force all the required prime and \(7p\) cores. To build any other odd height, split off every prime other than 7. An exponent of 7 at least two is a sum of twos and threes, supplied by the known heights 49 and 343. A lone 7 can be paired with a prime \(p\ne3,7\) using \(7p\), or with two factors 3 using the known height 63. The only leftover integers are the already excluded 7 and 21. This specializes the inherited two-lane endgame to one family of mixed divisor shapes. [Window and SCC-grade reduction](../chanpak_window_endgame_20261123/README.md) · [49, 63 and 343 witnesses](../chanpak_midpoint_ogg_chambers_20261118/README.md)

**FINITE-EXACT, \(p=11\).** Here \(B_{11}=z^5, B_{49}=7z^7+z^8, B_{77}=10z^7+6z^8\). This gives \(W_3=30z^{13}+18z^{14}\), \(W_7=42z^{15}+6z^{16}\) and \(W_{11}=60z^{15}+36z^{16}\). These are minimum orders and multiplicities *within grade three*. The replay uses the [canonical class rows through order eight](../../../session212/cycle-levels/MANIFEST.json), the [order-nine strong spectrum omitting 77 and beginning at 75](../strong_H_spectrum_m9_isoclass_monad_s6.out), and the [Busch bound](https://doi.org/10.37236/1141) \(H\ge5^{(n-1)/3}\), which exceeds 77 for \(n\ge10\). The ordered choices already explain the multiplicities 30, 42 and 60 at the first orders. A search for ordinary positivity of \(21p^2\) would discard the component information needed for this reduction.

## Coprime heights cannot be reached in one flip from these witnesses

**PROVED.** If a tournament \(T\) has at least three nonsingleton strong components, one a directed triangle, then every one-arc neighbor \(T'\) has \(\gcd(H(T),H(T'))>1\). This applies at every order to **every witness** in the three rows above. Both proposed children \(3H(T)\pm2\) are coprime to the odd height \(H(T)\). The central child \(3H(T)\) can always be realized by a one-arc flip after adjoining a transitive triple as an ordinal summand and turning it into a directed triangle.

For the proof, a changed arc inside one component leaves some other nontrivial factor unchanged. A changed arc from one component \(A\) to a later component \(C\) leaves the components outside the interval as ordinal factors in both completions. Coprimality would therefore force the interval to span every nontrivial component. The intervening ordinal sum \(M\) then contains at least one nonsingleton component. Suppress the exterior singleton factors and reverse \(u\to v\), where \(u\in A,v\in C\). No old Hamilton path uses that arc, since a forward path must visit \(M\). The exact gain is
\[
           H(T')-H(T)=s_A(u)\,Q_2(M)\,e_C(v).
\]
Here \(s_A(u)\) counts ordered two-path covers of \(A\), with the second path starting at \(u\) and the first allowed empty. The reflected \(e_C(v)\) requires the first path to end at \(v\). Either of the two slots for \(M\) may be empty. Cutting a new path at \(v\to u\) gives these three independent ordered covers: both sides of the cut traverse the blocks in ordinal order. Concatenation is the inverse. The unrestricted two-slot count \(Q_2\) multiplies across ordinal sums.

If the triangle is in \(M\), its two-slot factor is 12. If it is \(A\) or \(C\), its indicated endpoint cover count is 3: the distinguished path can have length one, two or three in exactly one way each. Thus the spanning gain and the old height share a factor of three. For edits with an unchanged nontrivial exterior factor, their common divisor need not be three. The number of nonsingleton components is essential: one reversal joining two consecutive triangles takes 9 to 17.

More generally, with \(g\ge3\) nonsingleton components, an edit spanning them all yields
\[
                   H(T')\ge(1+4^{g-2})H(T).
\]
Cutting any Hamilton path immediately before \(u\), or after \(v\), injects into the respective endpoint covers. Every Hamilton path in an interior block \(D\) yields \(|D|+1\) ordered covers by its cuts, all distinct since concatenation recovers the path. Every interior nonsingleton block has at least three vertices, and there are at least \(g-2\) of them. Multiplying the bounds proves the claim. A chain of \(g\) triangles with a flipped arc between the two outer blocks attains equality. This gives an obstruction to small coprime-height corrections depending only on the number of nonsingleton components.

## Exact edit cost of removing copies of an odd prime from the height

**PROVED.** Let \(p\) be 3, 7 or 11, and let \(R_p\) be the cyclic tournament on \(\mathbb Z/p\), directed from \(i\) to \(i+j\) for \(1\le j\le(p-1)/2\). Take any prefix tournament \(A\) with \(p\nmid H(A)\), or the empty prefix of height one. Start from the ordinal sum \(T=A\to R_p\to\cdots\to R_p\) with \(r\) copies. For each \(0\le s\le r\), the distance in the single-arc metagraph from \(T\) to a tournament with \(v_p(H)=s\) is exactly
\[
                         \left\lceil\frac{r-s}{2}\right\rceil .
\]
The formula holds both on labelled tournaments and on ordinary isomorphism classes. In particular a move to \(3H(T)+2\) or \(3H(T)-2\) would cost at least \(\lceil r/2\rceil\) flips whenever \(r\ge1\), regardless of the arithmetic adjacency of those integers.

For the lower bound, \(k\) edited arcs touch at most \(2k\) of these modules. The untouched ones retain independent cyclic rotations, even if their SCC boundaries have changed. Their rotation product group acts freely on the set of directed Hamilton paths: fixing an ordered path fixes every vertex. Hence \(p^{\max(0,r-2k)}\mid H(T')\). A sequence in the unlabelled metagraph lifts stepwise along relabellings to the same number of physical edits.

For equality it suffices that \(v_p(H(R_p))=1\), one internal arc edit gives a \(p\)-unit height, and one edit joining two consecutive copies gives a \(p\)-unit height. Replace pairs of copies in disjoint adjacent positions by the last operation and, for an odd number of removed factors, edit one remaining copy internally. The resulting ordinal block product has valuation exactly \(s\). These sufficient conditions apply to other prime-rotation modules whenever their two local moves are available.

**FINITE-EXACT local data.** Number the cyclic vertices starting at zero, use the internal pair \(0,1\) and the two zero vertices in the two-copy move. Classify paths starting at zero on a subset by even or odd subset size, weighting each by the Hamilton count of its complement. Reflection gives the same two weights for paths ending at zero. The [rank-two adjacent seam formula](../chanpak_whip_seam_20261115/BACK_CROSS.md) gives the two-copy height \(h^2+4EO\).

| \(p\) | Original \(h\) | \(E,O\) | Internal edit | Two-copy edit |
|---:|---:|---:|---:|---:|
|3|3|1, 2|1|17|
|7|175|117, 142|157|97,081|
|11|93,027|71,589, 80,046|88,537|31,575,675,105|

At \(p=3\) the constructed height is exactly \(H(A)3^s17^{\lfloor(r-s)/2\rfloor}\). At 7, the cyclic module has height 175 even though 7 itself is absent. At 11, the other odd bracket prime produces a second genuine instance of the same edit law. The prime 2 differs intrinsically: every Hamilton height is odd, and a tournament admits no automorphism of order two, since interchanging a pair would reverse its mutual arc.

Run `python3 -B divisor_core_words.py` in this directory to reproduce the finite polynomials, all possible coprime exits from the singleton-free grade-three cores at \(p=11\), and the local valuation certificates.
