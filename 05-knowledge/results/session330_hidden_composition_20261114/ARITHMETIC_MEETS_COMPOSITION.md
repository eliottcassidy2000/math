# 5,723 and 3,225: the same arithmetic, different executable moves

## Factorizations are not the same as strong seams

\[
5723=59\cdot97,\qquad 3225=3\cdot5^2\cdot43=25\cdot129=43\cdot75.
\]
Both can be made by the same true Fermat block (five Hamiltonian paths, \(j=2\)) on an order-eight strong source at root 7. The actual formula is \(H+4L\), not \(5H\):
\[
5723=283+4\cdot1360,\qquad3225=133+4\cdot773.
\]
The source rows are \((6,40,82,21,35,77,155,63)\) and \((2,8,67,149,39,15,59,119)\). Both outputs are strong of order eleven. In the single-seam grammar with \(j\in\{1,2,4\}\), the first possible strong source orders are eight for 5723 and seven for 3225. An order-seven 3225 source (root 6) has rows \((2,20,65,7,45,15,59)\), giving \(73+16\cdot197\). That latter output has *twelve* vertices: minimizing input order need not minimize output order.

For 5723, a realization that genuinely factors into multiple SCC counts greater than one requires one count 59 and one count 97, so it takes at least fourteen vertices: both prime factors are first realized at order seven, by the small catalogue. Singleton padding of a single realizer is not such a factorization. In contrast, 3225 has genuinely different same-size 6+7 ordinal splits: profiles with counts 25 and 129, or 43 and 75. Their literal masks are checked by the replay. These restrictions do not prove global minimal order for a strong 5723 tournament.

## Arithmetic geometry of the two-sided cut

If \(U,W\) have vectors \((K_U,H_U),(K_W,H_W)\), the triangle output is the real part of
\[
(K_U+H_U\sqrt{-2})(K_W-H_W\sqrt{-2})=N+D\sqrt{-2}.
\]
Here \(D=H_UK_W-K_UH_W\) is, up to direction, the exact change from moving one neutral singleton between the reservoirs. In particular,
\[
N^2+2D^2=(K_U^2+2H_U^2)(K_W^2+2H_W^2).
\]
This is an exact dot-product plus oriented-area description of the seam. It provides a filter or coordinate system, not a physical source selector.

There is a separate arithmetic connection for the two requested integers: both are norms in this ring, each with two positive-coordinate representations:
\[
5723=75^2+2\cdot7^2=45^2+2\cdot43^2,\qquad
3225=55^2+2\cdot10^2=5^2+2\cdot40^2.
\]
These follow by choosing the conjugate of one factor: use \((3+5\sqrt{-2})(5\pm6\sqrt{-2})\) for the first, \(5(1+\sqrt{-2})(5\pm3\sqrt{-2})\) for the second. This does **not** identify those vectors with physical tournament profiles. It is the same algebra acting on a different input. The 97 in the factorization is also the vertex of the old square-sum exchange, but the norm identity alone provides no reason that specific vertex should form a spanning repair.

## Where 3,225 came from in the old work

[THM-454](../../../01-canon/theorems/THM-454-blowup-product-and-lift-laws.md) already realized 3225 by doubling every vertex into a transitive pair. It was one of two outputs, 3225 and 2785, from base tournaments with the same length-typed odd-cycle inventory: three triangles and one 5-cycle. In fact the two entire original slot polynomials are both $z^5+6z^3+2z$, but their four-cycle counts are three and two, and twin inflation can turn that difference into an odd cycle. This is an earlier example of hidden composition becoming visible one operation later, as in the [new slot observer](HIDDEN_COMPOSITION.md).

**A new uniform obstruction.** If every base vertex is inflated to a transitive module of order divisible by \(2^t\), then the entire optional-slot polynomial satisfies \(Q(z)\equiv z^n\pmod{2^{t+1}}\) coefficientwise, for output order \(n\); in particular \(H\equiv1\pmod{2^{t+1}}\). For each \(k>0\), count ordered packings of \(k\) disjoint odd directed cycles by occupancies \(u_{ic}\) in module \(i\) and cycle \(c\). Given these occupancies, independent choices of disjoint subsets within a transitive module contribute a row multinomial; order-preserving relabeling makes the remaining directed-cycle count independent of the chosen subsets. Any odd cycle has some odd occupancy \(u_{ic}\). Its row multinomial has a factor \(\binom{m_i}{u_{ic}}\) divisible by \(2^t\), since
\[
u_{ic}\binom{m_i}{u_{ic}}=m_i\binom{m_i-1}{u_{ic}-1}.
\]
Keep the total number \(m\) of cycle vertices fixed along with \(k\), so each group contributes to the same exponent \(n-m+k\). Thus the ordered-packings count is divisible by \(2^t\). Odd-cycle packing weights mean one divides it by \(k!\) and multiplies by \(2^k\). Since \(v_2(2^k)-v_2(k!)\ge1\), every nonempty contribution is divisible by \(2^{t+1}\), and the empty packing contributes one. The first-level bound is sharp: a doubled directed triangle has 45 paths.

So 5723 is excluded from uniform even transitive inflation at every order (it is 3 mod 4), while 3225 is a literal twin-inflation witness and cannot be uniformly inflated by module sizes divisible by 8, since \(v_2(3225-1)=3\). Both still occur by the same Fermat block. This separates two seemingly interchangeable ways of adding a geometric piece.
