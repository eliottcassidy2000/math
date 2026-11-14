# The entire ordinal/outer-wrap closure is sparse

Start from singletons and all strong tournaments with the fewest
directed triangles at their order (the first triangle terrace). Close
under nonempty ordinal sums and the operation \(W(B)\) that surrounds
any nonempty middle \(B\) by two singleton ports and reverses their
outer arc. There is no bound on the number of operations, nesting
depth, or size of a middle.

**Result.** The entire closure has \(O(X^{0.952})\) distinct ordinary
Hamiltonian path counts at most \(X\); the strong outputs that are
outer wraps have only \(O(X^{0.928})\) distinct counts. The latter
includes arbitrary prior nesting. Both inverse-power series converge,
so the corresponding bounds can also be stated with little-o. Both
prime and seven-times-prime candidate spectra remain incomplete.

This sharpens the background researcher's sufficient exponent 0.975
from an ordered grammar: the key is replacing sequences by multisets
*at every nested node*, not just at the final ordinal product. It
uses [the color transfer and one-wrap count](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/WRAP_CHANNELS.md)
and [the underlying terrace pressure bound](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/PRESSURE.md).

## 1. Commutativity of full observations, not of the tournaments

Let \(F_m(B)\) be the sum over \(m\)-colorings of the product of the
Hamiltonian counts on the color classes, with empty classes counting
one. Then \(F_1=H\). Every \(F_m\) is multiplicative across ordinal
components. The whole profile \((F_m)_{m\ge1}\) updates under a wrap by
\[
 F_m(W(B))=m^2F_m(B)+m(F_{m+1}(B)-F_{m-1}(B)).
\]
Thus one may reorder children anywhere in a nesting expression
without changing any \(F_m\) at the top. No isomorphism between the
reordered tournaments is asserted. It would not suffice to know
commutativity of \(F_1,F_2\) alone, since subsequent wraps ask for
\(F_3,F_4,\ldots\).

Use atomic descriptions: either a base tournament or \(W\) of a
nonempty **multiset** of atoms. Their formal size is the actual
vertex count. Describe a general middle by a nonempty multiset of
atoms. Distinct descriptions may have equal profiles, harmless for
an upper bound.

## 2. A small description series at the color scale

The base series by order, counting all base classes, is
\[
 A(z)=z+z^3+\frac{z^4}{1-2z}.
\]
Indeed there is one singleton, one directed triangle, and
\(2^{n-4}\) normalized anchor words at every order \(n\ge4\).
Take \(s=116/125=0.928\) and \(z=3^{-s}<181/500\).
Let \(U_d(z)\) sum \(z^{|T|}\) over atoms of depth at most \(d\).
If \(U_d(z)\le1\), then their multiset series (including the empty
one) is at most
\[
 \exp\!\left(\frac{U_d(z)}{1-z}\right)<5.
\]
To see it, take the product of \((1-z^{|T|})^{-1}\) over atomic
descriptions and use \(-\log(1-x)\le x/(1-z)\) for \(0\le x\le z\).
The strict bound uses \(1/(1-z)<8/5\) and \(e^{8/5}<5\).
Therefore \(U_{d+1}(z)\le A(z)+4z^2<1\); the last rational inequality
already holds at \(181/500\). The induction starts at \(A(z)<1\).
Monotone convergence gives a finite atomic sum at unlimited depth,
and in particular a finite weighted sum over all reservoir
multisets. The integer replay verifies the small margins without
decimal assumptions.

## 3. Why the color cost pays for the grammar twice

For a wrapped atom,
\[
 F_2(W(B))=4F_2(B)+2F_3(B)-2F_1(B)\ge 2F_3(B)\ge2\cdot3^{|B|}.
\]
The first inequality uses \(F_2(B)\ge F_1(B)\), by coloring all
vertices the same; the second uses that each tournament has at least
one Hamiltonian path, even in the empty-class convention. Hence
the sum of \(F_2(\text{wrapped atom})^{-s}\) over all wrapped atomic
descriptions is finite by the previous section.

For base terrace atoms the same \(F_2^{-s}\) sum is finite. The
[earlier scanner proof](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/WRAP_CHANNELS.md)
greedily extracts disjoint three-site windows of type 1?1 in the
free anchor bits, giving a factor \((3/2)\) in \(F_2/2^n\) for each.
At the slightly sharper exponent \(s=0.928\), use the same seven
states \((\epsilon,0,1,00,01,10,11)\), with hit-weight
\((2/3)^s<6865/10000\). The positive test column
\[
 (10000,10000,9015,10000,9015,8869,8274)^t
\]
has one-step bound \(1.902\) times itself, strictly below \(2^s\).
This gives geometric convergence in base order. The singleton and
the order-three case add only finite terms.

Now all atoms have summable \(F_2^{-s}\) cost and \(F_2\ge2\).
Use a **second** multiset product, this time with weight
\(F_2(\text{atom})^{-s}\), for all reservoir profiles. Since
\(H(W(B))=H(B)+F_2(B)\ge F_2(B)\), it gives the asserted
summability of \(H^{-s}\) over distinct wrapped outputs.
Finally the unwrapped base counts \(H>1\) have summable
inverse powers at 0.952 from the terrace pressure argument.
One more multiset product on the **distinct count values** covers
all final ordinal products; singleton factors 1 can be ignored.
Convergence of \(\sum H^{-s}\) over distinct outputs implies
\(o(X^s)\) by discarding any fixed finite portion of the series
and bounding the tail above \(X_0\), then sending \(X_0\) to infinity.

In particular the closure cannot supply all but finitely many odd
primes, or all but finitely many \(7p\), since the sum of reciprocals
of the primes diverges. Unbounded linear memory is real, but it
does not by itself deliver coverage. A useful next construction must
escape the commuting outer-wrap grammar, not just add further
wraps or rearrange ordinal blocks.
