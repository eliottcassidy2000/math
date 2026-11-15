# What can and cannot amplify an exact carrier bank

These are statements about operations on exact supports, not statements that the tournament spectrum or any full carrier is sparse.

## Multiplicative closure preserves reciprocal summability

For a set A of distinct attained counts >1, let S be its multiplicative monoid. If \(\sum_{a\in A}a^{-s}<\infty\) for some s>0, then
$$
\sum_{m\in S}m^{-s}\le\prod_{a\in A}(1-a^{-s})^{-1}<\infty .
$$
The product counts factor multisets and can count one integer multiple times. Convergence follows by taking logarithms, using \(-\log(1-t)=O(t)\) near zero. Splitting the sum at a fixed D and then taking D to infinity gives \(\#(S\cap[1,X])=o(X^s)\).

Consequently a finite union of exact supports with **upper bounds** \(O(X^\theta(1+\log X)^C)\), θ<1, even together with any fixed finite census, stays density zero at every depth of ordinal-product joining. A power-sparse **lower bound** for a carrier does not qualify it for this obstruction.

## Even long exact windows and local factors do not suffice

For any finite bank of odd witnesses there exists a multiplicatively closed set of odds containing that bank and arbitrarily long intervals of consecutive odds, with relative density exactly 1/2 in every fixed compatible residue class. This is an arithmetic comparison set, not a model of the tournament spectrum.

Start with Liouville's completely multiplicative sign λ and flip the signs to +1 at a set R of primes. Initially include the prime factors of the bank. At stage k choose fresh odd primes \(p_0,\ldots,p_{k-1}>2k\) that divide no previous window, with their reciprocal sum at most \(2^{-k}\). By CRT choose odd N with \(N+2j\equiv p_j\pmod{p_j^2}\) for every slot j. Each p_j occurs to valuation one in its own slot and in no other. Flip p_j only if its slot would otherwise have negative sign. Future flips do not alter this window. The sign-positive odds are multiplicatively closed.

Ordinary Liouville has sum o(X) on every fixed progression: use \(\lambda(n)=\sum_{d^2\mid n}\mu(n/d^2)\), the [Möbius prime number theorem for progressions](https://doi.org/10.1002/cpa.3160020208) on each fixed d, and the tail bound \(X\sum_{d>D}d^{-2}\). Coprime residues suffice because dividing by the gcd of the residue and modulus reduces the rest by complete multiplicativity. For finitely many changes of prime signs, expand each \((-1)^{v_p(n)}=1+2\sum_{e=1}^{v_p(n)}(-1)^e\). The resulting divisor coefficients are absolutely summable with weight 1/d, so cancellation persists by truncation. Adding the remaining flips changes at most \(X\sum_{p\in R\text{ outside the truncation}}1/p\) integers up to X. Letting that tail vanish proves the assertion. In contrast, if one can prove cofinite exact coverage of one odd residue modulo each of a set of odd primes with divergent harmonic sum, CRT gives density one. This needs genuinely additive constructions, not ordinal products.

## A bounded H cannot provide a dense menu of rooted responses

For a rooted tournament T let h=H(T), e be its full Hamiltonian paths ending at the root, and u be the ordered covers (P,Q) of T by two paths with P possibly empty and Q starting at the root. These are the coordinates of the [one-arc join law](../session231_overnight_completion_20261109/descent/PROOF.md).

For fixed H₀ the triples with h≤H₀ form a finite union of exact dyadic rays in u, with separate finite base cases for e. To prove this, list the strong components in ordinal order. Every component before the root component goes wholly into P, contributing its ordinary count. The root component contributes its rooted u. Each component after it contributes its two-slot cover count with empties. Concatenation proves the product rule, and a nonempty suffix makes e=0. A singleton contributes 1 before the root and 2 after. Every nonsingleton component has H≥3 and bounded order by the [Busch floor](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3), and there are at most \(\log_3 H_0\) of them. Remove the singleton components other than a singleton root. The remaining finitely many profiles are bases, each generating exactly u multiplied by arbitrary powers of two by appending singleton components. Their endpoints are unchanged in the zero-power case and zero after appending. Thus any construction that wants arbitrary rooted responses must allow h to grow or use a different interface.

A quantitative form is
$$
\operatorname{oddpart}(u)\le h^{\log_5 40},
\qquad h\ge \operatorname{oddpart}(u)^{\log_{40}5}.
$$
For a b-vertex component, inject a two-path cover into a Hamiltonian path plus b assignment bits: greedily choose the head that beats the other available head, preserving each list order, then record which vertices came from each list. Restriction recovers the cover. The unrestricted count Q is even by slot swap and is at most \(2^b H\), so its odd part is at most \(2^{b-1}H\). For a rooted cover the root bit is fixed, giving the same bound directly for u. The SCC product rule and \(H\ge5^{(b-1)/3}\) yield the displayed estimate. Conversely Q≥(b+1)H by cutting a full Hamiltonian path at each of its b+1 possible positions and recovering it by concatenation. These are asymmetric costs: singleton padding can freely grow the dyadic part, not the odd part. The estimate is not claimed sharp.
