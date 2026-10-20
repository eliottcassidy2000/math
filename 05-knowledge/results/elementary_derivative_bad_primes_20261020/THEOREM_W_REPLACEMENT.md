# Drop-in elementary replacement for Theorem W(ii)

2026-10-20. **PROVED; independently audited.** This supplement gives an elementary replacement for the infinitude argument in [THM-4724, Theorem W(ii)](../../../01-canon/theorems/THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md), without modifying that incoming source. See the [independent audit](audit/REVIEW.md) and [proof details, controls, and scope](README.md).

**Theorem.** For every nonzero integer \(c\ne1\), infinitely many odd primes \(p\nmid c\) satisfy \(c\in-K_{\max}(p)\). More strongly, if \(|c|\ge2\), then for every odd prime \(\ell\nmid c+1\), every prime divisor of
\[
N_\ell(c)=(c^\ell+1)/(c+1)
\]
has \(\operatorname{ord}_p(c)=2\ell\).

**Proof.** For \(|c|\ge2\), the alternating quotient \(N_\ell(c)\) is a positive odd integer greater than one. Its polynomial expression gives
\[
\gcd(N_\ell(c),c)=1,\qquad
\gcd(N_\ell(c),c+1)=\gcd(\ell,c+1)=1.
\]
If \(p\mid N_\ell(c)\), then \(p\) is odd, \(c^\ell=-1\pmod p\), and \(c\ne-1\pmod p\). The order therefore divides \(2\ell\), does not divide \(\ell\), and is not 2. Since \(\ell\) is prime, it is exactly \(2\ell\). In particular \(2\ell\mid p-1\), so \(p\ne\ell\). Distinct allowable \(\ell\) yield disjoint prime-divisor sets, and infinitely many odd primes avoid the fixed integer \(c+1\).

In a cyclic group of even order, the negative of the maximal odd-order subgroup consists exactly of the elements whose order has 2-adic valuation one. This gives the stated coset condition. Finally \(c=-1\) works at every odd prime. The exception \(c=1\) has order 1 at every prime; \(c=0\) is outside the nonzero-dilation domain. ∎

**Derivative consequence.** Choose one \(K_{\max}(p)\)-coset from each antipodal pair to obtain an invariant tournament connection set \(S\). At every prime in the theorem, every such \(S\) satisfies
\[
1_S(cx)=1-1_S(x),\qquad D_c(S)=\mathbb F_p^\times.
\]
Thus the output is the sign-odd complete graph, whose vertex transpositions have edge-reversal sign \(-1\), even though its vertex degrees are even. Every fixed nontrivial nonzero derivative therefore has infinitely many failures of sign-even-valuedness.

**Prime-power extension.** The reduction kernel of \((\mathbb Z/p^r\mathbb Z)^\times\to\mathbb F_p^\times\) has odd order \(p^{r-1}\). Thus the order's 2-adic valuation is unchanged, and the bad coset condition lifts to every \(p^r\). For ring circulants, the maximal odd unit subgroup has no orbit containing both \(x\) and \(-x\), even for nonunit \(x\ne0\). Choosing one orbit from each antipodal pair therefore gives the same sign-odd complete output on all nonzero residues. This is a ring statement, not a finite-field extension.

**Scope.** No analytic prime-distribution theorem is used. The proof does not prescribe an arithmetic progression beyond \(p\equiv1\pmod{2\ell}\) for the constructed witnesses and does not give a density. Separately, the established [Boolean-affine obstruction](../affine_circulant_ideal_obstruction_20261020/AFFINE_CIRCULANT_IDEAL_OBSTRUCTION.md) already rules out every derivative type bijection at every admissible prime \(p\equiv1\pmod4,\ p>5\). Neither statement is a no-go theorem for arbitrary nonlinear local dilation rules.
