# Exact inverse and cyclotomic sieve for a cyclic-wrap family

## Wrap and its arithmetic restrictions

For a nonempty tournament \(R\) of order \(m\), write \(h=H(R)\) and let \(Q\) count ordered two-path covers of \(R\), allowing empty slots. Use the cyclic wrap \(W(R)=C_3[b,a,R\to x]\), with singleton vertices \(a,b,x\). Each Hamiltonian cycle of \(W(R)\) has the unique form \(b\to a\to P\to x\to b\) for a Hamiltonian path \(P\) of \(R\), so its unrooted cycle count is \(h\). Its path count is the inherited \(N=h+2Q\): writing \(B=R\to x\), the three rotations of the quotient triangle each yield \(h\) single-visit paths; the double visit to \(B\) must have the form \(B_1,b,a,B_2\), giving \(2Q-2h\) options because \(P_B(2)=2Q\). There can be no further visit.

In fact \(Q\ge(m+1)h\), with equality exactly for a singleton and the directed triangle. Split any Hamiltonian path at each position, including the ends, to count the covers that concatenate. A nonstrong \(R\) of order at least two has a nonconcatenating cover: split at an SCC boundary and swap the two parts. For strong \(R\) of order at least four, take a directed Hamiltonian cycle by Camion's theorem. If any consecutive triple \(a,b,c\) on the cycle has chord \(a\to c\), use the cover with first part \((b)\) and second part \(a,c,\ldots\), continuing around the cycle to the predecessor of \(a\); the join \(b\to a\) fails. Otherwise all forward distance-two chords point backward. Choose four consecutive vertices \(a,b,c,d\); use parts \((c,a,b)\) and \((d,\ldots)\), continuing the second around the cycle to the predecessor of \(a\). The join \(b\to d\) fails by the backward-chord rule. The singleton and triangle give \((m,h,Q)=(1,1,2),(3,3,12)\), respectively.

It follows that \(N\ge(2m+3)h\). If \(N=kp\), where \(p\) is prime and \(2m+3>k\), then \(\gcd(h,Q)\mid k\): it divides \(N\) but is at most \(h<p\). In the \(7p\) lane this holds for \(m\ge3\); the smaller sources have \(h=1\). Thus a prime target requires coprime \(h,Q\), and a \(7p\) target permits only a common divisor of seven. In particular a triangle SCC rules out both lanes for this wrap regardless of what other components are present.

Both \(h\) and \(Q\) multiply across ordinal sums: for the latter, intersect each labeled path slot with each component and concatenate the pieces in ordinal order. For a source with \(s\) singletons and \(a_i\) copies of nontrivial SCC types \(D_i\) with profiles \((h_i,Q_i)\), the wrap count is
\[
N=\prod_i h_i^{a_i}+2^{s+1}\prod_i Q_i^{a_i}.
\]
Prime \(N\) requires \(\gcd(s+1,a_1,\ldots,a_r)\) to have no odd prime factor: if an odd prime divides all these exponents, \(N\) is a sum of two like odd powers and has a proper factor. The coprimality guard on \(h,Q\) can also see a prime shared between the \(h\)-coordinate of one component and the \(Q\)-coordinate of another.

## A target-first decoder and its exact density loss

Use only singletons and copies of the strong four-vertex \(D\) with out-neighbor bitmasks \((8,1,3,6)\), where bit \(j\) means an arc to vertex \(j\). Its five Hamiltonian paths have vertex words \(0321,1032,2031,2103,3210\). For \(Q\), the two empty-slot cases contribute \(2\cdot5\); the sizes 1 and 3 contribute \(2(3+3+1+1)\) by the four single-vertex deletions; the size 2/2 case contributes 6. Hence \(Q=32\). If \(a\) of these components and \(s\) singletons are used, then
\[
N=5^a+2^k,\qquad k=5a+s+1.
\]
Since \(5^a<2^k\), the prescribed target determines \(k=\lfloor\log_2 N\rfloor\) uniquely. Test if \(N-2^k=5^a\) for an integer \(a\ge0\), then set \(s=k-5a-1\). Accept exactly when \(s\ge0\) and \(s+a>0\). This recovers a physical source profile, not a unique ordering of its SCCs.

For \(a\ge1\), a prime output needs
\[
a\equiv k\pmod2,\qquad \gcd(a,k)\text{ a power of }2 .
\]
The first condition excludes a proper factor 3, since \(5^a+2^k\equiv(-1)^a+(-1)^k\pmod3\); the second excludes all common-odd-power factorizations. They are not sufficient.

Let \(L=\lfloor\log_2 X\rfloor\). Up to the \(k=L\) row, of size \(O(L)\), counting distinct outputs at most \(X\) is counting the lattice triangle \(1\le a\), \(5a+1\le k\le L\); its count is \(L^2/10+O(L)\). For odd \(d\), its points with \(d\mid a,k\) and \(a\equiv k\pmod2\) number \(L^2/(20d^2)+O(L/d+1)\): divide both coordinates by \(d\), so parity agrees in half the points in each integer row up to error at most one. Apply Möbius inversion over odd \(d\le L/5\); the total error is \(O(L\log(2L))\). Since \(\prod_{p\ {\rm odd}}(1-p^{-2})=8/\pi^2\), the survivor count is
\[
\frac{2}{5\pi^2}L^2+O(L\log(2L)).
\]
Thus the two elementary restrictions retain a fraction \(4/\pi^2+o(1)\) of this family's values. This counts candidates, not primes.

## Precise open-problem bridge

Setting \(a=0\) gives the transitive-source wrap spectrum \(1+2^k\), \(k\ge2\). The odd-divisor factorization forces \(k=2^j\) if the output is prime. Conversely every Fermat prime other than 3 occurs as a prime output of this particular all-order family. Thus asking whether transitive-source wraps attain infinitely many primes is *exactly* the [Fermat-prime infinitude problem](https://arxiv.org/abs/2204.08302), not an analogy based on the names of the constructions. This does not reduce Chan–Pak to Fermat primes: other tournament families may also realize those counts or produce all the other primes. The wrap identities and SCC multiplication are inherited; the claims here are the sharp concatenation test and the exact restricted decoder and sieve.

## An infinite prime obstruction for every fixed finite palette

The sparsity of a fixed palette has a stronger, local cause: **there is a reduced arithmetic progression of primes \(\ell\) such that no count in the entire palette family is divisible by \(\ell\)**. This obstruction holds for arbitrary orders, arbitrary multiplicities, and arbitrary changes in the ordinal order of the SCCs.

For an odd prime \(\ell\nmid\prod_i h_iQ_i\), division by \(h\) gives
\[
N/h\equiv 1+2^{s+1}\prod_i(Q_i/h_i)^{a_i}\pmod\ell.
\]
As \(s,a_i\) vary, the product term hits every element of the subgroup of \(\mathbf F_\ell^\times\) generated by \(2\) and the \(Q_i/h_i\). Nonnegative exponents suffice because every generator has finite order; increasing the exponent of \(2\) by its order as needed enforces \(s\ge1\), hence nonemptiness. Thus some wrap count is divisible by \(\ell\) **if and only if** at least one generator has even order. Indeed a subgroup of the cyclic multiplicative group contains \(-1\) exactly when its order is even, and its order is the least common multiple of its generator orders. This is a sign-reachability condition, not the full-generation condition in the [Artin primitive-root problem](https://arxiv.org/abs/math/0412262).

To force the opposite condition uniformly, choose \(M=8\prod_q q\), where \(q\) ranges over distinct odd primes dividing the finitely many integers \(h_iQ_i\). For every prime \(\ell\equiv-1\bmod M\), the number 2 is a quadratic residue since \(\ell\equiv7\bmod8\); for each \(q\), reciprocity and \(\ell\equiv-1\bmod q,\ \ell\equiv3\bmod4\) give
\[
(q/\ell)=(-1)^{(q-1)/2}(\ell/q)=(-1)^{(q-1)/2}(-1/q)=1.
\]
All generators are therefore in the odd-order subgroup of squares. Such primes cannot divide any count of the family; infinitely many exist by Dirichlet's theorem. This excludes both \(\ell\) and \(7\ell\), as well as every other multiple of \(\ell\), as *values of the fixed family*. It is not a global obstruction on tournaments. Even combining finitely many such palettes does not help: take their union before forming \(M\).

For the explicit \((5,32)\) palette the generators are equivalently \(2,5\). Both are squares for primes \(\ell\equiv31,39\bmod40\); here \(\ell\equiv7\bmod8\) and \(\ell\equiv\pm1\bmod5\). Both progressions are invisible to the entire physical family, even with variable \(a\). An all-target proof of Chan–Pak via cyclic wraps must therefore either introduce infinitely many source SCC profiles or leave this form of wrap on some targets.

The obstruction also has an exact prime density on the \(\ell\equiv3\bmod4\) half. Let \(r\) be the rank of the square classes of \(2,Q_i/h_i\) in \(\mathbf Q_{>0}^\times/(\mathbf Q_{>0}^\times)^2\). On that half, a prime is invisible exactly when all \(r\) basis classes have Legendre symbol \(+1\), because the squares are precisely the elements of odd order. The sign \((-1/\ell)\) is an independent character: no product of positive square classes represents the class of \(-1\). On reduced residue classes modulo a common conductor, the \(r+1\) independent quadratic characters consequently take each sign vector equally often (expand its indicator and sum the nonprincipal characters). The prime number theorem in arithmetic progressions then gives relative natural density \(2^{-r}\) among primes \(\ell\equiv3\bmod4\), or a contribution \(2^{-r-1}\) among all primes. No density is asserted here for the other half, where membership among the squares need not force odd order.

Finally, any unit profile \((h,Q)\) that succeeds at a prime divisor \(\ell\equiv7\bmod8\) must have the exact Legendre fingerprint \((hQ/\ell)=-1\): the equality \(h\equiv-2Q\bmod\ell\) gives the general fingerprint \((hQ/\ell)=(-2/\ell)\). In particular a target-dependent choice of source must break the fixed-palette square-class restriction, not merely reorder the same components.

**Finite-rank strengthening.** For exact prime and \(7p\) targets, finitely many SCC types are not needed for the negative result. Even an infinite bank fails if all source ratios \(Q/h\) lie in a fixed finite-rank rational square-class subgroup. Form \(M\) from the prime support of a squarefree basis for that subgroup. For primes \(\ell\equiv-1\bmod M\), any ratio that is a unit has Legendre symbol \(+1\); reducing square factors in its rational expression cannot change this symbol. In a proposed realization of \(\ell\) or \(7\ell\), the bound above forces \(h<\ell\), while \(\ell\mid h+2Q\) then forces \(Q\) to be a unit as well. The necessary negative fingerprint cannot occur. Thus an all-target cyclic-wrap proof requires unbounded square-class rank, not just infinitely many profiles. This strengthening does not claim that infinitely many profiles with fixed square-class rank cannot realize some large *composite* multiples of \(\ell\) with both coordinates divisible by \(\ell\).

## Imported density theorem and a simultaneous-modulus warning

For a finite palette, the exact subgroup criterion identifies the set of
prime divisors of its wrap spectrum, off the finite bad set, with the
primes where the reduction of the rational group
\(G=\langle 2,Q_i/h_i\rangle\) has even order. This is exactly the input of
[Sgobba's unconditional reduction-order density theorem](https://arxiv.org/abs/2110.08911)
with divisor 2. It gives existence and a rational expression for the
**full** prime-divisor density, including the primes \(1\bmod 4\).
The elementary square-class count above gives an explicit obstruction
without doing the higher Kummer calculation. In the transitive case
this reduction recovers the classical \(17/24\) prime-divisor density
(the Hasse case of [Moree's order-divisibility result](https://arxiv.org/abs/math/0407421)).
This is a density of primes *dividing* some output, not primes *equal* to
an output. The latter infinitude question is still the Fermat problem.

For one odd good prime, divisibility by any power \(\ell^e\) has the same
criterion as divisibility by \(\ell\): reduction of the unit groups has
odd-order kernel, so parity of the subgroup order is unchanged. A product
of primes is not the same question. Already for transitive sources with
\(N=1+2^k\), divisibility by a product of distinct odd primes holds for some
\(k\ge2\) exactly when the \(2\)-adic valuations of all their orders
\(\operatorname{ord}_{\ell}(2)\) are equal and positive. Each factor
requires \(k\equiv \operatorname{ord}_{\ell}(2)/2\bmod\operatorname{ord}_{\ell}(2)\).
Unequal valuations are inconsistent; for a common valuation \(u\ge1\)
take \(k=2^{u-1}\) times the least common multiple of the odd parts,
and if necessary an odd multiple to make \(k\ge2\). Thus 3 and 5 are
each prime divisors of infinitely many transitive-wrap outputs, but 15
divides none. This is a literal failure of gluing local witnesses on
an exponent orbit, not a claim of failure for the full tournament spectrum.
