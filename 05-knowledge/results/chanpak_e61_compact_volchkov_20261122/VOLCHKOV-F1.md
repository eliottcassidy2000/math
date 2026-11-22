# The numbers in Volchkov's kernel, and a genuine toric carrier for the local seam law

**CITED + PROVED extension.** The displayed identity is exactly [Volchkov's 1995 criterion](https://umj.imath.kiev.ua/index.php/umj/article/view/5435). [Alpöge's November 20 announcement](https://alpo.ge/remark.pdf) and the [November 21 Lean/Mathlib verification](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/) give the classical Riemann hypothesis needed to use that criterion as an equality. The independent verification compared against Mathlib's own statement and checked an inert proof export with multiple kernels. The announcement describes a Stepanov-inspired use of multiplicity through iterated theta lifts and an integral-structure divisibility argument. That is inspiration to build a certificate that preserves signs and multiplicity; it is not a black box producing Hamilton paths or an absolute integer inverse.

## Read the constants as a kernel, not as a list of coincidences

Denote the integral with \(\log|\zeta|\) replaced by \(\log|F|\) by \(\mathcal V[F]\). Its even kernel is literally
\[
 K(t)=\operatorname{Re}(1+2it)^{-3}
     =\frac{d}{dt}\frac{t}{(1+4t^2)^2}.
\]
The \(4\) is the square of \(2\), the \(12\) is \(3\cdot4\) from the real part of the third power, and the initial \(1\) is the unit in that resolvent. The vertical boundary \(1/2\) is the critical line; the lower \(t\)-endpoint \(0\) picks the positive half of the even integrand. Since \(t/(1+4t^2)^2\) vanishes at both endpoints, the zero-frequency part vanishes. For \(u\ge0\),
\[
 \int_0^\infty K(t)\cos(ut)\,dt
     =\frac{\pi}{32}u^2e^{-u/2}.
\]
For a direct proof, express \((1+2it)^{-3}\) as the Fourier transform with negative sign of \(v^2 e^{-v/2}{\bf1}_{v\ge0}/16\). For positive frequency only one of the two symmetrized terms contributes. This gives half its full-line Fourier transform, and the half-line cosine integral is half again. Thus \(32=2^5\) in the cosine measure, not an empirical path count. The independent \(3-\gamma\) on the right is the constant in Volchkov's theorem; \(\gamma\) is the Euler–Mascheroni constant. One should not identify the \(3\) there with an exponent on a tournament module without a map.

## A Volchkov identity for the divisor-mass lanes

For an integer \(k\ge1\) and a finite set of primes \(P\), define
\[
 f_{k,P}(n)={\bf1}_{p\nmid n\ (p\in P)}\,\frac{\sigma_k(n)}{n^k},
 \qquad
 M_{k,P}(s)=\sum_{n\ge1}\frac{f_{k,P}(n)}{n^s}
          =\zeta(s)\zeta(s+k)\prod_{p\in P}(1-p^{-s})(1-p^{-s-k}).
\]
The Dirichlet identity follows by convolution for \(\Re s>1\), and then gives continuation. For \(k=1,P=\{2\}\), \(f\) is the normalized odd local mass of the repository's [F1 seam family (THM-5134)](../s92_prime_powers_quaternions_bases_20261111/README.md): the normalization removes \(8/\pi^2\). For \(k=3,P=\varnothing\), \(f\) is the normalized E8 shell read \(r_{E_8}(n)/(240n^3)\) from THM-5135 there. No tournament family is identified with the E8 read.

The following shifted Volchkov statement is equivalent to RH for every \(k,P\), and hence follows from the cited proof:
\[
 \mathcal V[M_{k,P}]
 =\frac{\pi}{32}\left[
    3-\gamma-\frac{\zeta'(k+1)}{\zeta(k+1)}
      -\sum_{p\in P}\log p\left(\frac1{p-1}
                               +\frac1{p^{k+1}-1}\right)\right].
\]
Indeed, write \(M=\zeta E\). The absolutely convergent log expansion for \(E\) on \(\sigma\ge1/2\) uses coefficient \(p^{-mk}/m\) outside \(P\) and \(-1/m\) inside \(P\) at phase \(u=m\log p\), with common factor \(e^{-u\sigma}\). Integrating in \(\sigma\) contributes \(e^{-u/2}/u\), and the kernel contributes \(u^2e^{-u/2}\pi/32\). Each resulting term is \((\pi/32)\log p\) times \(p^{-m(k+1)}\) outside \(P\) or \(-p^{-m}\) inside \(P\). Summing gives exactly the Euler correction. All exchanges are justified by absolute convergence, using integrability of \(K\); the finite removed primes are harmless. Subtracting it leaves the original Volchkov identity, giving the reverse implication. The odd F1 case has correction \(-\zeta'(2)/\zeta(2)-(4/3)\log2\). Also \(M\) has exactly the nontrivial zeros of \(\zeta\) in \(\Re s>0\), away from its pole at 1, since its other factors there are nonzero.

## F1 the join family versus \(\mathbb F_1\) the geometry

These names are distinct. The [earlier Borger/toric carrier](../ports_owners_f1_20261025/f1_bridge.md) implemented commuting Frobenius power maps on a conic, but proved that its consecutive Euclid-parameter tournament subfamily was not stable under powering. There is a different, very small honest carrier for the local divisor-mass law. Take
\[
 X=\operatorname{Spec}\mathbb Z[1/P]\ \sqcup\
              \mathbb A^k_{\mathbb Z[1/P]}.
\]
At every remaining prime and extension degree \(m\), \(X\) has \(1+p^{km}\) points. Its local zeta factor is \(1/((1-T)(1-p^kT))\), and its arithmetic zeta evaluated at \(s+k\) is precisely \(M_{k,P}(s)\). The polynomial point count splits into a point and a \(k\)-dimensional affine piece. The maps on affine coordinates \(x_i\mapsto x_i^r\), and identity on the point, commute under composition \(r,s\mapsto rs\) and reduce at each remaining prime to Frobenius. Thus the elementary carrier really has the multiplicative lift data in the sense of [Borger's \(\Lambda\)-rings paper](https://arxiv.org/html/0906.3146v1); it doesn't pretend there is an ordinary field of cardinality one. It geometrizes the shifted **local surrogate**, not the literal path-set fiber. This distinction is compulsory: the surrogate is positive at the absent counts 7 and 21.

## Signed cancellation for the inverse local mass, not for tournament counts

Let \(\mu_P(n)=\mu(n){\bf1}_{p\nmid n\ (p\in P)}\). The Dirichlet inverse of \(f_{k,P}\) is \(g=\mu_P*(\mu_P/\mathrm{id}^k)\). At retained primes its local coefficients are \(1,-1-p^{-k},p^{-k},0,\ldots\). The bound
\[
       \sum_{n\le x}g(n)=O_{k,P,\epsilon}(x^{1/2+\epsilon})
\]
for every \(\epsilon>0\) is equivalent to the classical [Mertens criterion, Titchmarsh Thm. 14.25(C)](https://sites.math.rutgers.edu/~zeilberg/EM18/TitchmarshZeta.pdf). To pass from \(\mu\) to \(\mu_P\), convolve with the product of geometric prime-power series at the finitely many primes in \(P\), whose weighted sums converge. Convolve next by \(\mu_P/\mathrm{id}^k\), whose absolute sum with weight \(n^{-1/2-\epsilon}\) converges. Conversely convolve \(g\) with \({\bf1}_{p\nmid n\ (p\in P)}/\mathrm{id}^k\) to recover \(\mu_P\), then with the finite factors \(1-p^{-s}\) to recover \(\mu\). Each operation preserves the exponent of the bound by absolute convolution estimates. The certified RH therefore yields signed cancellation in this inverse mass; it does not give cancellation of Liouville signs on actual tournament words or identify any desired arc flip.

A small cross-index isolates actual repo uses of these integers. Zero is the neutral-flip read, and one is the [transitive unit of the ordinal monoid](../residues_monoid_tower_tree_20261116/MONOID.md), not the first forbidden count. Two enters the literal seam and the two-slot product; three is the number of vertices in the cyclic wrapper and the prime that kills the [triangle receiver](ORDER-123.md). The separate [Leech-phase construction](../chanpak_leech_perfect_phases_20261122/EXPONENT-61.md) uses an actual 12-step exponent phase; the same integer in \(K(t)\) comes instead from \(3\cdot4\). Four is both the square in the denominator normalization and the size of the hard tagged receiver. Thirty-two is a kernel normalization, but is also a genuine binary place value in the receiver's forced five final ones, \(y=32z+31\), and the block size in the [first-32 odd address calculation](../residues_monoid_tower_tree_20261116/MONOID.md#2-the-first-32-odd-numbers-3--1--28-and-the-address-coordinate). These roles are all explicit; none licenses replacing one of their operators by another.

The concurrent [tournament slot-kernel calculation](../chanpak_decimal_kernel_20261122/README.md) uses the same analytic kernel on a **finite** slot counting zeta product and extracts directed-cover coefficients, rather than applying it to the **infinite** divisor-density series above. Its primary source for the term \(3-\gamma\) is [Sekatskii–Beltraminelli–Merlini, Theorem 5a and Remark 2](https://umj.imath.kiev.ua/index.php/umj/article/download/2568/1895/6147): the Laurent coefficient of \(\zeta'/\zeta\) at 1 supplies \(\gamma\) and its finite correction \(3/(2(1-b))\) at \(b=1/2\) supplies 3. The two kernel calculations are compatible but use different measures.
