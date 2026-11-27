# Fermat primes as exact readouts for odd-pair matrices and Gram lifts

The earlier [Fermat tournament census](../fermat_monster_hamiltonian_counts_prime_classes_20261102/README.md) records the minimum witness orders \(3,4,6,8,11\) for the five known Fermat primes. The [changed-reader experiment at \(7\cdot65537\)](../chanpak_changed_contacts_20261117/FERMAT_TARGET.md) showed why an exact secondary readout matters, and the [Fermat-neighbor session](../clocks_carriers_genus_20261102/FERMAT_NEIGHBOURS.md) used the dyadic unit groups rather than just a coincidence of five numbers. Here those groups give exact finite-field dynamics, while two of the primes provide no-wrap certificates for the [Moore Gram lift](../odd_pair_matrices_moore_gram_and_pollock_fibers_20261126/README.md). The same primes expose a failure mode for the Pollock Gram matrix.

## 1. Odd pairs split the Fermat field into exactly two long cycles

Let \(p=F_m=2^{2^m}+1\) be prime, \(m\ge2\). The centered all-odd Collatz branch is multiplication of \(x+1\) by \(g=3/2\), and its odd-pair branch by \(g^2=9/4\). Since \(p\equiv5\bmod12\), 3 is a quadratic nonresidue by quadratic reciprocity (the classical [Pépin test](https://mizar.uwb.edu.pl/version/current/html/pepin.html) uses the same fact). Since \(p\equiv1\bmod8\), 2 is a square. The group order is a power of two, so every nonsquare is primitive. Thus \(g\) generates the units and \(g^2\) generates exactly their square subgroup.

The affine odd-pair matrix on \(\mathbb F_p\) fixes \(-1\); away from it there are exactly two orbits, the two square classes of \(x+1\), both of length \((p-1)/2\). A single odd matrix interchanges them. This is the finite-field meaning of the odd pair, not a statement about integer Collatz cycles. For any prescribed finite all-odd length \(L\), the actual integer guard \(x\equiv-1\bmod2^L\) can be imposed independently of the phase mod \(p\) by CRT; no one positive integer follows the all-odd branch forever.

| \(p\) | order of 2 | order of \(3/2\) | order of \(9/4\) |
|---:|---:|---:|---:|
| 3 | 2 | not a unit | not a unit |
| 5 | 4 | 2 | 1 |
| 17 | 8 | 16 | 8 |
| 257 | 16 | 256 | 128 |
| 65537 | 32 | 65536 | 32768 |

This also gives a tournament readout. On \(2^m+2\) vertices take the total order with just the global endpoint arc reversed. Its height is \(F_m\), by the [one-back-arc calculation](../euler_zigzag_minimum_matching_and_hidden_flip_20261126/README.md). Its two-color count is \(2^{2^m+2}+2(3^{2^m}-1)\), so
\[
P_T(2)\equiv 2(3^{2^m}-3)\not\equiv0\pmod {F_m}
\]
for every prime \(F_m\) with \(m\ge1\). Indeed 3 is primitive for \(m\ge1\) by the same nonresidue argument, and \(2^m-1\) is odd. Thus one uniform probe sees a nonzero residue even when the Hamiltonian height vanishes in the Fermat field. This applies to the explicit, not necessarily minimum-order witnesses.

## 2. Moore: exact Gram polynomial before attempting the matching split

Use the preceding [rooted residual notation](../odd_pair_matrices_moore_gram_and_pollock_fibers_20261126/README.md). With \(k\) new branches still missing and all \(q=d-k\) old mutual matchings fixed, assemble the residual matrix \(K\); its diagonal blocks are \(kI_b\), \(b=d-1\), and the off-diagonal blocks must be nonnegative. In any completion \(K=QQ^\top\) with the old-to-new contact matrix \(Q\).

For \(k=1\) the columns of \(Q\) are orthogonal and each has squared norm \(q\). Hence
\[
K^2=qK .
\]
With nonnegative residuals this is equivalent to actual last-star completion: the polynomial forces \(\mathrm{rank}\,K=b\) by its trace, and the anchor Schur complement gives the cocycle from the prior result. At degree 57, \(q=56\), and each entry of \(K^2-qK\) has absolute value at most 56. One exact calculation over the Fermat field \(\mathbb F_{257}\) therefore proves or rejects this integer matrix identity. It is not a probabilistic rank proxy.

For \(k=2\), write \(Q=(Q_U,Q_V)\) for the new stars and set \(S=Q_U^\top Q_V\). This matrix counts routes from \(U\) to \(V\) through old branches. Since only the direct \(U\)-to-\(V\) matching can fill the remaining slots, a completion implies \(S=J_b-P\) for that permutation \(P\). Consequently
\[
Q^\top Q=\begin{pmatrix}qI_b&S\\S^\top&qI_b\end{pmatrix}
\]
has eigenvalues \(2q,0\) on the two branch-constant directions, and \(q-1,q+1\), each with multiplicity \(b-1\), on the zero-sum directions. Indeed the singular values of \(J_b-P\) are \(q=b-1\) and \(1\) on the two respective spaces. The corresponding necessary test on the old data is the cubic identity
\[
K(K-(q-1)I)(K-(q+1)I)=2(q-1)J .
\]
It follows from the Gram spectrum and the all-ones eigenvector, since the number of old leaves is \(qb\). With \(q=55\) the three nonzero eigenvalues are \(110,54,56\), and the right side is \(108J\). Each entry of \(K\) lies in \([0,2]\) and the row sum is 110. Entries of \(K^2,K^3\) are at most 220 and 24200. Thus the absolute difference between the two sides of the displayed identity is at most \(24200+2\cdot3024+108=30356\) per entry, bounding the positive and negative parts separately. **One calculation mod \(65537\) tests it exactly, without integer wraparound**.

Given a coherent two-star split after this spectral test, the candidate final matching is already determined: \(\|S\|_F^2=q^2+b-1=bq=\sum_{ij}S_{ij}\), by the spectrum and the row sums. The entries are nonnegative integers, so \(S\) is 0–1 and \(P=J_b-S\). The remaining old-to-new route verification is still required; a pass of the polynomial alone cannot construct the stars. These are tests for a future complete old branch bank, not claims about a degree-57 graph or arbitrary partial banks. The replay checks the polynomial for all one- and two-branch deletions at one root in each of the Petersen and Hoffman–Singleton graphs (34 total controls) and rejects the width-three false star from the previous note.

## 3. Pollock: a Fermat-prime Gram matrix can lose rank without losing liftability

Under offset scaling \(y\mapsto ay\), the moment Gram transforms by congruence with \(\operatorname{diag}(1,a,a^2)\); in particular \(Q,B,p_4\) scale by \(a^2,a^3,a^4\). Over any Fermat prime \(p\ge5\), squaring loses a sign and the fourth power loses a four-element phase, while cubing on the unit group is a bijection. Thus a nonzero cubic moment determines the phase of any proposed scaling; it does not determine the coordinates. The other moments still test whether the phase is compatible. Squaring the Collatz multiplier and squaring a Pollock phase use the same index-two subgroup and lose the same quadratic-character bit.

More sharply, consider the earlier integral collision \((-17,-5,-4,12,14)\) and \((-16,-10,3,7,16)\), with power sums \((0,670,-630,143554)\) through degree four. At the Fermat prime 5, the common moment Gram is \(\mathrm{diag}(0,0,4)\) modulo 5, only rank one. But each tuple takes every residue mod 5 once. The Jacobian of the **four-moment map** has rank four: choose any four distinct values to obtain a Vandermonde minor multiplied by \(1\cdot2\cdot3\cdot4\), a unit mod 5. Thus a singular Gram is not a singular moment map over a finite field. The loss of rank is from isotropy, not repeated coordinates.

The same rank-four argument works at \(17,257,65537\); these tuples have at least four distinct residues in each. The two unordered reductions are different mod \(p^2\) for every one of these four primes (at 5 their equal mod-5 multisets do not suffice). Holding one coordinate fixed in each branch and lifting the other four by the invertible Jacobian shows: **every four-moment target in the common mod-\(p^2\) cylinder lifts along two distinct unordered \(p\)-adic branches**. One can fix the first moment at zero to stay centered. The replay lifts one perturbation through the next digit at each prime. This is not integer positivity or a proof that the Pollock spectrum is complete.

The practical distinction is that the Moore polynomial has an independent size bound that makes a Fermat check exact over the integers. In Pollock, even a low-rank reduction can conceal two regular lifting branches, and a finite-field lift does not supply a bounded positive integral witness.

Run `python3 replay.py` for the controls.
