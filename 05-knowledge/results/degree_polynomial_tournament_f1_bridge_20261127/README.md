# Degree polynomials, characteristic one, and what a tournament spectrum forgets

2026-11-27. Read Long's [single-scale curve RH](https://github.com/long-mathematics/abelian-degree-polynomials/blob/main/papers/curve_rh_single_scale.tex),
[interpolation](https://github.com/long-mathematics/abelian-degree-polynomials/blob/main/papers/abelian_comparison_elementary.tex),
and [Rosati](https://github.com/long-mathematics/abelian-degree-polynomials/blob/main/papers/abelian_comparison_streamlined.tex)
notes alongside our [RH and absolute geometry](../chanpak_rh_absolute_geometry_20261124/RH-F1-READING.md)
and [count monoid](../residues_monoid_tower_tree_20261116/MONOID.md) work.
These yield an actual abelian and toric realization of tournament matrices,
a sharp spectral-variance certificate for the Paley quotient, and an
explicit limit on using degree polynomials to prove Chan–Pak.

## 1. What transfers from the RH argument

The curve argument controls a nonnegative characteristic polynomial on
the real axis using the nef classes belonging to the sum and difference
of maps. For a polarized endomorphism of weight d, use the same scale
2d^(n/2) on the nth iterate. The nth root of the absolute value of each
linear factor tends to max(sqrt(d), |alpha|). This bounds each root
modulus from above, and conservation of the product of the moduli
forces equality. The argument does not have to rule out simultaneous
cancellation in individual traces. The second paper's bounded-block
Vandermonde inversion is the alternative when a trace bound is what
one actually has. The third paper's positive Rosati adjoint explains
the stronger norm and composition statements, not just eigenvalue moduli.

This is an expository proof of classical curve/function-field RH.
The [Alpöge–Claude number-field result](https://alpo.ge/remark.pdf),
with [Lean/Mathlib statement and proof verification](https://leanprover-community.github.io/blog/posts/riemann-hypothesis/),
is separate. In our previous reading the number-field proof detects a
hypothetical outlier through a one-sided analytic window estimate;
the curve proof detects an outlier through an asymptotic product of
linear factors. Both require an independent positive upper bound.
Neither derives one from a formal Euler product or from the existence
of a spectral interpretation.

## 2. The same integer matrix in two geometries

For an integral square matrix U and an elliptic curve E, use U to
act on E^N by integer combinations of the factors. Product polarization
makes transpose the Rosati adjoint. On the Tate module the action is
U tensor the identity on the two-dimensional Tate module of E, so

    deg([a] - [b]U) = det(a I - b U)^2.

Here degree is the full geometric degree, zero for non-isogenies.
In particular U^T U=d I gives a genuinely polarized abelian
endomorphism, not just a resemblance of spectra. In the Rosati note's
normalization, the mixed-degree polynomial for this product is
det(s I + U^T U). It equals (s+d)^N in the polarized case. One way to
see the product identity is to compare the Tate characteristic
polynomial in the square-root formula from the Rosati paper.

There is also a literal multiplicative model over the pointed
monoid {0} ∪ Z^N (the split torus in [Deitmar's monoid-scheme sense](https://arxiv.org/abs/math/0404185)).
Here 0 is absorbing and the copy of Z^N is the Laurent-monomial
group. U acts on the exponent lattice; integer exponents and
composition require no sum of monomials. After base change, for
q a prime power, consider the monomial map whose lattice map is qU.
Whenever q^n U^n-I has nonzero determinant, the nth fixed-point
count over the algebraic closure is

    #Fix((qU)^n) = |det(q^n U^n - I)|.

Indeed the fixed-point characters are maps from the lattice cokernel
to the multiplicative group. The cokernel has order prime to the
characteristic, since its defining matrix reduces to -I there.
The integer torus map U itself need not be invertible and qU need
not be the Frobenius of an arithmetic variety; this counts fixed
points of a chosen map. Torus degrees use one copy of the determinant,
whereas the elliptic product uses two. The monoid construction
does not give the positive polarization for free, and does not
produce the zeta function of Spec Z.

This is a different use of characteristic one from
[Connes–Consani's scaling-site covers](https://arxiv.org/abs/2501.06560):
there a prime orbit has length log p and its monodromy encodes
Frobenius. Here we have a finite rank lattice of exponent coordinates.
In the earlier count-monoid note the characteristic-one linear-order
picture was a third, analogy-only use. None of the three meanings
should be substituted for another.

## 3. A sharp quartic deficit for the Paley-59 quotient

Let S be the skew sign matrix of an odd-order tournament of order m,
and set C = [0, 1^T; -1, S], of even order N=m+1.
The doubly regular condition is exactly C^T C=mI:
the off-border blocks give S1=0 and the remaining block gives
SS^T=mI-J. So a doubly regular tournament gives a true polarized
map on E^(m+1), with characteristic polynomial on the exponent
lattice (t^2+m)^((m+1)/2). In unfactored form the Tate degree
polynomial is its square. For example when m>1 the ordinary finite product
det(I-C m^(-s)) has a critical line Re(s)=1/2 in the elementary
root-modulus sense. This is not a prime Euler product. Taking the
reciprocal gives poles, not zeros; compare the correction in our
[historical tournament L script](../../../04-computation/tournament_L_function.py).

For any skew sign C of size N, write its eigenvalues as ±i sqrt(t_j),
j=1,...,N/2, with zeros paired. The t_j are the squared singular values
counted once per pair. Put G=CC^T and D=||G-(N-1)I||_F^2.
Since the diagonal of G is N-1, the coefficient of x^(N-4) in
det(xI-C) is exactly

    e_2(t) = choose(N/2,2)(N-1)^2 - D/4.

Indeed sum t_j=N(N-1)/2 and sum(t_j-(N-1))^2=D/2.
The same coefficient in the monic square root gives an observable
of the degree polynomial. Equivalently e_2 is the sum of principal
4-by-4 Pfaffian squares. The relationship to diamonds and related
extremal Gram arguments is studied by
[Belkouche–Boussaïri–Lakhlifi–Zaidi](https://arxiv.org/abs/1906.04672).

**Parity certificate.** For N congruent to 2 mod 4, the sharp
general lower bound when a two-row deletion from a skew conference
of order N+2 is available is

    D >= 2N(N-2).

Write p_i as the product of all off-diagonal entries in row i.
For distinct i,j, the product of row-wise signs on their N-2
common coordinates is -p_i p_j. The dot product G_ij is 2 mod 4
if the row products agree, and 0 mod 4 otherwise. The number
of ordered same-sign pairs is minimized by two equally sized
groups. Each such pair spends at least four on D.
Equality requires equal group sizes, absolute dot products two
within groups, and zero across groups. This does not classify all
equality cases or assert any determinant maximization.

**Deletion realizes equality.** Delete two rows and columns from a
skew conference of order N+2, keeping C, and let B be the two
deleted rows restricted to the retained columns. Orthogonality of
the full matrix, and the deleted skew 2-by-2 block, give

    C^T C=(N+1)I-B^T B,     BB^T=N I_2.

Thus the squared singular values are 1 twice and N+1 on the remaining
N-2 dimensions. So D=2N(N-2), as claimed. In particular the
57-vertex tournament obtained from Paley-59 has a 58-dimensional
border with

    det(tI-C) = (t^2+1)(t^2+59)^28,    D=6496.

It minimizes the quartic deficit (1,624) against the flat benchmark
among all bordered 57-vertex tournaments. If p=N+1 is prime, the
border's Smith form has (p+1)/2 ones and (p-3)/2 copies of p.
Indeed modulo p, B^T B has a two-dimensional nondegenerate image
as BB^T=-I there. C preserves its image and kernel because
C^T C=-C^2. On the image C^2=-I, and on the p-3 dimensional
kernel C^2=0, so the rank is at most (p+1)/2. Its determinant
is p^((p-3)/2), giving the opposite rank inequality and each
nontrivial Smith factor equal to p.

The quotient Gram entry for distinct retained slots i,j is
1-(S^2)_ij = -u_i u_j-v_i v_j for the two deleted sign columns
u,v, exactly the correction behind our
[off-branch defect classes](../moore_57_slotwise_balance_20261127/README.md).
This proves quotient optimality for a specific averaged spectral
criterion, not optimality for the dense encoding's maximum entry.
For a conditional degree-57 Moore adjacency A, the symmetry-centered
operator 2A+I on the sum-zero lattice instead satisfies U^T U=225I
there, with eigenvalues +15 and -15 of multiplicities 1729 and 1520.
This is a scalar Gram identity on an invariant sublattice, not a
construction of A. Dividing by 15 does not preserve that lattice:
apply 2A+I to a difference of adjacent vertex basis vectors and read
the coordinate -1. The associated normalized involution should
not be mistaken for a permutation of the natural vertex basis.

## 4. A no-go theorem for unmarked degree data in Chan–Pak

On an ordered tournament, flip just the arc from i to j with j-i>=2.
It has 1+2^(j-i-1) Hamiltonian paths. Every path not using the
reversed arc is the increasing order (for j-i >=2). A path using
j->i has one increasing prefix ending at j and one increasing
suffix starting at i. Each intermediate label can independently
go on either side; labels outside [i,j] have a forced side.

On the even cyclically ordered set of N labels 0,...,N-1,
start from the skew matrix of the transitive order with only
the arc 1->3 reversed. With border root 0 the surviving tournament
has gap two and hence 3 Hamiltonian paths. Re-root at 2 and switch
each old coordinate with label <2 by minus one; list the remaining
labels in order 3,4,...,N-1,0,1. The new root row is again all +1;
all ordinary transitive arcs still agree with the new order.
The same single reversal is now at gap N-2. The second tournament
has 1+2^(N-3) Hamiltonian paths. For every even N>=6 the two
bordered matrices are conjugate by a signed permutation.
They therefore agree on all unmarked determinant degrees,
on all spectral moments, and on the above torus fixed-point
counts wherever finite. The discrepancy is exponential even
though the conjugating maps exist in the monoid geometry.
The border is a choice of point; transport that forgets it
does not preserve Hamiltonian counts.

Even supplementing the bordered characteristic polynomial with
the polynomial of the root-deletion does not recover the count.
On labels 0,...,6 reverse respectively {02,04,06} and {02,06,13}
from the transitive order. Both are strong (the arc 6->0 closes
the common increasing Hamiltonian path). The replay computes
Hamiltonian counts 31 and 41. Both unbordered polynomials are

    t (t^6 + 21t^4 + 99t^2 + 71),

and both bordered polynomials are

    t^8 + 28t^6 + 182t^4 + 252t^2 + 49.

These together determine all scalar resolvent moments at the
border by the cofactor formula, so even that pointed spectral
observable is insufficient. This is a counterpart of our
[append-interface example](../tournament_computational_interfaces_fractran_20261127/README.md),
where equality of H hid the future response. Here equality of the
tested spectral data hides H itself. The natural next refinement
is to retain multiple marked contacts and a no-repeat condition,
rather than pushing the unmarked zeta comparison harder.
None of these arguments resolves the spectrum conjecture or
strengthens the unrestricted entrywise Moore-encoding bound.

Run node replay.js for exact integer checks of the deletion
spectra and ranks, all labelled matrices of size six for parity,
and direct Hamiltonian counts for the collision families.
