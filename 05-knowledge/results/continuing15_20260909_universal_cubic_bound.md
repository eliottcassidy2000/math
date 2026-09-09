# The universal weighted cyclic bound from a real spectral inequality

**Status: PROVED ANALYTICALLY + FINITE-EXACT; INDEPENDENTLY AUDITED.** This closes the
arbitrary-support extension left unproved in the incoming THM-4466.
It is an algebraic matrix theorem in every finite dimension. It does not
assert an Euler trajectory remains in the oriented cone or settle a PDE.
No external priority claim is made.

The [independent referee](continuing15_20260909_universal_cubic_audit.md) accepts this result;
this filed copy changes only audit status and adds this link.

## 1. Inheritance and the changed representation

The incoming [THM-4466 / sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes](../../01-canon/theorems/THM-4466-sharp-tournament-cubic-bound-for-common-edge-and-substitution-classes.md)
proves the sharp cyclic bound for common-edge supports and a substitution
grammar. It explicitly leaves arbitrary overlapping cyclic supports
unproved. The [THM-4465 / weighted-tournament-shear-production-and-contact-kernel](../../01-canon/theorems/THM-4465-weighted-tournament-shear-production-and-contact-kernel.md)
identity is F=3C-T=4 tr(sym(M) skew(M)^2), with C,T the weighted cyclic
and transitive triple sums and E=||M||_F^2.

The closest proved mechanism is the common-edge Cauchy estimate; its
canonical equality hostile is one edge of weight two with four return
paths of unit weights. Equality is not restricted to one triangle.
The underused operation is passage to the full complex spectrum of a
REAL matrix. It preserves both trace moments and gives a one-sided
Frobenius energy bound. Conjugate pairing, lost by treating the spectrum
as an arbitrary complex list, is the decisive extra coordinate.

## 2. A stronger matrix theorem

**Theorem.** For any real n by n matrix M with tr(M^2)=0,

    |tr(M^3)| <= ||M||_F^3 / sqrt(3).                    (1)

If M!=0, equality on the positive side holds exactly when M is normal
and its eigenvalue multiset is

    r, r*omega, r*omega^2, 0,...,0,
    r>0, omega=exp(2*pi*i/3).                           (2)

In particular such a matrix has rank three. Negative equality is obtained
by replacing M with -M. The zero matrix is the separate equality case.
Neither nonnegativity, zero diagonal, nor tr M=0 is assumed in (1).

Let the eigenvalues, with algebraic multiplicities, be lambda_i=a_i+i b_i.
They occur in conjugate pairs because M is real. Put

    S=sum a_i^2=sum b_i^2, E=||M||_F^2.

The equality defining S follows from tr(M^2)=0. Complex Schur
triangularization and unitary invariance of the Frobenius norm give

    E >= sum |lambda_i|^2 = 2S,                         (3)

with equality exactly when M is normal. This can be seen directly by
writing E as the squared diagonal norm plus the sum of squared strictly
upper-triangular Schur entries.

If S=0, all eigenvalues vanish and tr(M^3)=0. Otherwise

    tr(M^3)=sum a_i^3-3 sum a_i b_i^2.                  (4)

If every nonreal eigenvalue has nonnegative real part, the second sum
is nonnegative, and (4) is at most sum a_i^3<=S^(3/2), which is strictly
smaller than (2S)^(3/2)/sqrt(3).

In the remaining case let -z<0 be the smallest real part of a nonreal
eigenvalue. At least TWO eigenvalues have this real part, since its
conjugate is also present. Their cubes contribute -2z^3. The cube sum
of the remaining real parts is at most (S-2z^2)^(3/2): bound each signed
cube by its absolute cube and use the l3-versus-l2 norm inequality.
Also a_i>=-z whenever b_i!=0. Thus

    tr(M^3) <= (S-2z^2)^(3/2)-2z^3+3zS.               (5)

Set w=z/sqrt(S), so 0<w<=1/sqrt(2). The right side divided by S^(3/2)
is f(w)=(1-2w^2)^(3/2)-2w^3+3w, and

    f'(w)=3 sqrt(1-2w^2) [sqrt(1-2w^2)-2w].

It increases up to w=1/sqrt(6) and decreases afterwards. Its maximum is
2sqrt(2/3). Combining with (3) proves the positive side of (1).
Apply the same result to -M for the absolute-value statement.

For nonzero positive equality, all these inequalities must be equalities.
The remaining cube norm has precisely one nonzero real part, positive.
At w=1/sqrt(6) it equals 2z. Equality in the b_i^2 bound puts every
nonzero imaginary part at real part -z; hence only the already selected
conjugate pair can be nonreal. Its imaginary parts are plus/minus sqrt(3)z.
All other eigenvalues are zero. Equality in (3) is normality. This gives
(2) with r=2z, and direct substitution proves the converse.

## 3. The entire weighted oriented cone is now closed

Let M have nonnegative real entries, zero diagonal, and M_ij M_ji=0.
Missing edges are allowed. The intrinsic pair relation is its positive
support orientation; no ties are artificially oriented for a proof.
Then tr(M^2)=0 and tr(M^3)=3C: a positive directed triangle has its
three possible starting vertices in the trace. Therefore

    F=3C-T <= 3C <= E^(3/2)/sqrt(3).                    (6)

This holds for EVERY dimension, EVERY nonnegative choice of amplitudes,
and EVERY cyclic overlap pattern. It removes the common-edge and
substitution restrictions of THM-4466. The actual shear identity is
inherited from THM-4465; it is not reinterpreted as tr(M^3)=F.

The spectral characterization (2) is an exact equality criterion for
3C. Equality for F additionally requires T=0. The zero matrix again
has equality. We do not assert here that every equality support has a
particular combinatorial partition; the matrix criterion is complete.

The incoming six-vertex example has E=12,C=8,T=0 and F=24, with
3F^2=E^3. Its matrix is normal of rank three, so it satisfies the new
criterion despite its multiple return paths. Disjoint equal-weight
triangles are strict unless only one has nonzero energy. A diagonal
similarity can preserve trace moments while increasing Frobenius energy;
normality is load-bearing, and similarity is not a lawful energy gauge.

## 4. Scope, tests and next boundary

The companion program exhausts all oriented matrices of orders one
through four with each unordered pair absent or assigned weight one or
two in either direction. It separately tests larger fixed-seed weighted
supports, direct matrix/triple identities, equality families, nonnormal
similarity, and rational orthogonal images lying outside the nonnegative
cone. All comparisons are exact; the finite sample does not prove (1).
The proof in Section 2 supplies the unrestricted statement.

Run `python 04-computation/continuing15_20260909_universal_cubic_bound.py`
and the corresponding optimized command after repository relocation.

The source-to-target map is an oriented amplitude matrix to its spectrum.
Trace moments are preserved; Frobenius energy is bounded from below,
with nonnormal departure retained as a strictness witness. The spectrum
forgets the actual support and boundary contacts, so it does not replace
THM-4465's quadratic contact kernel for nonuniform gluing. Nor does (6)
supply physical pressure, trajectory invariance, LRC ownership, or the
Hamiltonian/Pfaffian inequality H>=disc. Those problems remain **OPEN**.
