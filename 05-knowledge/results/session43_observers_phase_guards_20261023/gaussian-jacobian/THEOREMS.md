# Finite frames, Gaussian charge, and the exact anti-aliasing compiler

Status: elementary proofs below; finite checks are regression certificates, not proofs of the general Gaussian moments or planar Jacobian conjectures. No priority claim.

## 1. The finite-frame power-null obstruction

Let a finite group G of size N act by automorphisms on a reduced commutative Q-algebra A. Write R(f)=N^{-1} sum_g g(f). If R(f^j)=0 for j=1,...,N, then f=0.

Proof. Put a_g=g(f). The power sums p_j=sum_g a_g^j vanish for j<=N. Newton identities, valid because 1,...,N are invertible, recursively give every elementary symmetric polynomial e_j(a_g)=0. Consequently product_g (T-a_g)=T^N. Substitute T=f, one of the roots: f^N=0. Reducedness gives f=0. The bound N is sharp uniformly: C_N acting on C[z] by z->exp(2 pi i/N)z has R(z^j)=0 for 1<=j<N and R(z^N)=z^N.

Pointwise refinement: for functions on a finite orbit with d distinct nonzero values, vanishing of its first d positively weighted moments forces all values zero. Indeed group equal values lambda_i with positive masses w_i. The matrix (lambda_i^j)_(1<=j<=d) is an invertible Vandermonde times a nonzero diagonal. Thus it cannot annihilate the nonzero mass vector.

This applies to D4 polynomial Reynolds averaging, finite observer-orbit averages, and positively weighted functions on a finite Walsh cube. It does NOT apply to signed averaging, nonreduced algebras, positive characteristic dividing N, or infinitely supported Gaussian expectation. Controls: in Q[epsilon]/(epsilon^2), trivial-group averaging has epsilon^j=0 only for j>=2, so the hypothesis j=1 still fails; a cleaner reducedness control is a nontrivial C2 action epsilon->-epsilon, where all averaged positive powers vanish although epsilon is nonzero.

For circular Z=(X+iY)/sqrt(2), E[Z^m]=0 for every m>=1 while Z is nonzero. A finite frame cannot reproduce this property at every degree. This is a structural obstruction to transporting finite signed-Walsh cancellation as if it were Gaussian power-nullness, not an obstruction to useful bounded-degree exact compilation.

## 2. A repaired finite-frame compiler, with explicit bandwidth

On C[Z,W], let integer charge of Z^a W^b be a-b, and P_L average the C_L action (Z,W)->(zeta Z,zeta^{-1}W). Let P_infinity be the circle Reynolds operator. Then P_L keeps exactly charges divisible by L; P_infinity keeps exactly charge zero.

For every polynomial F whose nonzero charges have absolute value <L, P_L F=P_infinity F. More generally equality holds exactly when every nonzero charge component at a multiple of L vanishes. Distinct charge components cannot cancel as polynomials.

If P has charges between -D and D and Q between -E and E, L>mD+E certifies P_L(Q P^m)=P_infinity(Q P^m). This is a sufficient bound, not a minimal one. Exact support propagation gives a stronger adaptive test: (m*S(P)+S(Q)) intersect (L Z) must be contained in {0}; coefficient cancellation can improve the test further. Shared phases must be retained until multiplication: P_L(FG) generally differs from P_L(F)P_L(G), e.g. F=Z,G=W,L>1.

Combining cyclic guards is exact: P_a P_b=P_lcm(a,b). This follows monomialwise from simultaneous divisibility. In particular P_8 P_9=P_72, not P_8 or P_9. With P=Z, all powers before 72 are killed but the 72nd survives as Z^72. The number 8/9 here is NOT the spectral gap 8/9 from the separate tournament coupling theorem.

After P_infinity, Gaussian expectation is the radial linear functional L((ZW)^k)=k!. For bounded-degree exact P_L=P_infinity the same radial contraction is valid. If instead a finite-angle radial sampling model is used, aliases are real: Z^72 becomes r^72 under 72-angle averaging, and its radial mean is 36!, although the true Gaussian moment is zero. There is therefore no fixed 8,9,72 or other finite angular guard that certifies all powers of arbitrary polynomials.

## 3. Common-frame relative phase is a necessary multiplicative datum

For a simultaneous angular action on two vectors Z_1,Z_2, mixed monomials are selected by TOTAL charge. Independent averaging selects zero charge in EACH vector. Z_1 W_2 survives simultaneous averaging but vanishes under independent averaging. Setting Z_2=Z_1 produces |Z_1|^2 with Gaussian expectation 1. Thus independent orbit compression cannot preserve shared incidence contractions. The exact repair is the total-charge fiber (or equivalent relative phase), followed only then by the radial functional. Session42's mixed D4 separating invariants address finite-frame reconstruction; the infinite charge layer addresses a different loss and neither replaces the other.

## 4. Planar Jacobian checkpoint: what symmetry actually proves

Inherited from THM-4727/THM-1345, with a short independent derivation: for w>=1, a polynomial map equivariant for weights (1,-w) has F=(x A(t), y B(t)), t=x^w y. A Keller map has nonzero A,B. In orbit coordinates its quotient is phi(t)=t A(t)^w B(t), and det JF=A(t)^(1-w) phi'(t). If det JF=c nonzero then phi'=c A^(w-1). Degrees yield w deg A+deg B=(w-1)deg A, hence deg A+deg B=0. Thus A,B are constants and F is linear.

This covers EXACT equivariance, not arbitrary planar Keller maps, arbitrary common-frame invariants, or finite congruence symmetry. Finite cyclic symmetry only imposes charges modulo L, permitting extra terms. For L>=3, example F=(x+y^(L-1),y) is nonlinear, has determinant 1, and is C_L-equivariant for weights (1,-1), but is not C*-equivariant. For L=2 the same shear is linear but still not C*-equivariant. These maps are automorphisms, not counterexamples. This gives a concrete warning: replacing the infinite grading by mod-8 or mod-9 guards loses precisely the rigidity needed in the degree proof.

## Outcome

The useful common mechanism is preservation of a grading before contraction, not an equality of the objects. Finite observer frames are exactly computable at bounded bandwidth; infinite Gaussian power questions require an unbounded charge ledger or a proved finite reduction. Planar torus-equivariant rigidity likewise needs integer weights, not merely their residues. These are compatible with native port-preserving compilers but do not by themselves supply native tournament wiring or solve JC2.
