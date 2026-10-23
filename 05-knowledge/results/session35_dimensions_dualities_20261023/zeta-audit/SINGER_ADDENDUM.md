# Independent addendum: Singer augmentation, Heisenberg construction, and E7 obstruction

Verdict: PASS in characteristic zero, with one important wording qualification: this constructs a finite-group action on a Heisenberg Lie algebra (or group), not an embedding of C57 into the unipotent Heisenberg group itself.

## Positive construction

Work over Q, R, or C. Index the standard basis of V=k^57 by Z/57. Let P e_j=e_{j+1}, F e_j=e_{7j}, and W={v:sum v_j=0}. Then dim W=56 and V=k·1 direct-sum W. The operators P and F preserve W; F^3=I and FPF^{-1}=P^7. They give a faithful action of C57 semidirect C3, order171.

The simple alternating form x^T(P-P^{-1})y is nondegenerate on W: P has distinct 57th-root eigenvalues, and lambda-lambda^{-1}=0 only at lambda=1 since57 is odd. It is P-invariant but is not generally F-invariant.

The improved form uses
Omega=(P-P^56)+(P^7-P^50)+(P^49-P^8).
It is skew-symmetric and invariant under P and F. Its symbol is
f(x)=x+x^7+x^49-x^56-x^50-x^8.
Exact rational polynomial gcd(f(x),x^57-1)=x-1 shows its sole ambient kernel is k·1. Omega maps W to W, so its restriction is nondegenerate. Independently verified this by both polynomial arithmetic and the determinant of B^T Omega B for B=(e_0-e_56,...,e_55-e_56). Exact matrix tests also verify P^T Omega P=Omega and F^T Omega F=Omega.

Define heis(W)=W direct-sum k z with [u,v]=(u^T Omega v)z and [z,u]=0. Jacobi follows because every bracket is central, and nondegeneracy gives center exactly k z. Dimension is57. The finite group acts by Lie automorphisms, fixing z. Thus the precise maps are
C57 semidirect C3 -> Sp(W,Omega) -> Aut(heis(W)).
Equivalently the same linear action preserves the characteristic-zero Heisenberg group multiplication (u,a)(v,b)=(u+v,a+b+Omega(u,v)/2).

Do not abbreviate this as a nontrivial group homomorphism C57->H57: this unipotent group has no nontrivial finite-order elements, since (u,a)^m=(mu,ma). A combinatorial-to-Lie construction called “Singer57 to H57” is acceptable only with its category and retained action stated.

## Independent E7 obstruction

Let h=P^19, of order3. On V it has nineteen3-cycles, so fixed-space dimension19 and trace0. Removing the constant line gives dim W^h=18 and trace_W(h)=-1. Its eigenvalue multiplicities are (18,19,19) at (1,omega,omega^2).

Let G be the simply connected complex group E7 with its minuscule56-dimensional representation. Its maximal-rank A7 subgroup is SL8/{+I,-I}, with restriction of the56 equal to exterior-square(k^8) direct-sum its dual. Any order3 element of G is semisimple and conjugate into a maximal torus. Since the A7 subgroup has full rank, every such conjugacy class has a representative in it. Its order3 representative has an order3 lift to SL8: if a lift has cube -I, negate it. Negation preserves determinant because8 is even and does not alter its quotient image.

For that lift, write eigenvalue multiplicities (a,b,c) at (1,omega,omega^2). Then a+b+c=8 and b+2c=0 modulo3. The fixed dimension on exterior-square is binom(a,2)+bc. Doubling for the dual gives d=a(a-1)+2bc. Exact exhaustive enumeration of these small nonnegative triples yields
d in {2,14,20,32,56}.
The identity gives56. The nonidentity order3 traces are (3d-56)/2, namely {-25,-7,2,20}. In particular18 and trace-1 are impossible.

Therefore no homomorphism C57->E7(C) can make the minuscule representation restrict to this distinguished augmentation module W. Changing the symplectic form cannot repair the character mismatch. Consequently the full C57 semidirect C3 action cannot extend through that E7 action either. The same obstruction applies to real forms by complexification. This statement does not exclude abstract C57 subgroups of E7, different56-dimensional C57 modules, or an unmarked abstract isomorphism of heis(W) with another57-dimensional Heisenberg Lie algebra, including a suitable parabolic nilradical.

## Primary-source verification

Griess and Ryba, “Classification of finite quasisimple groups which embed in exceptional algebraic groups,” coauthor-uploaded full text:
https://www.researchgate.net/publication/243106241_Classification_of_finite_quasisimple_groups_which_embed_in_exceptional_algebraic_groups

The G2(3) discussion on printed page39 explicitly lists {-7,-25,2,20} as the possible order3 traces on the minuscule56, referring to [CG]. It corroborates the independent maximal-torus/A7 derivation above. The author-hosted URL supplied by the parent currently redirects to an unavailable path; the coauthor-uploaded full text was actually opened and its relevant passage inspected. No conclusion rests solely on a search snippet.

## Reproduction

Run python3 session35-zeta-audit/check_singer_addendum.py. Output is check_singer_addendum.out; every assertion passed. These are exact rational/integer matrix and polynomial tests, not floating-point tests. No Git operations or repository edits were performed.
