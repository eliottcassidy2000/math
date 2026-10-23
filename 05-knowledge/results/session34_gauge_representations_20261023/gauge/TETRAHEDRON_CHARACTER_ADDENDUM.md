# Irreducible blocks and tetrahedral gauge gluing

Status: PROVED finite-group character formula, with an explicit finite S3 certificate. Independently derived from elementary character orthogonality. No4-dimensional Yang–Mills or vertex-algebra consequence is asserted.

Let G be any finite group, with normalized averaging dg=|G|^(-1)sum_g. Let w:G->R be a class function satisfying w(g^-1)=w(g). Write its irreducible-character expansion

w(g)=sum_(rho) d_rho a_rho chi_rho(g),
a_rho=(1/d_rho) integral w(g) conjugate(chi_rho(g)) dg.        (T1)

All a_rho are real, and dual representations have equal coefficients. The normalization d_rho is essential. Orthogonality gives central convolution

chi_rho * chi_sigma = delta_(rho,sigma) chi_rho/d_rho,
(f*h)(g)=integral f(x) h(x^-1 g) dx.                          (T2)

Thus w^(convolution k)=sum_rho d_rho a_rho^k chi_rho. Irreducible blocks diagonalize gluing: a_rho multiply. This is an actual algebraic content for the user's 'representations analogous to primes' intuition; it is not unique prime factorization, and tensor products of irreducibles still decompose rather than remaining irreducible.

**Tetrahedral partition theorem.** Give each of the six unoriented edges of a tetrahedron a G link, with reversed links inverse. Form four face holonomies. Then the normalized-Haar partition function is

Z_tetra(w)=integral_(G^6) product_(four faces) w(F_face) dU
         =sum_rho d_rho^2 a_rho^4.                           (T3)

Choice of face orientation is harmless because w is inverse symmetric. Gauge-fix the three root edges to1 using a root-fixed vertex gauge transformation. Normalized averaging loses no factor. The three residual links a,b,c give faces a,b,c,abc^-1, up to face orientations. Therefore

Z=integral w(a)w(b)w(c)w(abc^-1) da db dc.

Since w is central, integrating c gives (w*w)(ab). Integrating a,b against the inverse-symmetric w gives (w*w*w*w)(1). Formula(T2), and chi_rho(1)=d_rho, prove(T3).

Equivalent geometric interpretation: four triangular faces form the boundary2-sphere of a tetrahedron. This is a finite2-dimensional surface gluing calculation, not a four-dimensional bulk construction. Four here counts faces, not edge degree or Fano Walsh degree. It agrees with the main note's Z2 surface expansion: for G=Z2 and w(g)=exp(beta g), a_trivial=cosh(beta), a_sign=sinh(beta), so normalized Z=cosh(beta)^4+sinh(beta)^4.

**Closed exact control for any finite group.** Put q=|G| and w(g)=z if g=1, otherwise1. Then

a_trivial=(z+q-1)/q; a_rho=(z-1)/q for every nontrivial rho.

Since sum_rho d_rho^2=q,

Z_tetra(z)=[(z+q-1)^4+(q-1)(z-1)^4]/q^4.                    (T4)

After root gauge fixing, the UNNORMALIZED sum over q^3 assignments is

z^4+6(q-1)z^2+4(q-1)(q-2)z+(q-1)(q^2-3q+3).                (T5)

For S3, this is z^4+30z^2+80z+105. In particular three trivial faces force the fourth trivial, giving the absent z^3 term; the identity-edge assignment gives the sole rooted fully flat configuration. The main note's full216 flat configurations on all six S3 edges account instead for the root-fixed gauge group size6^3, so there is no count discrepancy.

Certificate: tetrahedron_characters.py checks the full character convolution table on all group arguments, all216 rooted S3 assignments, and three arbitrary real central weights including negative values. Normal and optimized outputs agree, with58 explicit checks. The polynomial result is exact, not sampled in z. No external paper statement is required for the derivation; a literature-priority claim is not made.

For compact groups the same finite-character-expansion calculation holds termwise. Infinite expansions require convergence/integrability justification before summing or gluing; this addendum does not silently assume it. Heat-kernel specializations need an explicitly chosen Laplacian and normalization, neither supplied here. No PDF or attachment was read.
