# An arithmetic obstruction to self-dual Gaussian Hessian degeneracy

Status: PROVED using the cited Nesterenko algebraic-independence theorem, the inherited exact self-dual scalar identity, and an elementary rational-coordinate argument. This result does not need the new interval minimum certificate.

## Theorem

For every extremal even unimodular lattice of rank88, the determinant-one Gaussian Hessian at the self-dual temperature t=1 is nonsingular. Thus no actual angular eigenvalue realizes the abstract double-zero-at1 case in the scalar branch-sign diagram.

Nonsingular does not mean positive definite. The theorem removes an exact degeneracy, not the possibility of a saddle.

## Proof

1. The angular operator is defined over the rationals.

Choose a lattice basis with integral Gram matrix G. Model traceless Euclidean self-adjoint endomorphisms using

V_G={H in Mat_88(Q): H^T G=GH, tr(H)=0}.

Its realification is isomorphic to Sym_0(R^88), and the positive definite inner product is tr(HK). Choose any rational basis H_i. For each minimal lattice coordinate vector z in Z^88, its deformation functional is z^T G H_i z, a rational number. Therefore the matrices

S_ij=tr(H_i H_j),
B_ij=sum_{z minimal}(z^T G H_i z)(z^T G H_j z)-2723200*S_ij

are rational; S is invertible. The angular operator matrix is S^-1 B, so its characteristic polynomial belongs to Q[X]. In particular every eigenvalue of A is algebraic. The orthonormal real coordinate change does not alter those eigenvalues.

2. The self-dual scalar is transcendental.

Nesterenko's Theorem1 says that for any complex q with0<|q|<1, at least three of q,E2(q),E4(q),E6(q) are algebraically independent over Q. Here Eisenstein series are written in the q-variable. At q=exp(-2*pi), equivalently tau=i, E6(i)=0. Therefore the other three numbers are algebraically independent; in particular E2(i) and E4(i) are algebraically independent.

Using E2(i)=3/pi, put

y=4*pi^2 E4(i)=36 E4(i)/E2(i)^2.

This y is transcendental. Otherwise its minimal polynomial, after substituting36E4/E2^2 and clearing denominators, would give a nonzero rational polynomial relation between the independent E2(i),E4(i).

The Session28 exact reduction gives

r(1)=8192*(304325-17489844/y).

This is a nonconstant rational fractional-linear function of y and hence is transcendental.

3. Algebraic spectrum and transcendental scalar cannot cancel.

For every eigenvalue lambda of A, lambda+r(1) is nonzero, since lambda is algebraic and r(1) transcendental. Since b(1)>0, all Hessian eigenvalues b(1)(lambda+r(1)) are nonzero. This proves the theorem.

## Consequences and scope boundaries

* The self-dual point lies inside a fixed-sign interval for every actual angular branch, although it is a local maximum of the universal scalar ratio.
* Any actual Hessian degeneracy occurs away from t=1 and therefore has a distinct reciprocal partner.
* Combined with the scalar-shape theorem, a particular angular eigenvalue has0,2, or4 distinct zero temperatures. The abstract c=M three-zero case cannot occur for an integral lattice. The two-zero tangent case c=m is not excluded: no transcendence claim for m has been proved here.
* Gamma(1/4) is unnecessary for the proof. The classical identity y=3*Gamma(1/4)^8/(16*pi^4) supplies a complex-multiplication interpretation, but the direct Eisenstein argument avoids needing a separate special-value theorem.

The source-target map is precise: the finite rational shell operator supplies algebraic eigenvalues; modular special values supply a transcendental scalar translation. The preserved predicate is whether zero is in the translated spectrum. No equation between geometrically different spectra, no inference of positivity, and no assumed equality of lattices is used.

## Primary source verification

Yu. V. Nesterenko, Modular functions and transcendence questions, Sbornik: Mathematics187:9(1996),1319-1348, Theorem1 on printed p1320. DOI: https://doi.org/10.1070/SM1996v187n09ABEH000158 .

Original English PDF: https://www.mathnet.ru/links/f3406fc7df8e07b1490855a647fef28b/sm158_eng.pdf . Alternate indexed copy: https://www.mathnet.ru/links/438f6434f3a44bc6d872ed416438b2cd/sm158_eng.pdf .

The statement was verified from the primary PDF's indexed p1320 text on2026-10-22. Direct PDF open timed out, but the primary-source search extraction exposed the theorem statement and its quantifiers; no secondary paraphrase is needed as the mathematical dependency. The external theorem is CITED, not reproved or formally verified in this packet.

## Generalization: every extremal dimension congruent to16 modulo24

The same nondegeneracy theorem holds for every actual extremal even unimodular lattice of rank n=24m+16, m>=0. This is independent of whether its minimum shell is perfect.

Put k=n/2=12m+8 and a=k/4=3m+2. The degree2 harmonic theta series vanishes because, after division by Delta^(m+1), its weight is-2. Degree4 harmonic theta is a constant multiple of Delta^(m+1), its quotient having weight0. Thus the determinant-one Gaussian Hessian has the same fixed-operator form

Hess_t=b(t)(r(t)I+A), b(t)=(pi*t)^2 Delta(it)^(m+1),

where A is defined over Q by the first nonzero shell after subtracting its rational spherical average. In particular its spectrum is algebraic, by the same Gram-coordinate proof.

Write the rational modular theta polynomial as

F=c0 E4^a+c1 E4^(a-3)E6^2+terms divisible by E6^4.

For m=0 the E6^2 monomial is absent and c1=0. At tau=i, E6=0 and E4>0, hence c0=F(i)/E4(i)^a>0. Set C=5*a*c0/36+c1/2. Ramanujan differentiation gives, exactly,

D^2 F(i)=k(k+1)c0 E4(i)^a/(16*pi^2)+C E4(i)^(a+1),

DF(i)=k F(i)/(4*pi).

For clarity, the two contributing algebraic identities before substitution are

D^2(E4^a)|E6=0 =a(4a+1)E2^2 E4^a/36+5a E4^(a+1)/36,

D^2(E4^(a-3)E6^2)|E6=0 =E4^(a+1)/2.

Meanwhile

r(t)=[2D^2F/(k(k+1))-DF/(k*pi*t)]/Delta^(m+1).

Using Delta(i)=E4(i)^3/1728 and a=3m+2 yields

r(1)=2*C*1728^(m+1)/(k(k+1)) - c0*1728^(m+1)/(8*pi^2*E4(i)).

The second coefficient is nonzero. Therefore r(1) is a nonconstant rational affine function of1/(pi^2 E4(i)), hence transcendental by the same direct Nesterenko argument. No algebraic eigenvalue of A can cancel it. The self-dual determinant-one Gaussian Hessian is nonsingular in every such dimension for which an extremal lattice exists.

This general theorem concerns only self-dual nondegeneracy. The exactly-three-critical-temperature theorem and numerical minimum interval in REPORT.md remain specific to rank88; they have not been extended to other dimensions.
