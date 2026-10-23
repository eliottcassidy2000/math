# Equal closed simplicial amplitudes can hide opposite boundary positivity

Status: PROVED finite-group construction and exact rational certificate. This is an explicit obstruction to reconstructing boundary/time structure from closed partition numbers. It is not a claim about the attached papers or a solution of four-dimensional Yang–Mills theory.

## 1. Two strictly positive S3 face weights

Use normalized counting Haar measure on G=S3. Its irreducible characters are1, the sign character s, and the two-dimensional standard character chi, whose values on identity/transpositions/three-cycles are(2,0,−1). Define

    w_plus=1+(1/4)s+(1/4)chi,
    w_minus=1−(1/4)s+(1/4)chi.

Their three class values are respectively

    (7/4,3/4,1),   (5/4,5/4,1/2).

Both are strictly positive, inversion-symmetric, conjugation-invariant densities with Haar integral1. In the convention w=Σ_r d_r a_r chi_r, their coefficients are

    (a_triv,a_sign,a_std)=(1,±1/4,1/8).

## 2. Every closed orientable triangulated surface gives the same scalar

For a finite connected closed orientable triangulated surface with V vertices, E edges, F faces and genus g, integrate the product of these face weights over independent normalized edge Haar variables. Character expansion and matrix-element orthogonality give

    Z_g,F(w)=Σ_r d_r^(V−E+F) a_r^F
            =Σ_r d_r^(2−2g) a_r^F.                      (1)

Proof: each face contributes d_r a_r. Integration at an edge forces its two neighboring representation labels to coincide and contributes1/d_r. Connectedness forces a single common representation label, while the remaining index loop at each vertex contributes d_r. The Euler exponent follows. This is the elementary finite-group version of the classical two-dimensional character method; the accompanying tetrahedron note proves its sphere case separately.

Every closed triangulated surface has3F=2E, so F is even. Consequently our two weights give exactly the same number on EVERY such mesh and every genus:

    Z_g,F=1+4^(−F)+2^(2−2g)8^(−F).                      (2)

In particular the tetrahedral sphere has Z=1029/1024. The checker independently enumerates its216 tree-gauge assignments for each weight. Equality is not limited to a finite tested list: evenness of F and (1) prove it for all closed orientable triangulated surfaces.

This is not an assertion that the models are topological invariants or subdivision-consistent: their values depend on F. It is equality of the two assignments for each identical mesh.

## 3. Their reflected-boundary kernels have opposite positivity

The boundary transfer operator is

    (K_w f)(x)=(1/6)Σ_y w(x^(-1)y)f(y).

For real inverse-symmetric w this operator is self-adjoint. Its eigenvalue on each irreducible matrix-coefficient sector is a_r, with multiplicity d_r². Thus

    spec(K_plus)={1, 1/4, 1/8 [multiplicity4]},
    spec(K_minus)={1,−1/4, 1/8 [multiplicity4]}.

For the normalized sign function, the reflected quadratic form is

    (1/36)Σ_(x,y) s(x) w_minus(x^(-1)y)s(y)=−1/4.

Hence K_minus is not positive semidefinite, despite every density value being positive. K_plus is positive definite. This is failure of this explicitly specified one-step reflection/cylinder positivity test, not a claim about every possible reflection or a two-step coarse theory. Squaring K_minus removes its negative eigenvalue but changes the time resolution.

## 4. A sharper time-refinement distinction

K_plus lies in a continuous symmetric Markov convolution semigroup. Let each of the three transpositions occur as a jump with rate(log2)/3, and each of the two nonidentity three-cycles with rate2(log2)/3. The generator has spectral values0,−2log2,−3log2, so its time-one exponential has exactly the coefficients1,1/4,1/8 above. The explicit generator has nonnegative off-diagonal entries and zero row sums.

By contrast det(K_minus)=−(1/4)(1/8)^4<0. For every finite real matrix Q, det(exp Q)=exp(tr Q)>0. Thus K_minus cannot be the time-one transition of any finite-dimensional time-homogeneous real continuous semigroup, even without requiring self-adjointness. In particular it has no reflection-symmetric continuous-time refinement on this same six-state space.

This does not rule out a model with added hidden states, a different observation map, or a different time-step convention. Those are new carriers with additional data, not reconstructions from the closed numbers alone.

## 5. Why this connection matters

The representation coefficients act as diagonal convolution channels, not unique multiplicative primes. An even number of faces erases one channel's sign from every closed scalar amplitude. That sign remains decisive for boundary positivity and time refinement.

This is the gauge-theory counterpart of Session33's electrical and port examples: closed scalar summaries can agree while the data needed for lawful gluing or physical interpretation differ. It explains why a simplicial partition identity, however elegant, must retain boundary operators and positivity before it can be used in a constructive quantum-field argument.

Literature context: Witten, [On quantum gauge theories in two dimensions](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-141/issue-1/On-quantum-gauge-theories-in-two-dimensions/cmp/1104248198.pdf), discusses the two-dimensional character framework. The finite example, normalization and proof above are self-contained; no infinite-group interchange or continuum construction is used.

Reproduction: `python3 check_boundary_positivity.py > boundary_positivity.json`, also under `-O`. Exact rational matrices, direct tetrahedral sums and generator spectral projectors supply the finite checks.
