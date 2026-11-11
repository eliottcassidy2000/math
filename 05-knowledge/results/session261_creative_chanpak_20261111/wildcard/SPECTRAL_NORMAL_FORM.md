# Universal integer spectrum of fixed-apex spine transfers

**PROVED elementary theorem.** Every periodic block in the fixed-k-apex grammar has exponential bases among 1^L,2^L,...,(k+1)^L, where L is the block length. In particular irrational or complex spectral bases cannot arise by changing the periodic word or the core tournament.

## 1. Common local flags

Use the positive automaton from `PROOF.md`, in column-vector convention. Fix one permutation of the core. For each of the first k segments, the local status space has standard basis e=empty, g=legal last endpoint, b=illegal last endpoint. Reading a letter and assigning its vertex to this segment acts by a reset operator R:

* R(e) is0 if the first-endpoint test fails, otherwise either g or b.
* R(g)=R(b), equal to g or b according to the new last-endpoint test.

Put d=g-b. Then R(d)=0 for every letter. The common flag

    span(d) subset span(d,g) subset span(d,g,e)

is invariant under every reset operator, and their diagonal in the ordered basis(d,g,e) is always(0,1,0). The first segment simply has its first-endpoint test always true; its flag and diagonal are the same.

For the final segment, the two states are empty e and nonempty g. A reset has R(g)=g and R(e) equal to0 or g. Its common flag span(g) subset span(g,e) has diagonal(1,0).

The transition for reading a letter is the sum over the k+1 choices of assigned segment, each acting as the corresponding local reset on that tensor factor and as identity on all other factors. A tensor-product ordering of the local flags simultaneously triangularizes these sums. Their diagonals are independent of the letter: each is a sum of k+1 zeros and ones.

The different core permutations are direct-sum blocks. Core orientation affects acceptance, not this conclusion.

## 2. Exact eigenvalue multiplicities and block spectra

Let c_j be the coefficient of x^j in

    (2+x)^k(1+x).

For each fixed permutation, every single-letter matrix has eigenvalue j with algebraic multiplicity c_j. Over all permutations the multiplicity is k! c_j.

For any nonempty word block u of length L, its matrix M_u is a product of simultaneously triangular matrices. Its diagonal is therefore j^L at the same positions. Hence its characteristic polynomial is exactly

    product_(j=0)^(k+1) (t-j^L)^(k! c_j).

For k=2 the per-permutation multiplicities are4,8,5,1 at eigenvalues0,1,2,3. Thus the36-dimensional block characteristic polynomial is

    t^8 (t-1)^16 (t-2^L)^10 (t-3^L)^2.

This statement concerns the full representation, not the minimal recurrence seen by a particular seed/output. Cancellations can remove eigenvalues and lower multiplicities; they cannot introduce new bases.

For any prefix p and suffix q, the sequence H(T_A(p u^m q)) therefore has, after a finite initial transient, the form

    sum_(j=1)^(k+1) P_j(m) j^(L m),

with rational-coefficient polynomials P_j. A universal degree bound is deg P_j<k!c_j. The0-eigenspace contributes only finitely many initial terms. This is an exact formula-shape theorem, not an asymptotic claim.

## 3. A smaller uniform annihilator

Each tensor factor has a flag of length3, except the last of length2. A strictly flag-lowering map decreases the total tensor grade by at least1; the largest grade is2k+1. Thus any product of2k+2 strictly flag-lowering maps vanishes.

For a block matrix M=M_u, put

    p(t)=product_(j=0)^(k+1)(t-j^L).

The matrix p(M) has zero diagonal in the common tensor flag and is strictly flag-lowering. Therefore

    p(M)^(2k+2)=0.

This gives a uniform annihilating polynomial of degree(k+2)(2k+2), independent of k! and independent of the word and core. It implies the alternative polynomial-degree bound deg P_j<=2k+1 and transient length at most2k+2. One may use the better of this bound and the multiplicity bound above.

The proof uses the tensor partial order, not the assertion that every strictly triangular D-dimensional matrix has nilpotence index2k+2. The transition matrices and all their products are triangular for this specific tensor filtration.

## 4. Consequences for creative constructor search

* The bases9,4,1 in the two-letter periodic families are forced by k=2 and L=2. They are not evidence of a new algebraic-number carrier.
* Fibonacci/Pell/golden-mean growth cannot appear as the exact periodic-word growth in this fixed-k grammar. To transfer those mechanisms physically, change the core with each step, use state-dependent word selection, or use a genuinely different operation.
* This theorem strengthens a finite-ray fit into a universal spectral restriction, while still allowing rich nonperiodic arithmetic and integer polynomial factors.
* The parent's cycle-packing bound independently limits the number of root-module visits in a single clone ray. The two bounds agree in spirit but concern different operations; a length-L block may have base j^L rather than j.

No target-coverage conclusion follows from the allowed spectrum alone.

## 5. Sharper Jordan bounds from the tensor partial order

The parent's partial-order observation sharpens the preceding safe annihilator. Use the equivalent reversed local basis(empty,legal,legal-illegal), with grades0,1,2 and diagonal labels0,1,0; the last factor has grades0,1 and labels0,1. Matrix entries move only upward in the product partial order.

Let a tensor basis position have diagonal label j, the number of grade1 factors, and total grade g. Then g-j is even. For j=0 its possible grades lie between0 and2k. For j>=1 they lie between j and2k-j+2. Therefore a strictly increasing chain of distinct positions with the same label has length at most

    e_0=k+1,
    e_j=k-j+2,       1<=j<=k+1.

To turn this positional fact into a Jordan bound, let P_j be the spectral projector of a fixed block matrix M onto its generalized j^L eigenspace. It is a polynomial in M, hence belongs to the same incidence algebra of the tensor partial order. For each original basis vector e_s of diagonal label j, replace it by P_j e_s. These columns form a unit-diagonal triangular change of basis S in that incidence algebra; S inverse belongs to it as well. In the new basis M connects only positions of the same eigenvalue and only comparable positions. Consequently each nilpotent Jordan chain has at most e_j steps, by the grade bound above. Direct sums over core permutations do not increase the bound.

Thus the improved universal annihilating polynomial is

    t^(k+1) * product_(j=1)^(k+1)(t-j^L)^(k-j+2).

Its degree is (k+1)(k+4)/2. In particular, after a transient of at most k+1 terms,

    H(T_A(p u^m q)) = sum_(j=1)^(k+1) P_j(m) j^(L m),
    deg P_j <= k-j+1.

The largest exponential base has a constant coefficient, the next has at most a linear coefficient, and so on. For k=2 every periodic word has the exact eventual shape

    c*3^(Lm) + (a*m+b)*2^(Lm) + d*m^2+e*m+f,

with rational constants, after at most3 initial terms. A universal degree9 annihilator is

    t^3(t-1)^3(t-2^L)^2(t-3^L).

`check_spectrum.py` checks both the earlier safe annihilator and this sharper one on all letters and the selected one-/two-letter blocks for k=0,1,2,3. These controls corroborate the argument; the all-k proof is the incidence-algebra change of basis above.
