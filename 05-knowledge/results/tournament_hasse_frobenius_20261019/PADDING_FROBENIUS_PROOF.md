# Deletion, Hasse digits, and an exact p-complex on tournament functions

2026-10-19. **PROVED elementary identities and integral classification; independently reviewed.** This is a native-operation extension of the characteristic-two padding complex. It makes no external-priority claim and proves no Rule-30 or Gilbreath conjecture.

## Inheritance and concept board

Closest proved mechanism: `05-knowledge/results/tournament_padding_complex_20261019/PADDING_COMPLEX.md`, especially the source-adjunction contraction. Closest cellular-automaton carrier: THM-4210, `01-canon/theorems/THM-4210-rule30-lossless-dyadic-block-current-cartier-tree.md`, where Frobenius/Cartier transport needs a physical quadratic-current sidecar. Canonical hostile: THM-4210's ambient even-channel kernel and arbitrary post-prefix forcing; algebraic transport alone cannot decide a fixed-seed trace. Corrected near miss: the odd Cartier current is center-null directly but necessary for physical nonlinear current generation. Least-used sidecar pursued here: unordered deletion multiplicity before reduction, rather than only repeated one-vertex deletion.

Concept board: (1) deletion/source Weyl relation; (2) characteristic-p nilpotence; (3) Hasse p-adic digits and carries; (4) automorphism-invisible embedding counts; (5) physical versus ambient cellular-automaton carriers. The actual transfer is between (1)--(4); (5) is a sharp limitation, not a claimed dynamical equivalence.

## 1. Definitions, grading, and integral identities

Let C_n(R) be the R-module of isomorphism-invariant functions on n-vertex tournaments, with the unique empty tournament included at n=0. Set C_n=0 for n<0, and work on the graded direct sum C=direct_sum_(n>=0) C_n. All maps below retain the multiplicities of vertex subsets, including isomorphic deletions.

For m>=0 define the degree-m subset-deletion operator

    (D_m f)(T) = sum_{S subset V(T), |S|=m} f(T-S),   f in C_n, |T|=n+m.

Thus D_0=I and d=D_1. Define the degree-minus-one operator h by

    (h f)(T)=f(T+source),

where the new vertex points to every old vertex. Define h=0 on C_0. These operators are defined over Z before any reduction. For a,b>=0,

    D_a D_b = binom(a+b,a) D_(a+b),                  (1)
    h D_m - D_m h = D_(m-1)  (m>=1),               (2)
    [h,d]=I,                                       (3)
    d^m=m! D_m.                                    (4)

For (1), each deleted (a+b)-subset is partitioned into the two deletion batches in binom(a+b,a) ways. For (2), after adding the new source, split the deleted m-subsets according to whether they contain it. The terms not containing it are D_m h; the remaining terms are D_(m-1). This proof also handles bottom degrees, using the zero convention. Equation (4) follows by ordering the m deleted vertices.

More generally, induction on r gives

    h^r D_m = sum_(j=0)^min(r,m) binom(r,j) D_(m-j) h^(r-j).  (5)

No sign or automorphism factor is suppressed in these identities.

## 2. All-prime contraction and explicit primitives

Fix a prime p and work over F_p. From (4), d^p=0. Define H=h^(p-1). Then

    - sum_(j=0)^(p-1) d^(p-1-j) H d^j = I.          (6)

Proof: [h,d]=I gives the normal-order identity

    d^a h^b = sum_(i=0)^min(a,b)
        (-1)^i binom(a,i) (b)_i h^(b-i)d^(a-i),

where (b)_i=b(b-1)...(b-i+1). Apply it with b=p-1 and sum over a=0,...,p-1, adding the remaining power d^(p-1-a). The coefficient indexed by i is

    (-1)^i (p-1)_i sum_(a=0)^(p-1) binom(a,i)
      = (-1)^i (p-1)_i binom(p,i+1).

It vanishes modulo p for i<p-1. For i=p-1 the scalar is (-1)^(p-1)(p-1)!=-1 in F_p, including p=2. This proves (6).

For every 1<=k<p and every degree n,

    ker(d^k:C_n -> C_(n+k))
       = im(d^(p-k):C_(n-p+k) -> C_n).              (7)

Indeed, d^p=0 gives image containment. If d^k f=0, terms j>=k in (6) vanish, so the explicit primitive is

    g = - sum_(j=0)^(k-1) d^(k-1-j) h^(p-1) d^j f,
    d^(p-k)g=f.                                    (8)

Every term of g has degree n-p+k, so negative-degree cases are included. In particular, d^k is injective when n<p-k. For k=1, the primitive is simply -h^(p-1)f.

This is an exact p-complex with a displayed p-contracting homotopy. Ordinary exactness is the special case p=2. It is not a statement that a primitive is represented by a single tree or even a single oriented graph.

Equivalently, the nilpotent operator d on C has only Jordan strings of length p: no shorter nonzero string survives. One may omit this equivalent module-language consequence without affecting any proof above.

## 3. Every Hasse digit has the same exact contraction

For e>=0 set s=p^e and

    delta_e=D_s,      eta_e=h^s.

Then

    delta_e^p=0,               [eta_e,delta_e]=I.    (9)

The first equality follows from (1): the integer scalar multiplying D_(ps) in D_s^p is (ps)!/(s!)^p. Its p-adic valuation is one, because v_p((p^(e+1))!) - p v_p((p^e)!) = 1. For the second, set r=m=s in (5). All binomial coefficients binom(s,j), 0<j<s, vanish modulo p; the endpoints are D_s h^s and I.

Applying the algebraic proof of (6) to this pair gives

    -sum_(j=0)^(p-1) delta_e^(p-1-j) eta_e^(p-1) delta_e^j=I.  (10)

Therefore, for 1<=k<p,

    ker(delta_e^k:C_n -> C_(n+ks))
       = im(delta_e^(p-k):C_(n-(p-k)s) -> C_n),      (11)

with the explicit primitive

    g=-sum_(j=0)^(k-1) delta_e^(k-1-j)
          h^((p-1)s) delta_e^j f.                   (12)

This acts separately on each degree residue class modulo s. In particular delta_e^k is injective for n<(p-k)s. The all-digit theorem is stronger than merely noticing d^p=0.

### Digit reconstruction and the exact carry rule

Write m=sum_e m_e p^e, with 0<=m_e<p. Repeated use of (1), or Lucas's elementary binomial digit rule, gives

    D_m = product_e delta_e^(m_e)/(m_e!).            (13)

Each factorial in the denominator is invertible in F_p. A self-contained verification is to use

    (1+z)^m = product_e (1+z^(p^e))^(m_e) in F_p[z],

which follows by Frobenius; coefficient comparison proves the needed binomial digit rule. For a+b, equation (1) is nonzero precisely when addition of a and b in base p has no carry. When there is a carry, the binomial scalar is zero; without a carry it is the product of the local digit-binomial coefficients. The composition operator is nonzero exactly in the no-carry case, since D_(a+b) is nonzero on 1_0; this does not mean every input has nonzero output. The deletion composition support is therefore computed by a two-state carry transducer. This is the precise cellular-automaton/Pascal connection: it is an operation on actual subset-deletion functions, not an analogy based on a triangle picture.

The higher digits cannot be discarded. On the empty-tournament function 1_0,

    d^p 1_0=0,        D_p 1_0=1_p !=0.              (14)

Thus no positive iteration of d can recover all Hasse deletions. More generally the graded operator algebra generated by delta_0,...,delta_(e-1) has nonzero homogeneous degrees only from 0 through p^e-1, whereas delta_e has degree p^e and is nonzero on 1_0. Each new digit is genuinely new operator data.

The digit contractions are not a family of mutually independent Weyl pairs. From (5),

    [h^(p^e), D_(p^f)] =
      0                    if e>f,
      I                    if e=f,
      D_(p^f-p^e)          if e<f.                  (15)

The lower-triangular interaction is a borrow sidecar. Claiming independent digitwise source coordinates would lose (15).

## 4. Embedding observability: the special role of characteristic two

For an oriented graph S on n vertices, let phi_S(T) be its number of spanning injective embeddings into the n-vertex tournament T, reduced modulo p. Let B_n be the linear span of all such functions. Then

    B_n = {f in C_n : f(T)=0 whenever p divides |Aut(T)|}.  (16)

Proof of containment: Aut(T) acts freely on Emb(S,T) by postcomposition. A spanning embedding is a bijection, so only the identity automorphism can fix it. Therefore |Aut(T)| divides every embedding count. Conversely, if p does not divide |Aut(T)|, take S=T. Then phi_T is |Aut(T)| times the indicator of T's isomorphism class, giving that indicator after multiplication by an inverse scalar. These indicators span the right side of (16).

For p=2, every tournament automorphism group has odd order: an involution would exchange a pair and reverse its directed edge. Hence B_n=C_n, recovering the special full-function-space statement in the inherited padding theorem. For every odd p, the cyclic tournament on Z/pZ with arcs i->i+a for a=1,...,(p-1)/2 has a translation automorphism of order p, so B_p is proper. At p=3 the directed triangle is the smallest blind host. All its spanning embedding counts are zero modulo 3, while its class indicator is nonzero.

The ordinary pair d,h preserves B:

    d phi_S = phi_(S disjoint union one isolated vertex),
    h phi_S = sum_(u:indegree_S(u)=0) phi_(S-u).       (17)

The same direct image/preimage counting as in the inherited theorem proves these formulas over Z. Thus the p-contraction (6)--(8) restricts to B, and also descends to C/B. In particular its exactness is not an artifact of adding non-embedding functions.

But higher Hasse digits need not preserve B when p is odd. The witness is already (14): 1_0 belongs to B_0, whereas D_p1_0=1_p does not vanish on the cyclic p-vertex tournament, hence is outside B_p. More concretely,

    phi_(S disjoint union m isolated vertices) = m! D_m phi_S.  (18)

For m>=p, the right side is zero modulo p even when D_m phi_S is nonzero. Unordered subset-deletion before reduction retains precisely the information erased by labelled-isolate permutations. Division by m! after reduction is invalid.

### Integral embedding lattice: reduction image is not abstract reduction

Let E_n be the Z-lattice spanned by all phi_S before reducing their values. The same free postcomposition action and the source choice S=T prove exactly

    E_n=direct_sum_[T] |Aut(T)| Z 1_[T].             (18a)

Thus the image of E_n modulo p inside C_n(F_p) is B_n, but the abstract tensor product E_n tensor F_p need not inject into C_n(F_p). A basis vector |Aut(T)|1_[T] becomes invisible under this map precisely when p divides |Aut(T)|. The index [C_n(Z):E_n] is product_[T]|Aut(T)|.

This distinction affects the differential itself. The integral d and h preserve E, but d^p need not vanish after abstract reduction of E. For p=3,

    d^3 1_0 = 6 1_3.

Its coefficient at the directed-triangle basis vector 3·1_[C3] of E_3 is 2, nonzero modulo 3. In the ambient pointwise lattice C_3 the same function has all coefficients divisible by 3. Therefore the p-complex statements for B are statements about pointwise-reduced observables, not about E tensor F_p. The first failed implication in the naive transfer is that a nonsaturated lattice inclusion would remain injective modulo p.

## 5. Integral reconstruction basis and complete deletion cokernels

The preceding field statements lift to a stronger exact integral decomposition. Return to C_n(Z), put t_n=rank C_n, t_(-1)=0, and define

    K_n=ker(h:C_n -> C_(n-1)),       k_n=t_n-t_(n-1).

A nonempty tournament has at most one universal source. Source adjunction therefore bijects the (n-1)-vertex classes onto the n-vertex classes with a universal source. Consequently h is a split surjection, and K_n has the canonical basis of indicators of n-vertex tournaments with no universal source. The empty tournament is included as a source-free core at n=0; K_0=Z.

### Integral Taylor projector

Define the degree-preserving operator

    P=sum_(j>=0) (-1)^j D_j h^j.                    (19)

Each sum is finite on a homogeneous function, so no convergence or completion is involved. Equation (2) telescopes to hP=0, and P fixes ker h. Furthermore,

    f=sum_(j>=0) D_j P h^j f.                       (20)

Indeed, inserting (19) into the right side gives coefficient

    sum_(j=0)^ell (-1)^(ell-j) binom(ell,j)

on D_ell h^ell. This is one for ell=0 and zero otherwise. If a belongs to K_r, then h^i D_j a=D_(j-i)a for i<=j and is zero for i>j; also P D_l a=0 for l>0. Therefore

    P h^i D_j a = a if i=j, and 0 otherwise.         (21)

Equations (20)--(21) prove an integral direct sum, not merely a rational splitting:

    C_n(Z) = direct_sum_(r=0)^n D_(n-r) K_r.         (22)

The maps D_(n-r) are injective on K_r and the displayed basis is unimodular. In algebra notation C is a free graded module over the divided-power algebra, with primitive lattice K. No division by a factorial has occurred.

### The basis is counts of source-free induced subtournaments

For each source-free tournament A of order r<=n, let

    b_A^(n)(T)=#{U subset V(T): |U|=r, T[U] is isomorphic to A}.

These are exactly D_(n-r) applied to the indicators of A, so all the functions b_A^(n), over all source-free cores A of order at most n, form a Z-basis of C_n(Z). Thus every integer-valued isomorphism-invariant tournament statistic has a unique integral expansion in source-free induced-subtournament counts. The coefficients are explicit:

    f(T)=sum_A a_A(f) b_A^(n)(T),

    a_A(f)=sum_(S subset V(A)) (-1)^|S|
         f((A-S) with n-r+|S| successive universal sources).  (23)

The augmented host in every summand has exactly n vertices. This is a finite inclusion-exclusion reconstruction formula. It is not a claim that proper subtournament counts reconstruct every tournament: source-free cores of order n remain among the basis coordinates.

### Diagonal form and the integral cokernel

From (1), for a in K_r,

    D_m D_(n-r)a = binom(n-r+m,m) D_(n-r+m)a.        (24)

Thus in the unimodular bases (22), D_m:C_n(Z)->C_(n+m)(Z) has the rectangular diagonal presentation

    binom(n-r+m,m), repeated k_r times, r=0,...,n.   (25)

The unused codomain basis has cardinality t_(n+m)-t_n. In particular D_m is injective over Z and

    coker D_m is isomorphic to
      Z^(t_(n+m)-t_n)
        direct_sum [direct_sum_(r=0)^n
          (Z/binom(n-r+m,m)Z)^(k_r)].               (26)

The valuation of each entry at a prime p is exactly

    v_p binom(j+m,m)=sum_(a>=1)
      [floor((j+m)/p^a)-floor(j/p^a)-floor(m/p^a)],   j=n-r.

Each summand records a carry across the corresponding base-p digit boundary. This is the exact carry multiplicity, not merely the zero/nonzero Lucas test.

Formula (26) is a complete integral torsion classification. The diagonal entries in (25) need not divide one another, so they are not asserted to be the Smith invariant factors in that order. To obtain the actual invariant factors, for each prime sort the valuations of all t_n diagonal entries increasingly, including zeros, and recombine the primes position by position. Equivalently, compute successive determinantal divisors of the diagonal list. The finite audit also computes these divisors directly from the original deletion matrices.

For ordinary deletion d_n=D_1, the diagonal is n-r+1 repeated k_r. For example, t_0,...,t_5=1,1,1,2,4,12 and k_0,...,k_5=1,0,0,1,2,8. Thus d_4 has diagonal 5,2,1,1, actual Smith invariant factors 1,1,1,10, and

    coker d_4 = Z^8 direct_sum Z/10Z.               (27)

More generally, for n,m>=1 the image of D_m is never saturated: the r=0 summand in (26) is Z/binom(n+m,m)Z and is nontrivial. For n=0 the image is saturated. This is the exact integral obstruction that fieldwise exactness alone does not reveal.

### Finite-field ranks and string multiplicities

Reduction of the unimodular bases gives, over F_p,

    rank D_m on C_n = sum_(r=0)^n k_r
          * indicator[p does not divide binom(n-r+m,m)].      (28)

Lucas carries therefore classify every rank loss, not only nilpotence. When 1<=m<p, D_m is injective in all degrees n<p-m and first fails injectivity at n=p-m, with the constant function furnishing a kernel vector. This last threshold is only for 1<=m<p; equation (28) is the all-m statement.

For d, each divided-power string attached to K_r splits modulo p into blocks of p successive indices j=lp,...,lp+p-1. All intermediate coefficients are nonzero and the final coefficient vanishes. Hence the number of p-strings starting in degree n is exactly

    b_n=sum_(l>=0) k_(n-lp)
       =sum_(l>=0) [t_(n-lp)-t_(n-lp-1)],          (29)

where negative subscripts mean zero. For T(z)=sum_n t_n z^n,

    sum_n b_n z^n = T(z)/(1+z+...+z^(p-1)).         (30)

For the digit delta_e=D_s, s=p^e, the corresponding p-string starts have multiplicity

    b_n^(e)=sum_(l>=0)[t_(n-lps)-t_(n-s-lps)],
    sum_n b_n^(e)z^n=T(z)/(1+z^s+...+z^((p-1)s)).  (31)

To see this directly, (24) and Lucas show that delta_e on D_j K_r multiplies by the e-th base-p digit of j plus one. A string starts when that digit is zero, has p members spaced s apart, and ends when that digit reaches p-1. All other digits remain fixed within the string. Formulas (29)--(31) are consequences of the integral decomposition, not separate conjectural enumerations.

## 6. Scope and exact connection contract

Source: tournament functions with literal vertex-subset deletion and canonical source adjunction. Target: an explicit p-contractible graded module, with Hasse digit generators and a carry transducer. Preserved: every function value, deletion multiplicity, source operation, kernel/image equality, and the embedding-visible subspace for the ordinary pair. Destroyed by ordinary mod-p iteration: deletion levels at multiples of p; restored by divided-power operators. Destroyed by spanning-embedding observables at odd p: hosts whose automorphism order is divisible by p; the class-indicator sidecar restores them. Higher Hasse operations need the full function carrier or an enlarged observable class.

No map from physical Rule-30 nonlinear currents to tournament functions is constructed. THM-4210's quadratic physical-admissibility constraint remains essential. No prime-distribution input, absolute-difference absorption theorem, modular-period relation, or fixed-seed asymptotic is implied by this p-complex. The transferable result is the native operator mechanism and its exact failure boundary.

## 7. Exact verification and source boundary

Run from this packet directory:

    python3 verify_padding_frobenius.py
    python3 -O verify_padding_frobenius.py

The normal and optimized receipts match byte-for-byte on 46,839 checks. The class counts at orders zero through five are 1,1,1,2,4,12. Integral subset/source identities retain all multiplicities. The contraction is checked on every class-basis function at orders zero through five for (p,e)=(2,0),(3,0),(5,0),(7,0),(2,1),(2,2),(3,1). There are 52 direct kernel/image rank cases where the target order stays at most five, and 12 direct nonzero-kernel basis primitive checks. The carry atlas checks all 0<=a,b<=64 for p=2,3,5,7. The source-embedding hostile exhausts all 27 oriented graphs of order three at the cyclic host. The integral audit constructs the source-free induced-count basis and its explicit two-sided inclusion-exclusion inverse through degree five, checks their determinants are units, and diagonalizes every positive D_m with source and target orders at most five. For all 15 such matrices it independently computes every determinantal divisor from the original matrix and checks the predicted Smith invariant factors. These bounded computations are controls for the proofs above, not their universal justification.

The framework of N-complexes and the homology groups ker(d^k)/im(d^(N-k)) is classical; see M. M. Kapranov, [On the q-analog of homological algebra](https://arxiv.org/html/q-alg/9611005v1), Definitions 0.1 and 1.1 (preprint v1, 1996; inspected 2026-10-19). This is a terminology/framework attribution only. The characteristic-p tournament and Hasse identities above are proved directly here, without invoking a theorem of that paper. No priority for the abstract framework, divided powers, or Lucas arithmetic is claimed.
