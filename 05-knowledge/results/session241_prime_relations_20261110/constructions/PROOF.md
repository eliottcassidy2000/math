# Faithful multiplicative arithmetic in actual tournament path counts

**Status: elementary all-order proofs, independently checkable candidate for canon; finite controls are separately labelled. Chan–Pak remains OPEN.**

## Inheritance, novelty boundary, and concept board

The count H is the number of directed Hamiltonian paths, including one empty path for the empty tournament. Sources recovered before construction:

- `05-knowledge/results/arithmetic_seams_20260921_tournaments.md`, section 5, already proves H(E_n)=1+2^(n−2) by a path bijection. `01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md` explicitly records Fermat path counts. Neither the gadget nor that observation is new here.
- `01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md` is the inherited SCC-product/prime-lane mechanism. Its product identity and missing 7/21 are not replaced by this construction.
- `01-canon/theorems/THM-5092-four-test-tail-content-and-invisible-reroot-escape.md` shows fixed-port profiles lose reroot behavior even when order is retained. Its states (9,9,3,3) demonstrate a real gcd obstruction and its possible escape, not arbitrary primality control.
- `01-canon/theorems/THM-5000-coprime-interface-composition-consecutive-odd-windows.md` composes coprime moving cores into long consecutive odd intervals. Its baseline/placement limitation remains: coprimality alone is not target selection.
- `01-canon/theorems/THM-4493-collatz-gcd-height-rank-certificate.md` already separates multiplicative independence from private primes, using (6,10,15). Here (15,51,85) is an odd, actually tournament-realizable analogue, and the private-prime condition characterizes divisibility reflection rather than independence.

The proposed synthesis is the **all-positive-integer multiplicative divisibility-lattice embedding into actual H-values**, its exact loss ledger, and the **sharp cost theorem within the inherited endpoint-reversal bank**. No literature-priority claim is made.

Concept board: (i) SCC multiplication versus additive graft increments; (ii) independent valuation columns versus private primes; (iii) Fermat coprimality versus numerical primality; (iv) induced containment versus scalar divisibility; (v) bank cost versus spectrum surjectivity. The tests below were chosen to separate these five distinctions.

## 1. Native tournaments and multiplication

For a positive integer a, let E_(a+2) have vertices 0,...,a+1, all arcs i→j for i<j except the endpoint arc a+1→0. It is strong, because the increasing path closes to a Hamiltonian cycle. A Hamiltonian path not using the backward arc is the unique increasing path. A path using it is uniquely

    (S increasingly), a+1, 0, (the remaining interior vertices increasingly),

for an arbitrary subset S of the a interior vertices. Hence

    H(E_(a+2)) = 2^a + 1.                                      (1)

This is exactly the inherited bijection, repeated to make the present proof self-contained.

For tournaments A,B, let A⇒B mean their disjoint union with every cross arc pointing from A to B. Every Hamiltonian path visits all of A before all of B and is uniquely a path of A concatenated with a path of B. Thus

    H(A⇒B)=H(A)H(B).                                           (2)

The empty tournament is the strict identity for this operation; its H is 1. If empty tournaments are disallowed, use any transitive tournament for the scalar encoding of 1, with no claim of a strict tournament-level identity.

## 2. Exact scalar criteria for any gadget bank

Let (w_i) be an arbitrary finite or countably infinite family of integers >1, and write F(e)=product_i w_i^(e_i) for finitely supported vectors e of nonnegative integers. The source uses coordinatewise divisibility e≤f, gcd=min, and lcm=max.

**Proposition A.**

1. F is injective iff the prime-valuation columns (v_q(w_i))_q are linearly independent over Q with respect to finite relations.
2. F reflects divisibility, meaning F(e)|F(f) iff e≤f, iff each w_i has a **private prime**: some prime dividing w_i that divides no w_j, j≠i.
3. F preserves gcd for every pair iff the w_i are pairwise coprime. The same equivalence holds for lcm preservation, and for preservation and reflection of coprimality.

**Proof.** A scalar equality is an integer relation among valuation columns. Split any such relation into its positive and negative parts; rational relations can be cleared of denominators. This proves (1).

For (2), a private prime reads its coordinate inequality from valuations. Conversely, suppose w_i has no private prime. For each of its finitely many prime divisors q, choose another w_(j(q)) divisible by q. Choose sufficiently large positive powers of these finitely many other weights so that their product is divisible by w_i. Then F(e_i) divides the image of a vector with i-coordinate zero, contradicting reflection. This argument covers infinite banks as well.

For (3), applying gcd preservation to two distinct unit vectors forces gcd(w_i,w_j)=F(0)=1. With pairwise coprimality, valuations belong to unique columns and min/max of exponent vectors commute with F. Lcm preservation similarly forces lcm(w_i,w_j)=w_iw_j. Coprimality preservation forces the same pairwise condition; under that condition disjoint source supports are equivalent to disjoint prime supports of their images. ∎

**Hostile.** Weights 15,51,85 are realized by the ordinal products of the gadgets of H-values (3,5), (3,17), and (5,17). Their valuation matrix at primes 3,5,17 is

    [1 1 0]
    [1 0 1]
    [0 1 1],

with determinant −2, so the scalar encoding is injective. No column has a private prime. Indeed 15 divides 51·85=4335 with quotient289, although the first source generator does not divide the product of the other two. Moreover gcd(15,51)=3 although the source generators are coprime. This is the precise failure of the implication “injective multiplicative encoding implies faithful divisibility.”

## 3. All positive integers embed into actual H-values

Let p_0=2,p_1=3,p_2=5,... be the ordinary primes in increasing order. Let

    F_i = 2^(2^i)+1,
    B_i = E_(2^i+2).

The Fermat identity product_(j<i) F_j=F_i−2 implies pairwise coprimality: any common divisor of F_i and an earlier F_j divides 2 and is odd, hence equals1. This does not require any F_i to be prime.

Given n=product_i p_i^(e_i), form T(n) by taking e_i copies of B_i, in increasing i order, with ordinal joins between all blocks. Use T(1) empty. Define Phi(n)=H(T(n)). Equations (1) and (2) give

    Phi(n) = product_i F_i^(e_i),
    |T(n)| = sum_i e_i(2^i+2).                                (3)

**Theorem B.** For all positive integers m,n:

    Phi(mn)=Phi(m)Phi(n),
    m|n iff Phi(m)|Phi(n),
    Phi(gcd(m,n))=gcd(Phi(m),Phi(n)),
    Phi(lcm(m,n))=lcm(Phi(m),Phi(n)),
    gcd(m,n)=1 iff gcd(Phi(m),Phi(n))=1.

In particular Phi is injective. These follow directly from Proposition A and Fermat coprimality. Thus the multiplicative monoid and its entire divisibility lattice are faithfully present inside the tournament H-spectrum, without assuming Chan–Pak.

**Decoding.** For a promised image value x, repeatedly divide by F_i to recover its exponent. Composite F_i cause no ambiguity because no prime divisor of one F_i occurs in another. Only finitely many i can satisfy F_i≤x, so this is a terminating exact procedure. To return n, the fixed prime-index dictionary i↔p_i is needed. No factorization of F_i is needed. This is a computability statement, not an efficient factoring or primality algorithm. Computing T(n) from n uses its prime factorization or a supplied factorization certificate.

**Tournament-level operation and quotient.** The raw ordered tournament T(m)⇒T(n) need not be isomorphic to T(mn): SCC words have different block order. Sorting the strong components by the fixed gadget index gives the canonical T(mn). At the scalar H level no sorting is needed. Multiplication is preserved exactly by H, not asserted as a commutative identity of ordered tournaments.

## 4. Sharp cost within the endpoint-reversal gadget bank

For all positive integers a,b, let g=gcd(a,b). Then

    gcd(2^a+1,2^b+1) = 2^g+1 if v_2(a)=v_2(b),
                       1     otherwise.                       (4)

**Proof.** Write a=gA,b=gB with gcd(A,B)=1 and x=2^g. For a common divisor d of x^A+1 and x^B+1, d is odd and x is invertible modulo d. If A,B have different parity, x^(AB) is simultaneously1 and−1 modulo d, so d=1. If A,B are both odd, choose rA+sB=1. The sum r+s is odd, so x≡(−1)^(r+s)=−1 modulo d. Thus d divides x+1; conversely x+1 divides both odd powers plus1. The two odd reduced exponents are equivalent to equal original 2-adic valuations. ∎

**Theorem C.** For every integer k≥1, a bank of k endpoint-reversal gadgets E_(a_i+2) has pairwise-coprime H-values iff the v_2(a_i) are all distinct. Every such bank satisfies

    sum_i |E_(a_i+2)| ≥ 2^k−1+2k,
    max_i |E_(a_i+2)| ≥ 2^(k−1)+2.                             (5)

Equality in the total-order bound holds exactly when the unordered exponent set is {1,2,4,...,2^(k−1)}. Thus the first k Fermat gadgets are the unique minimum-total-order coprime bank within this family, up to permutation.

**Proof.** Sort the distinct valuations as t_0<...<t_(k−1). Then t_j≥j and a_j≥2^(t_j)≥2^j. Summing gives (5). Equality in the sum forces equality in every term. The maximum bound follows from the last term. ∎

This is a **family-scoped** exponential cost barrier, not a lower bound for arbitrary tournament encodings. Gcd/lcm-faithful endpoint-reversal banks must pay it; divisibility-only banks need private primes and need not be pairwise coprime. No global minimality or polynomial-size impossibility is claimed.

## 5. Lost information and failed transfers

- **Numerical primality is not preserved.** The input prime13 maps to F_5=4294967297=641·6700417. Primes correspond to irreducibles of the *image submonoid*, not primes in all positive integers. This distinction prevents an unwarranted Chan–Pak prime-lane conclusion.
- **Addition is not preserved.** Phi(2)+Phi(3)=3+5=8, while Phi(5)=17. In fact any map with values all odd cannot preserve ordinary addition, since the sum of two odd image values is even. An operation on the image can be defined by decode/add/re-encode, but that supplies no local tournament or graft mechanism.
- **Induced containment is not divisibility.** E_3 occurs as an induced subtournament of E_4 (vertices0,1,3), yet3 does not divide5. More generally E_(a+2) embeds in E_(b+2) for a≤b by choosing endpoints and a interior vertices. The scalar lattice is not the raw induced-subtournament order.
- **Order comparison is not the target relation.** For example Phi(5)=17>Phi(6)=15 although5<6. No preservation of the usual linear order is asserted.
- **Injectivity does not imply all H-values are used.** For example11 is not in the Fermat-generated image: every nonunit image has a factor F_i, and no F_i divides11. The image is deliberately sparse and contains no argument that every allowed odd target is realized.
- **The chosen gadgets are explicit, not an intrinsic tournamentization of the divisibility relation.** Vertices are gadget positions, pairwise orientations are increasing order plus the specified endpoint reversal and SCC order. There are no ties. The faithful arithmetic predicate is the divisibility of actual Hamilton-path counts. The fixed prime-index dictionary carries the semantic identification of the source primes; it is not recoverable without specifying the construction.
- **Scalar profiles and graft continuations.** The embedding neither determines endpoint-response profiles from H nor repairs the inherited reroot/profile collision of THM-5092. It creates a faithful multiplication carrier, not a successor operation on arbitrary tournament states.

## 6. Exact controls and reproduction

Run from the repository root:

    python3 05-knowledge/results/session241_prime_relations_20261110/constructions/verify.py

Output is `exact.json` in this directory. The verifier is self-contained, standard library only. It computes H by endpoint subset DP on the actual adjacency lists and independently checks near-transitive gadgets by permutation enumeration for a≤6. It tests (4) for all4096 pairs1≤a,b≤64; all arithmetic claims for the961 pairs1≤m,n≤31; actual T(n) for1≤n≤16 whenever the order is≤15; and the exact minimum over every k-subset of{1,...,16}, k≤5, with distinct2-adic valuations. It prints the actual consequence values, the false-divisibility witness, the induced containment witness, and the composite-prime image.

These are **FINITE-EXACT controls**. The unbounded statements rest on the proofs above. The finite injectivity box for (15,51,85) has343 vectors with coordinates0,...,6; its all-order injectivity rests on determinant−2, not that box. No theorem identifiers are reserved here.
