# Fixed cyclic caps, arbitrary seed eigenmodes, and physical-context tomography

Status: elementary deductions PROVED by the arguments below, independently audited; see ../audit/FIXED_CAPS_AUDIT.md; computations FINITE-EXACT. No Chan–Pak completeness claim. All tournaments and substituted blocks below are finite and nonempty. H counts directed Hamiltonian vertex words.

## Inheritance and scope

The source snapshot is `0f8750677dd20b594a0e032a3b56aa23f58a950c`.

- THM-1975, `01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md`: scalar H does not compose; path-cover data do.
- THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`: audited exact quotient-walk substitution and balanced directed-triangle kernel.
- THM-3134, `01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`: ordered covers as permutation-polynomial jets and Newton transforms.
- THM-3166, `01-canon/theorems/THM-3166-falling-factorial-order-join-path-colour-transform.md`: Q is multiplicative under ordinal joins; earlier recovery via evaluations is inherited.
- THM-3181, `01-canon/theorems/THM-3181-tournament-half-grid-reciprocity-and-repeated-join-recurrence.md`: order parity and half-grid recovery.
- THM-3202, `01-canon/theorems/THM-3202-c3-repeated-join-moving-jet-formula-and-cfinite-obstruction.md`: all three growing factors require a moving jet and have a factorial-cube lower bound.

The contribution is an alternative **symmetric transitive cyclic-cap probe family**, with diagonal (m!)², and the explicit fixed-cap mixed-seed compiler/non-dyadic pump. **Actual positive-context profile recovery and contextual equivalence are already proved in Session176**, `05-knowledge/results/session176_tournament_information_20261104/information/PROOF.md`, Sections 2–3. Its probes P_r use a transitive exterior a_1,b_1,...,a_r,b_r with T->a_i and b_i->T and diagonal r!. The known-order parity compression and exact realized affine rank are already in `05-knowledge/results/proof_principles_20261103/profiles/PROOF.md`, R1–R4. These sparse-tree sources were recovered by the independent auditor via pinned `git show`; overlooking them initially would have overstated novelty. Section 3 below specializes their observational principle to a different symmetric family, and does not claim generic physical tomography or its compressed probe count as new.

The universal kernel, polynomial transforms, parity, ordinal products, finite recurrences at fixed depth, and all generic recovery/equivalence principles are inherited. THM-3166 already recovers profiles from negative Q evaluations; Session176 already upgrades recovery to actual positive tournament probes. The specific cap coefficient formula and arbitrary-seed specialization are offered as explicit compositional packaging and applications, not a new universal theory or literature-priority claim.

Concept board: scalar H / bounded path-cover jets / ordinal products / cyclic interleaving / permutation descents / actual physical contexts. Anchor: construct exact odd H values. Niche: arbitrary seed pumping instead of singleton pumping. Wildcard: recover a complete compositional invariant using only H of actual tournaments. The hostile is confusing fixed caps with balanced growing caps, or confusing an algebraic evaluation with a physical tournament.

## 1. Exact cap compiler

Write F_T(k)=k! pc_T(k), where pc counts unordered spanning directed path covers, and set F_T(k)=0 for k outside 1..|T|. Let

    Q_T(t)=sum_k binom(t,k) F_T(k).

THM-3121 gives, for the cyclic substitution A -> B -> C -> A,

    H(C3[A,B,C]) = sum_{x,y,z>=1} w(x,y,z) F_A(x) F_B(y) F_C(z),

where w=3 if x=y=z, w=1 if max(x,y,z)-min(x,y,z)=1, and w=0 otherwise. This identity is inherited; its direct bijection splits a global path into maximal monochromatic block runs and assigns ordered local covers to the occurrences of each block.

Fix A,C and put m=min(|A|,|C|). Define the nonnegative integers

    c_k(A,C)=sum_{x,z>=1} w(x,k,z) F_A(x) F_C(z),  1<=k<=m+1.

Then

    H(C3[A,B,C])=sum_{k=1}^{m+1} c_k(A,C) F_B(k).              (1)

Indeed, nonzero w requires k<=min(x,z)+1<=m+1. Thus fixed caps see a bounded jet regardless of B's size. The cap sizes need not be equal and the caps need not be transitive.

Define possibly signed integers

    beta_j(A,C)=sum_{k=j}^{m+1} (-1)^(k-j) binom(k,j) c_k(A,C).

Newton inversion F_B(k)=sum_{j=0}^k (-1)^(k-j)binom(k,j)Q_B(j), with Q_B(0)=0, proves

    H(C3[A,B,C])=sum_{j=1}^{m+1} beta_j(A,C) Q_B(j).          (2)

The sign of beta has no positivity interpretation: (1), not (2), is the positive cover count.

## 2. Arbitrary seed-word composition and a real non-dyadic pump

Let B=D_1 => ... => D_r be an ordinal sum of arbitrary nonempty tournaments (all arcs from earlier to later blocks), with r>=1. The inherited ordinal product identity yields

    H(C3[A,B,C])=sum_{j=1}^{m+1} beta_j(A,C) product_i Q_{D_i}(j).  (3)

This is a bounded-register compiler for arbitrarily long words over any finite seed bank. The j-th register multiplies by Q_D(j) on reading seed D; the final linear functional has coefficients beta_j. Hence the register multiplications commute. Order of seed blocks can change the tournament but not this profile or any cyclic substitution built from it. This is an intentional loss, not a claim that tournaments become isomorphic.

For a repeated seed B=D^{=>r}, (3) is an exact exponential sum in r. It satisfies the constant-coefficient recurrence with annihilator

    product_{j: beta_j != 0} (E-Q_D(j)),

where E shifts r by one. Repeated equal eigenvalues may be merged; the displayed operator is valid, not asserted minimal. The statement concerns r>=1. An algebraic r=0 extension is not an empty-block C3 construction.

For C a singleton, direct collection in (1) gives

    H(C3[A,B,1])
      =(Q_A(2)-H(A))Q_B(2)+(3H(A)-Q_A(2))H(B).             (4)

In particular

    H(C3[1,D^{=>r},1])=Q_D(2)^r+H(D)^r.                   (5)

A concrete five-vertex seed D is given by adjacency rows

    00011
    10000
    11000
    01100
    01110

Its ordered-cover profile is (11,78,198,240,120), so H(D)=11 and Q_D(2)=100. Consequently, for every r>=1 there is a **strong tournament of order 5r+2** with

    H=100^r+11^r.

The first values are 111,10121,1001331,100014641. The recurrence is u_(r+2)=111u_(r+1)-1100u_r. Strongness holds for every C3 substitution of nonempty blocks: each block reaches either other block around the directed triangle, and two vertices inside a block can be connected by a detour through the other two blocks. This genuinely replaces singleton eigenmodes j by arbitrary physical seed eigenmodes Q_D(j); it is not a claim that this sparse family fills an odd tail.

Connection to the concurrent multi-exit compiler: setting every D_i to a singleton makes Q_D(j)=j and yields its integer exponential modes. For A=C3 and C=TT_3 the independently computed vector is

    beta=(87,-189,72,36),

giving H(C3[C3,TT_r,TT_3])=87-189*2^r+72*3^r+36*4^r for r>=1. This is the same family as the companion multi-exit descent-polynomial EGF, not an independent new family. The descent distribution and F are related by the inherited permutation-cut identity F_T(k)=sum_b A_T(b) binom(|T|-1-b,k-1-b), with invalid binomials zero. Formula (3) extends these modes to arbitrary tournament seeds.

## 3. Actual positive contexts separate exactly the path-cover profile

Let TT_m be the transitive m-vertex tournament and define

    Y_m(T)=C3[TT_m,T,TT_m],  m>=1.

Its two caps have F_(TT_m)(m)=m!, so equation (1) is triangular:

    H(Y_m(T))=sum_{k=1}^{m+1} c_mk F_T(k),
    c_m,m+1=(m!)^2>0.                                   (6)

An optional closed coefficient formula, with f_k=F_(TT_m)(k)=k! S(m,k) and f outside 1..m set to zero, is c_mk=3f_k^2+2f_k(f_(k-1)+f_(k+1))+f_(k-1)^2+f_(k+1)^2. It follows by listing the seven allowed (x,z) pairs around k.

To prove the leading coefficient, k=m+1 forces x=z=m. That kernel cell has weight one. No other cell contributes. Hence H(T), followed by H(Y_1(T)),...,H(Y_(n-1)(T)), reconstructs every F_T(k) for an n-vertex T by

    F_T(m+1)=[H(Y_m(T))-sum_{k=1}^m c_mk F_T(k)]/(m!)^2. (7)

Each observation is H of an actual strong tournament, with no negative graph, signed tournament, or virtual vertex. This family shares the physical-observer advantage already established by Session176; its particular symmetric caps and diagonal differ from Session176’s transitive-exterior probes.

**Context equivalence.** For arbitrary nonempty tournaments S,T (even of unequal order), the following are equivalent:

(a) their complete F profiles agree after zero extension;
(b) H(S)=H(T), and H(Y_m(S))=H(Y_m(T)) for every m>=1.

The forward direction is (6); the reverse is induction in (7). Since F_T(|T|)=|T|! and higher entries vanish, full profile agreement also forces equal order. Thus this small explicit class of positive contexts is a complete observer for exactly the invariant needed by all substitution kernels, although not a complete observer for tournament isomorphism.

**Known-order compression.** If n is known, put d=floor((n-1)/2). H(T) and the first max(0,d-1) cap observations suffice. For n=1,2 the profile is forced. Otherwise equations (6)-(7) give F_T(1..d), hence Q_T(1..d). The inherited parity Q_T(-t)=(-1)^n Q_T(t), monicity of degree n, and Q_T(0)=0 leave exactly d unknown coefficients: powers 1,3,...,n-2 for odd n, and 2,4,...,n-2 for even n. Their values at distinct positive j=1..d determine them by a Vandermonde matrix in j^2, after dividing rows by j (odd n) or j^2 (even n). The determinant is nonzero. This is an interpolation upper bound, not a claim of minimality on the realizable tournament set. THM-3181 gives the parity mechanism, and proof_principles/profiles, R3 already proves compressed physical-probe recovery. Thus the compression itself is inherited, applied here to the symmetric cap family.

## 4. Hostiles, boundaries, and what remains open

- Stronger same-order hostile (inherited Session176 Section 5, relabelled): five-vertex masks 69 and 76 both have H=15, but F=(15,78,198,240,120) and (15,90,210,240,120). Singleton cyclic caps give H=123 and H=135. The verifier checks these directly.
- Scalar-H replacement fails: TT_1 and TT_2 both have H=1 but singleton cyclic caps give H=3 and H=5. Their profiles and orders differ.
- Positive cover coefficients do not imply positive spectral coefficients: caps TT_2 and TT_1 give beta=(-1,3), so (2) is 3Q_B(2)-H(B).
- The coefficient in (4) is **3H(A)-Q_A(2)**, not 4H(A)-Q_A(2); TT_2 caps catch this off-by-one error immediately.
- Empty caps are excluded. Removing a block breaks the strong cyclic quotient and the formula's positive-content premise.
- Fixed m bounds the needed jet. If both caps grow, m grows; the THM-3202 factorial-cube obstruction and moving modes are not contradicted.
- Strongness of outputs and infinitely many explicit odd counts are proved. Cofinite coverage, target inversion for all odd N>=23, and any improvement over published finite coverage remain OPEN. The new spectrum family may consist entirely of values achievable by other constructions; no numerical first-occurrence claim is made.
- The map is maximal-block-run splitting. It preserves full cover counts; scalarization discards run multiplicities. The restoring sidecar is F or equivalently Q, and the physical cap tests prove that sidecar cannot be discarded if every such composition's H is to be preserved.

## 5. Deterministic evidence

Run `python3 verify.py` and `python3 -O verify.py` from this directory. Outputs are byte-identical to `receipt.json` and `receipt.optimized.json`.

The explicit universe is all 75 labelled tournaments of orders 1..4, plus named deterministic larger cases. A subset-Hamiltonian/set-partition engine computes F; independent permutation-and-cut enumeration validates all 75 profiles and the five-vertex seed. A separate forward Held–Karp engine counts H of 160 actual cyclic compositions, including independent physical cap tomography through 13 vertices, two-copy seed pumping, and asymmetric mixed-word compositions. Profile reconstruction has 76 tests, and parity-compressed reconstruction has 79 tests, including named orders 6,7,8. No assertions are used: checks remain active under Python -O. The 87,-189,72,36 bridge is checked exactly. These are finite controls of general proofs, not a bounded search substituted for a proof.
