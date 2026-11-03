# Slot representations turn moving guards into Young-lattice dynamics

**Status.** Exact representations and constructive transfer rules, followed by an exponential asymptotic theorem for a literal strong-tournament family. The companion [uniform theorem](UNIFORM_PLANCHEREL.md) proves quadratic Plancherel convergence for all growing tournaments, with optimal order. These are refinements of the full-cover calculus, not a solution of the odd Hamilton-count spectrum. No claim of priority against the complete literature is made.

The input is [S164's pointed profile algebra](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session164_emergent_20261103/profile_algebras/POINTED_PROFILE_ALGEBRA.md), with its ordinary [matching product](https://github.com/eliottcassidy2000/math/blob/89416555c4d95c6fa20c9ebafca1a180ad3fcaf3/05-knowledge/results/session164_emergent_20261103/profile_algebras/PROFILE_ALGEBRA.md). The odd-cycle expansion used in Section 5 belongs to established Rédei–Berge theory: [Grinberg–Stanley, arXiv:2307.05569](https://arxiv.org/abs/2307.05569), with a matrix-algebra treatment by [Irving–Omar](https://arxiv.org/abs/2412.10572). The constructive contribution here is to combine the **slot** action with physical guard moves and obtain positive Young-lattice transfer and a spectral limit. Slot representations must not be confused with a Schur-positive vertex symmetric function.

All representations below are finite-dimensional over the complex numbers. Direct sums and tensor products are actual representations, not merely scalar count operations.

## 1. The correct symmetric group acts on slots

Let X be a tournament. Let C_X(r) be its set of ordered r-slot path covers, allowing empty slots, and let V_X(r)=C[C_X(r)]. The group S_r permutes slots. It does not relabel vertices or claim to act by automorphisms of X. Write p_k(X) for the number of unordered k-path covers, and

    F_X(t) = sum_k p_k(X) (t)_k.

Let R_r be the permutation representation on the r slots. For every sigma in S_r,

    character(V_X(r))(sigma) = F_X(fix(sigma)).                 (1)

Indeed, a nonempty path cannot occupy two distinct slots of the same cover. Thus a fixed cover uses only individually fixed slots. In particular, every derangement has trace zero when X is nonempty.

For each unordered k-path cover, the injections of its k distinct paths into r slots form one transitive S_r-set, with stabilizer S_(r-k). Consequently

    V_X(r) ≅ direct sum_(0<=k<=r) p_k(X) Ind_(S_(r-k))^(S_r) 1.  (2)

The omitted terms k>r contribute zero. Under the Frobenius characteristic map,

    ch V_X(r) = sum_(0<=k<=r) p_k(X) h_(r-k) h_1^k.             (3)

If lambda partitions r, its Specht multiplicity is therefore

    m_lambda(X;r) = sum_(0<=k<=r) p_k(X) f^(lambda/(r-k)),       (4)

where f^(lambda/(r-k)) counts standard tableaux of the skew diagram obtained by removing a first row of length r-k; it is zero if that row does not fit. Formula (4) follows directly by repeated one-box branching from the trivial S_(r-k) representation. For r-k=0, the inner diagram is empty.

This already gives a precise representation-theoretic reading of every cover coefficient. For example the sign multiplicity is

    multiplicity(sign, V_X(r)) = p_(r-1)(X)+p_r(X).             (5)

There is one caveat: at a fixed rank, these modules need not distinguish all profiles. The injection modules at k=r-1 and k=r are both regular representations. We do not claim that a single small S_r-module determines an arbitrary unbounded profile.

### Composition becomes tensor product

Every slot path in X⊕Y splits uniquely at the block boundary. Slot permutations commute with this splitting, so there is an actual equivariant bijection

    C_(X⊕Y)(r) ≅ C_X(r) × C_Y(r),

and hence

    V_(X⊕Y)(r) ≅ V_X(r) ⊗ V_Y(r).                             (6)

At the injection-module level the decomposition is

    M_(r,k) ⊗ M_(r,l)
      ≅ direct sum_t binom(k,t) binom(l,t) t! M_(r,k+l-t),      (7)

with M_(r,j)=0 when j>r. Here t is the number of slots shared by the two injections. Choosing which k- and l-paths share slots and matching them proves (7). Thus the matching product is exactly a positive tensor-product rule, not a resemblance to one. The tensor symmetry does not imply that X⊕Y and Y⊕X are isomorphic tournaments.

## 2. Positive pointed functors and the triangle reset

Suppose q has unique in-neighbor rho. Refine V_X(r) into the representations S,U where q respectively starts or does not start its path. Let d be the slot-cover representation of X−q restricted to covers where rho starts its path. These are invariant subsets, so all three representations are genuine permutation modules.

Let W_r be the augmentation-zero subspace of R_r. For r>=2 it is the standard representation S^(r-1,1), and R_r=1⊕W_r. For r=1 set W_1=0.

Add a vertex a with q→a and a→X−q. Mark a for move C, or retain q for move D. S164 proves the corresponding scalar insertion bijections. Keeping their slot action gives the following stronger identities:

    C(S,U,d) = (W_r⊗S ⊕ R_r⊗U, S⊕U⊕d, S),

    D(S,U,d) = (R_r⊗S ⊕ d, (R_r⊕1)⊗U, W_r⊗d).             (8)

Every operation in (8) is positive in the representation semiring.

Here is why the factors W_r are legitimate, rather than formal subtraction. In a cover in S, let i be the slot beginning at q. The allowed initial insertion slots are j≠i. Send the basis element (j,cover) to

    (e_j−e_i) ⊗ cover.

For each fixed cover these r−1 vectors form a basis of W_r; the map is equivariant under simultaneous slot permutations. The same argument applies to d with the marked slot beginning at rho. The remaining insertion classes in S164 are equivariant without any basis change. For example the r choices counted in the first coordinate of D are the r−1 other initial slots and the insertion immediately after q, identified with the marked slot i.

Ordinal grafting on the right becomes the componentwise module action

    (S,U,d)_(X⊕Y,q) ≅ (S⊗V_Y(r), U⊗V_Y(r), d⊗V_Y(r)).       (9)

The graft is not generally strong. By contrast both C and D preserve strength of a strong input: rho→q→a→rho attaches the new vertex to the existing strong component.

### A representation-level triangle lift

Let G(X)=X→a→b→X and mark b. Its three modules satisfy

    d_G(r) ≅ Ind_(S_(r-1))^(S_r) V_X(r-1),
    U_G(r) ≅ R_r ⊗ Res_(S_r)^(S_(r+1)) V_X(r+1).              (10)

For d, a is a singleton in a distinguished slot and X occupies the others. For U, the segment a,b lies in some distinguished ordinary slot i. Split that slot's X-vertices before a and after b into slot i and an extra, distinguished slot r+1. This is an equivariant bijection proving the second identity.

There is also a positive concrete complement description of S_G:

    R_r⊗R_r⊗V_X(r) ≅ S_G(r) ⊕ d_G(r).                       (11)

Given an X-cover, append a in slot i and prepend b in slot j. Every pair is valid except i=j with that old slot empty. The excluded pairs are exactly the distinguished-empty-slot model for d_G. Thus (11) is a split inclusion of actual permutation modules, not a subtraction-based existence claim.

Equations (6), (8)–(11) give one constructive category of slot modules covering ordinal composition, strong triangle resets, physical moving guards and forgetting the guard. They do not make arbitrary representation triples realizable as tournaments.

## 3. A nonnegative Young-lattice engine

Fix r>=1. Index vectors by partitions lambda of r. Define

    K_(lambda,mu) = number of partitions nu of r−1
                   obtainable by removing one corner from both lambda and mu.

The standard tensor identity R_r⊗V ≅ Ind Res V and one-box branching imply that K sends the Specht multiplicity vector of V to that of R_r⊗V. Moreover K−I is entrywise nonnegative: off-diagonal entries already are, and K_(lambda,lambda) is the positive number of removable corners of lambda. It sends multiplicities to those of W_r⊗V.

Therefore the literal pointed moves have the **nonnegative integer** transfer rules

    C(s,u,d) = ((K−I)s+Ku, s+u+d, s),
    D(s,u,d) = (Ks+d, (K+I)u, (K−I)d).                       (12)

This is a finite multitype branching system on three copies of the partitions of r. A K-step removes one box and adds one box. It computes every Specht multiplicity of every physically constructed C/D word without enumerating path covers of the growing tournament. The state retains the endpoint information needed to compose the operations.

The same semiring calculation yields a representation-valued recurrence. If V_j is the full cover module after j consecutive C-moves, then

    V_(j+3) ≅ R_r⊗V_(j+2) ⊕ V_(j+1) ⊕ R_r⊗V_j.            (13)

To verify it before any specialization, put R=1+W in the positive 3×3 block matrix in (8). Its cube equals R times its square plus itself plus R times the identity, entry by entry as polynomials in W. Thus (13) is an isomorphism of genuine representations, and its multiplicity vectors obey

    m_(j+3)=K m_(j+2)+m_(j+1)+K m_j.                         (14)

### Explicit strong example

Start with a directed edge and mark its sink. Its full module is R_r^2, and one C-move gives a directed triangle. The initial full modules are

    V_0 ≅ R_r^2,
    V_1 ≅ R_r^3 ⊕ 2R_r,
    V_2 ≅ R_r^4 ⊕ 4R_r^2.

All j>=1 are strong. If e is the vector of the trivial Specht representation, their multiplicity generating function is

    sum_(j>=0) m_j z^j
      = [I−zK−z²I−z³K]^(-1) [(1+z²)K²+2zK] e.             (15)

The inverse is formal, so (15) is an exact finite transfer formula. The tournaments here are S164's C-chain, including its identified S159 subsequence; this is not a new graph family.

For r=3, in the order trivial, standard, sign, the first strong examples are

| Vertices | Hamilton paths | Full 3-slot dimension | Specht multiplicities |
|---:|---:|---:|---|
| 3 | 3 | 33 | (7,11,4) |
| 4 | 5 | 117 | (22,39,17) |
| 5 | 9 | 411 | (73,137,64) |
| 6 | 17 | 1449 | (250,483,233) |
| 7 | 31 | 5109 | (867,1703,836) |

The standard summand has dimension two. Indeed the exact S3 formula for every nonempty X is

    m_triv = (F_X(3)+3H(X))/6,
    m_standard = F_X(3)/3,
    m_sign = (F_X(3)−3H(X))/6.                               (16)

### A six-coordinate exact observer at rank three

There is an especially small positive quotient that keeps the Hamilton count visible. Let Q be the regular representation of S3 and R=R_3. For any nonempty tournament,

    V_X(3) ≅ H(X) R ⊕ (p_2(X)+p_3(X)) Q.                    (17)

This follows directly from (2): the one-path injection module is R, and the two- and three-path injection modules are both Q. In particular

    H(X)=m_triv(X;3)−m_sign(X;3).

The Hamilton count is the exact imbalance between the two one-dimensional sectors. The full module is never exactly regular, even though its normalized character can approach the regular one.

Each of S,U,d has the same two-generator form. Write its coordinates as (h,b), where h is its one-slot count and b its total unordered two- and three-path count. The identities

    R²=R⊕Q,   R⊗Q=3Q,   W⊗R=Q,   W⊗Q=2Q

turn (8) into six **nonnegative** integer updates. For input ((h_s,b_s),(h_u,b_u),(h_d,b_d)), C gives

    S: (h_u, h_s+2b_s+h_u+3b_u),
    U: (h_s+h_u+h_d, b_s+b_u+b_d),
    d: (h_s,b_s),

while D gives

    S: (h_s+h_d, h_s+3b_s+b_d),
    U: (2h_u, h_u+4b_u),
    d: (0, h_d+2b_d).                                        (18)

Thus nine coefficients through cover-degree three compress to six integers without losing any S3-module information. A right ordinal graft also closes here: tensoring coordinates (h,b) and (h',b') gives

    (hh', hh'+3hb'+3h'b+6bb').                              (19)

This observer does not determine the entire unbounded profile or support every triangle reset at rank three. It is an exact smaller positive observer for the stated moving-guard and right-graft language.

## 4. Literal guard growth forces a Plancherel limit

**Theorem.** Start with any valid guarded tournament and apply C repeatedly. Fix r>=3. Let m_lambda(j;r) be the multiplicity of S^lambda in the r-slot cover representation, and let D_j=F_(X_j)(r) be its dimension. If d_lambda=dim S^lambda, then

    m_lambda(j;r)/D_j
      = d_lambda/r! + O((rho(r−2)/rho(r))^j),                (20)

where rho(t) is the positive dominant root of

    z³−t z²−z−t=0,     t>=1.

Thus the fraction of the entire representation in the lambda-isotypic component tends to

    d_lambda²/r!,                                          (21)

the Plancherel weight. The constants can depend on the fixed seed and r. This is a limit in a growing family of finite permutation modules, not a statement that the finite modules are exactly regular.

**Proof.** At a permutation with t fixed slots, the character is F_(X_j)(t) by (1). For every t>=1 its three pointed traces evolve by

    C_t = [[t−1,t,0],[1,1,1],[1,0,0]].

This nonnegative matrix is primitive, including at t=1: its directed graph is strongly connected and has a self-loop. The seed is nonzero, so Perron–Frobenius gives F_(X_j)(t) asymptotic to a positive constant times rho(t)^j. Comparing the irreducible nonnegative matrices at two positive t shows that rho(t) is strictly increasing. Their characteristic polynomial is the displayed cubic.

The identity has r fixed slots, whereas a nonidentity permutation has at most r−2. Derangements contribute zero. The character inner-product formula therefore gives

    m_lambda(j;r)
      = d_lambda D_j/r! + O(rho(r−2)^j).

Divide by D_j, which is asymptotic to a positive multiple of rho(r)^j. This proves (20), and multiplying by d_lambda gives (21). ∎

For r=2 the conclusion is exact at every step and for every nonempty tournament: the unique nonidentity permutation is a derangement, so the 2-slot module is a sum of copies of the regular representation. The r=1 case is tautological.

This supplies a concrete answer to what the representation viewpoint adds: the same literal guard operation that produces tribonacci Hamilton counts also drives every fixed-r irreducible sector toward a prescribed distribution, with an explicit spectral rate. It does not produce an affine ternary observer return or establish any new exact Hamilton-count coverage.

The companion [uniform Plancherel theorem](UNIFORM_PLANCHEREL.md) extends convergence to every growing tournament sequence at an optimal uniform quadratic order. The physical C-chain above has the stronger exponential rate.

## 5. The odd-cycle tensor lift, and a distinction that matters

For an odd directed cycle packing C, let c(C) be its number of cycles and

    b(C)=n−sum_(gamma in C)(|gamma|−1).

The established odd-cycle expansion, with singletons allowed outside C, specializes to

    F_X(t)=sum_C 2^c(C) t^b(C).                              (22)

One can see the specialization without introducing slot representations: apply the Hamilton odd-cycle formula separately to every ordered vertex subset assigned to a slot. A fixed global packing is assigned to slots by assigning each cycle and each unused singleton as a whole, giving t^b(C) choices.

Combining (1) with (22) gives the genuine representation isomorphism

    V_X(r) ≅ direct sum_C 2^c(C) R_r^(tensor b(C)).            (23)

Consequently

    p_k(X)=sum_C 2^c(C) StirlingSecond(b(C),k).                (24)

This connects the positive tensor category to the existing cycle-packing theory. It implies that F_X is monic of degree n, has only powers of parity n, and has nonnegative even coefficients below its leading term. These are established cycle-expansion consequences; they are not asserted as a new odd-cycle theorem. Equality of characters proves (23); we do not claim to have constructed a canonical cover-to-cycle bijection.

The **vertex-degree** Rédei–Berge symmetric function is a different object. Its odd-cycle expansion retains cycle lengths in power sums, rather than replacing every cycle/singleton by R_r. It is not automatically Schur-positive. For the cyclic regular tournament on Z/5Z with i→j when j−i is 1 or 2, direct enumeration gives five directed triangles and two directed 5-cycles. Hence

    U_X = p_1^5 + 10p_3 p_1² + 4p_5.

The S5 character of shape (3,2) on cycle types 1^5, 3·1², 5 is respectively 5,−1,0. Its Schur coefficient is therefore −5. For example these character values follow by subtracting the 1-subset permutation character from the 2-subset permutation character: binom(fix(sigma),2)+number_of_2_cycles(sigma)−fix(sigma).

There is no contradiction with (3), (8), or (23): those are genuine S_r representations of **slot permutations**, whereas U_X is a degree-n symmetric function recording **vertex cycle sizes**. A hypothetical positive representation lift of U_X requires additional hypotheses or a different construction.

## Reproduction and scope

Run `python3 check.py` and `python3 -O check.py`. Both runs produce identical JSON. The checker uses independent subset/path dynamic programming and direct path-and-cut enumeration, with explicit conditions active under `-O`.

- 363 direct fixed-cover character checks on all tournaments through order three, across all slot permutations through S4.
- 1,100 pointed Specht-transfer checks: both literal C/D moves on all 110 eligible guarded labeled tournaments through order four, at ranks 2 through 6.
- 450 independent triangle induction/restriction/complement controls.
- 220 independent checks of the six-coordinate compact transfer.
- 1,099 independent odd-cycle/Stirling profile checks on all labeled tournaments of orders 1 through 5.
- 726 ordinal character-product checks, and literal C-chain recurrence/decomposition controls through nine vertices.
- K−I nonnegativity through rank 10 and the regular T5 distinction above.

The equivariant bijections, branching identities and Perron–Frobenius argument prove the general assertions. The finite tests are controls, not extrapolated proofs. The construction remains inside the exact physical guard language, and the limit concerns its representations rather than universality of its scalar counts.
