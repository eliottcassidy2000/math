# Two uniform probes recover aggregate start skew

## 1. Universal identity under profile-preserving quotient converse

Let Q be any tournament on m quotient vertices and let

    T=Q[A_1,...,A_m],    U=Q^op[B_1,...,B_m].

Each A_i and B_i uses the same nonempty set of physical names. Assume their complete ordered path-cover profiles agree: q_k(A_i)=q_k(B_i) for every k. Their interiors need not agree. Let I be any union of whole child sets. Write S_I(T) for the number of Hamilton paths in T whose first physical vertex lies in I. Use fixed ports x->y, with a_u=[u->x] and c_u=[u->y], and Delta(a,c)=F_T(a,c)-F_U(a,c).

Then the following exact identity uses two actual uniform physical queries:

    2[S_I(T)-S_I(U)] = Delta(1_I,0)-2 Delta(0,1_I).       (1)

The masks 1_I set the port incidence of every vertex of I, not a single selected vertex. No sum over singleton queries is needed. Formula (1) applies to arbitrary quotient order and arbitrary child orders and profiles. It includes ordinary T versus T^op by using singleton children. It is an identity under the specified comparison, not a formula for the absolute start count of one arbitrary tournament from those two queries alone.

### Complete physical-word proof

Associate a variable z_i with each child and let K_ij be the formal generating series of all nonempty quotient walks starting at i and ending at j, weighted by the product of variables of their successive vertices. Consecutive vertices differ and follow a quotient edge. Reversing the quotient transposes K by reversing the walks.

Let h be the sum of all entries of K, r the sum of entries whose starting child lies in I, and c the sum of entries whose ending child lies in I. Transposition fixes h and exchanges r,c. Apply the common linear functional

    L(product z_i^k_i)=product q_k_i(A_i),   q_0(A_i)=0.

Every profile has finite support. A Hamilton path decomposes uniquely into maximal nonempty child runs, filled by ordered covers in their order of occurrence. Therefore L(h)=H(T)=H(U), L(r)=S_I(T), and L(c)=S_I(U).

Include both ports in all possible positions, with their actual edge x->y. For query (1_I,0), the two possible port orders give the full generating series

    W_(I,0)=(1+c)(1+h)+c(1+h-r).

For x before y, the middle old segment must be empty because no old vertex points to y. For y before x, the middle segment must be nonempty, end in I, and the suffix starts outside I. These are exactly the two displayed terms.

For query (0,1_I), x is a universal source relative to the old names. Only the x-before-y term is possible, and

    W_(0,I)=(1+c)(1+h-r).

Taking the difference under transposition gives polynomial identities in the commuting formal series h,r,c:

    W_(I,0)(K)-W_(I,0)(K^t)=2(c-r)(1+h),
    W_(0,I)(K)-W_(0,I)(K^t)=(c-r)(2+h).

Their first expression minus twice the second is 2(r-c). Applying L proves (1). Products under L are not replaced by products of numerical evaluations; the cancellation is performed before L is applied. Thus repeated child runs and all physical cubic and quartic contributions are retained.

This refines the singleton start-count identity in `covered_quotient/PROOF.md`, Section 1. The physical-word method follows `covered_response/PROOF.md` and `coarse_mod8_boundary/PROOF.md`; the new point is the exact cancellation for an entire union of children under profile-preserving quotient converse.

## 2. Two unequal transitive children and an arbitrary third child

**All-order separator.** Let a,b be distinct positive integers, let C be any tournament with at least two vertices, and compare

    T=C3[TT_a,TT_b,C],    U=C3^op[TT_a,TT_b,C]

on the same physical names with interiors unchanged. The two queries

    (a_mask,c_mask)=(1_C,0), (0,1_C)

cannot both have equal responses. The same conclusion holds if the compared children are replaced by any tournaments on the same physical child sets with the same respective complete profiles.

For a nonempty tournament D put p_D(k)=q_k(D), extended by zero outside its support. The cyclic-run identity from `covered_quotient/PROOF.md`, Section 2, gives

    S_C(T)-S_C(U)
      = sum_(k>=1) p_C(k)[p_a(k)p_b(k-1)-p_a(k-1)p_b(k)].       (2)

Indeed, a quotient run word beginning in C cycles through C,A,B. Of its three possible terminal positions, only the run-count pattern (k,k,k-1), in that order, changes under converse. Filling its runs by ordered covers proves (2).

For a transitive tournament,

    p_m(k)=k! S(m,k),

where S(m,k) is a Stirling number of the second kind. If a>b>=1, all minors in brackets in (2) are nonnegative and their k=2 value is

    2^a-2^b > 0.                                    (3)

For completeness, consecutive Stirling rows have increasing ratios. Each row is log-concave with consecutive positive support, as follows by induction from s_(m+1,k)=k s_(m,k)+s_(m,k-1): if b_k=k a_k+a_(k-1), then

    b_k^2-b_(k-1)b_(k+1)
      = k^2(a_k^2-a_(k-1)a_(k+1)) + a_(k-1)a_(k+1)
        +(k+1)(a_k a_(k-1)-a_(k-2)a_(k+1))
        +(a_(k-1)^2-a_(k-2)a_k) >= 0.

The wider minor is nonnegative by successive log-concave ratios, and zero-extension handles boundaries. Thus

    s_(m+1,k)/s_(m,k)=k+s_(m,k-1)/s_(m,k)

is strictly increasing wherever defined. Multiplying positive consecutive-order ratios proves the required two-row minor inequality on common support; at and beyond the smaller row's upper endpoint it is immediate. Factorials do not change its sign. Equation (3) also covers the smaller order one.

Every tournament of order at least two has q_2>0: partition it into two nonempty parts, choose a Hamilton path in each part, and order those paths. Consequently all terms in (2) have one weak sign and the second term is strictly signed. In fact

    |S_C(T)-S_C(U)| >= q_2(C)|2^a-2^b| > 0.

Formula (1) now proves separation by the two advertised physical queries. No log-concavity assumption on the arbitrary third child's profile is used. This strictly extends the earlier unequal-three-transitive-child result by allowing the third child to be arbitrary and reducing its singleton battery to two uniform probes.

If a=b, this start-skew certificate is identically zero; the equal-profile and repeated-TT2/TT3 results remain separate methods. More strongly, C3[TT2,TT2,C3] has all three aggregate start skews zero despite unequal child profiles: the only nontrivial minor is 2(2*3-6)=0. Nevertheless the physical query (3,12) has values 1959 and 2043. Start-skew blindness is therefore not coarse or full response blindness. If C is a singleton, the strict second summand is absent and this argument alone gives no conclusion for this choice of queried child.

## 3. Original-module root-family corollary

Assume the same unequal transitive pair and |C|>=2. Select any family of modules of the original T, allowing overlaps or nesting within children, with the whole root selected an odd number of times. Reverse every selected interior by edge XOR, giving V.

The fixed three-query battery

    (0,0), (1_C,0), (0,1_C)

distinguishes T from V. This is original-versus-every-root-containing-family separation, not an all-pairs orbit identification claim.

A proper module of the cyclic blowup lies inside one child: a set meeting two children is distinguished by the third; meeting all three forces every child to be included. Thus V has converse root quotient. The independently established original-module monotonicity and hereditary equality theorem in `module_decomposition/PROOF.md` says H(V)>=H(T), and equality forces all corresponding child profiles to agree. If strict, F(0,0)=H separates. If equal, apply Sections 1–2 to the profile-preserving root converse. No restriction that the selected family be laminar is used.

This is the same root-family reduction used in `covered_root_family/PROOF.md` and `arbitrary_third_root_family/PROOF.md`, now combined with the two-uniform-query identity and the arbitrary-third unequal-transitive-pair theorem. The whole-core scope is essential here. Nonzero integer gaps are not claimed to survive arbitrary strong outer ambient substitution.

## 4. Exact reproduction and controls

The self-contained standard-library checker verifies the two polynomial identities coefficient by coefficient. It checks all labeled tournaments through order four against their ordinary converses for every subset of physical names, then checks unions of child modules in three nontrivial substitution examples, including a prime-five quotient and changed child interiors with identical profiles. All checks use literal augmented-tournament Hamilton-path dynamic programming and direct physical start counts.

It also checks the Stirling minors through order sixteen and arbitrary-third examples with transitive, cyclic, and strongly connected four-name third children, including reversed size order and a singleton member of the unequal transitive pair. The complete two-query values and start counts are retained. Sixteen original root-containing families on C3[TT2,TT3,C3], including crossing TT3 intervals, check the three-query corollary. The unequal-profile all-three-start-skew obstruction above is also checked literally. These finite controls validate the implementation; the all-order statements follow from the proofs above.

    python3 verify.py > normal.json
    python3 -O verify.py > optimized.json
    cmp normal.json optimized.json

All mathematical comparisons use explicit errors and remain enabled under optimization. The four-file source packet is `PROOF.md`, `verify.py`, `normal.json`, `optimized.json`.
