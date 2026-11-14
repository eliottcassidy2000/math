# A factor-free reservoir inverse and commutative profile capacity

Session 324. **PROVED elementary bounds and algorithms; independent audit pending.** Chan–Pak remains OPEN. The new inverse is complete only for the specified selective-quotient grammar. Its exponent is a guaranteed upper bound, not a claim of optimal complexity or dense output.

## 1. The retained physical state

The Session 321 selective quotients have rows

    R13=(24,9,35,4,14,27),
    R15=(108,97,2,6,79,92,12).

U occupies vertex zero and V the last vertex. Their exact strong output, for any nonempty U,V, is

    N=A(U) Q_V(2)-H(U)H(V),
    A(U)=-3H(U)+6Q_U(2)+gamma Q_U(3),  gamma=2 or 6.

Here Q_T(z) counts ordered z-tuples of directed paths covering T, with empty slots allowed. Both quotients are strong, so every such inflation is strong. We use the proved layout identity from `s321/SELECTIVE_RAYS.md`, not a new fit to examples.

The useful invariant is the complete source triple (H,Q2,Q3) and reservoir pair (H,Q2). On an ordinal sum all these coordinates multiply. Their multiplication is commutative, even though the tournaments in different component orders need not be isomorphic. Thus the output forgets SCC order exactly. Counting ordered SCC words pays for information the chart does not use.

Closest mechanism: `s323/PAIRED_END_DESCENT.md` gives an exact fixed-source inverse after factoring its forced reservoir count. Hostile: the all-modulus aperture is absolutely sparse, and residue agreement leaves the absolute slot cost unpaid. Corrected near miss: an ordered-word generating function is a valid bound but needlessly requires letter mass below one. The underused coordinate is the multiset of physical SCC profiles.

## 2. An Euler-product bound without the letter-mass restriction

Let an alphabet consist of explicitly realizable strong SCC letters and the singleton. Give each letter a weight w>1. Suppose it is locally finite by weight and

    sum_letters w^(-s) < infinity.

Then the sum of the same weight over **all finite multisets** is finite:

    Z(s)=product_letters (1-w^(-s))^(-1)<infinity.

Indeed all weights in our applications are integers at least two, so the logarithm of the product is at most the convergent letter sum divided by 1-2^(-s). Multiple presentations of one SCC or profile may be retained as different letters: this only overcounts. Choose one actual representative when constructing a witness.

Consequently at most Z(s)Y^s multisets have product weight at most Y. This needs convergence of the letter sum, not that its value be less than one. A finite additional alphabet changes Z(s) by a finite factor and never changes a convergence threshold.

For the selective charts, Q2>=2H and gamma>=2 give

    N >= 2 Q_U(3) Q_V(2).                              (1)

Apply the product argument to the union of source letters weighted by Q3 and reservoir letters weighted by Q2, keeping the two roles as distinct labels. It bounds all candidate pairs, including duplicates, by O_s(X^s). More strongly, the sum N^(-s) over distinct positive outputs is finite whenever both letter series converge: select one multiset pair for each output and use (1).

This directly explains why a sufficiently small alphabet cannot cover either the ordinary-prime or 7-times-prime atom lane. Choose s<1 in its convergence range. If all but finitely many odd primes, or all but finitely many allowed 7p, were outputs, their reciprocal sum would converge, contrary to Euler's prime reciprocal divergence. This uses no assumption on distribution of primes in the explicit families.

## 3. Three quantitative consequences

For a strong one-apex-over-chain letter with chain length r>=2, the first cut bit is one and the last zero. There are at most 2^(r-2) words. The extreme inversion in the exact cycle formula gives

    Q(2)>=2*3^(r-2),     Q(3)>=3*4^(r-2).              (2)

Therefore its Q2 letter series converges for s>alpha=log(2)/log(3), and its Q3 series for s>1/2.

If the letter additionally has dyadic Q2, the S323 paired-end proof supplies all candidates in O(24^(r/8) poly(r)) arithmetic work, with only O(24^(r/8)) candidate words. This is a proved necessary pruning followed by an exact power-of-two test, not an assumption that candidates are actual solutions. By (2), the dyadic Q2 letter series converges for

    s>beta=log(24)/(8log(3)) < 1/2.                    (3)

The exceptional chain length three adds only finitely many letters. The count at each r can be used as an upper bound even if no actual dyadic letter exists there.

For every epsilon>0, the following bounds now hold for distinct selective-chart outputs through X:

1. Arbitrary ordinal sums of arbitrary strong one-apex cores on both sides, plus any fixed finite additional SCC bank: O_epsilon(X^(alpha+epsilon)). This sharpens the S322 ordered-word estimates 199/200 and 1999/2000 for these same grammars.
2. Arbitrary such sources, but only dyadic one-apex cores and a fixed finite dyadic bank in the reservoir: O_epsilon(X^(1/2+epsilon)).
3. A fixed source, or a source restricted to transitive chains, and the same dyadic reservoir: O_epsilon(X^(beta+epsilon)).

The finite bank in (1) need not be dyadic. In (2)-(3), every additional reservoir letter must have a verified dyadic Q2. Source finite-bank letters can be arbitrary. For transitive sources the weights are 3^ell, whose series converges for every positive s. None of these exponents is asserted sharp.

## 4. Exact fixed-source inversion without factoring D

Fix a literal nonempty source, and put h=H(U), A=A(U)>h. As in S322/S323, a dyadic reservoir has

    N=A*2^m-hD,       D=H(V)<=2^(m-1).

The target determines the unique possible m: the first integer with A*2^m>N, with m>=1. Check the strict shell and require

    D=(A*2^m-N)/h

positive, integral, odd and at most 2^(m-1). Rejection is exact for the grammar.

Now perform the following **factor-free** procedure.

- Generate every dyadic strong one-apex letter of cost c=log2 Q2<=m using the S323 paired-end compiler. Its chain length is bounded by 2+(m-1)/log2(3) using (2). Add any chosen finite verified dyadic bank. All cataloguing is finite.
- Keep one witness of least cost for each H if desired. Higher-cost copies of the same H are unnecessary because singleton padding supplies every unused unit of dyadic cost.
- Fix a total order on the retained letters. Enumerate their multisets of total cost at most m, maintaining the product of their H values. Stop successfully at product D, and add m minus the used cost singleton components. The empty multiset is allowed and gives D=1.
- If none matches D, reject. No factorization or divisor enumeration is performed. Optional tests that the current H product divides D are safe but not needed for the complexity bound.

Every enumerated object is physical. Conversely any successful reservoir gives such a multiset, because each SCC Q2 divides the dyadic product and hence is dyadic, and the named grammar specifies the permitted SCC alphabet. Replacing a letter by a cheaper witness of the same H preserves the target and adds only neutral singleton padding. Thus this inverse is necessary and sufficient for the full grammar, not just a finite factor subbank.

### Complexity including discovery

By (2)-(3), letter discovery through cost m costs 2^(beta m) poly(m). The Euler product bounds the number of multisets of total cost at most m by O_epsilon(2^((beta+epsilon)m)). A depth-first enumeration in nondecreasing letter order explores exactly those prefixes; each successful edge creates another admissible multiset, and the first failed cost test terminates that sorted loop. The number of components and the bit lengths of every accumulated integer are O(m). Thus, with polynomial factors absorbed into epsilon, the complete inverse takes

    O_epsilon(N^(beta+epsilon))                       (4)

arithmetic operations for a fixed source and fixed finite extra bank. This includes discovering the infinite alphabet below the target bound. There is no separate O(sqrt D) trial-factoring cost. The constant may depend on epsilon, the source, and the finite bank.

The algorithm enumerates physically possible products instead of attempting to factor an arbitrary D. It may be much slower than factorization on particular inputs; (4) is a worst-case guarantee for this restricted construction problem, not an integer-factoring algorithm.

## 5. A source-free exact compiler

Allow U itself to be an arbitrary ordinal sum of strong one-apex letters, singleton letters, and a fixed finite source bank. Retain the dyadic reservoir grammar above. Given an absolute target N, enumerate all multisets of the two-role letter alphabet satisfying

    product(source Q3) * product(reservoir Q2) <= N/2.

The source catalogue needs only Q3<=N/2 and is generated from every word allowed by (2); its cost is O(N^(1/2) poly(log N)). The reservoir catalogue is generated by the carry compiler in O(N^beta poly(log N)). Enumerate the multiset pairs by nondecreasing tagged-letter order and the product bound, discarding empty roles. For every pair calculate the exact triples, test the two values of gamma, and output a literal quotient inflation exactly when its count is N.

By (1), every realization in the grammar appears; no successful profile pair is excluded. The multiset bound in section 3 gives the complete compiler an arithmetic-operation bound

    O_epsilon(N^(1/2+epsilon)).                       (5)

The catalogue and traversal use integer bounds; real logarithms above only express their size. Actual witness order is O(log N), since the size of each allowed letter is bounded by a constant times its log weight. This algorithm does not assume a supplied source, supplied factorization, or an a priori list of all strong dyadic cores.

## 6. Physical halving works on either module

This statement is broader than the algorithmic alphabet: U,V may be **arbitrary** nonempty tournaments in either selective quotient. If one module has at least two SCCs, remove any one entire SCC and keep the other module unchanged. The resulting chart is still strong and its target N' satisfies

    0<N'<N/2.                                        (6)

For a reservoir component with H=d and Q2=q, write the old reservoir coordinates as D,B. Then

    N=A B-hD,        N'=A B/q-hD/d,
    N-dN'=AB(1-d/q)>0.

A nontrivial strong component has odd H>=3 and q>=2d, so N>3N'. For a singleton d=1,q=2, one has N-2N'=hD>0.

For a source component put a=H(C), b=Q_C(2), c=Q_C(3); let the remaining source have coordinates (u,q,r), and the unchanged reservoir have (w,y)=(H,Q2). Then

    N=6byq+gamma cyr-a(3y+w)u,
    N'=6yq+gamma yr-(3y+w)u,
    N-aN'=6(b-a)yq+gamma(c-a)yr>0.

For a nontrivial component a>=3, using b>=2a,c>=3a gives N>3N'. A singleton has (a,b,c)=(1,2,3), and direct subtraction gives

    N-2N'=(3y+w)u+gamma yr>0.

Stop before a module becomes empty. Every deletion is a literal induced deletion of a whole SCC inside the same named strong quotient. This extends S323's one-sided reservoir descent to either module, without imposing dyadic responses. It does not give a scalar operation determined by H alone: the selected physical components and their profiles remain part of the state, exactly as required by the S323 scalar barrier.

## 7. Verified examples and limits

The replay constructs the H=19,Q2=256 reservoir from target 3821 with singleton source in the gamma=2 chart, without factoring its required D. The gamma=6 chart selects the same reservoir plus two singleton pads for

    27629=7*3947,

with 3947 prime. It also selects the H=5,Q2=32 reservoir for the prime 859. Trial division verifies these small primality statements. They demonstrate both infinite atom lanes at particular targets; they are not new claims of membership or all-prime coverage.

The replay checks the exact inverse against a direct multiset enumeration on bounded budgets, tests source-free selection, compares literal quotient path counts, and verifies source and reservoir deletion on actual graphs. Its checks do not assume the conjecture or infer an infinite dyadic-core classification from the observed small bank.

A separate exploration found no further one-apex dyadic cores through chain length37 beyond the profiles (5,32) and (19,256). That finite observation is not a dependency of any theorem or algorithm above and is not an all-order classification. The exact compiler must continue to admit new larger letters if the carry test finds them.
