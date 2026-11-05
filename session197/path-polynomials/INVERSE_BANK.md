# A target-bounded inverse algorithm for finite physical context banks

Session 197. **PROVED complete finite decision algorithm for a specified bank; FINITE-EXACT D/E/F/G prime coverage and explicit all-order omissions.** This does not prove the full `{7,21}` spectrum and does not claim a polynomial-time algorithm in the bit length of the target.

## 1. A finite bank, not arbitrary positive matrices

Fix finitely many actual rooted tournament templates (Q_i,r_i). Start from K1 and permit:

- ordered joins A->B;
- the whole-module substitutions C_i[T]=Q_i[r_i:T].

Let M=max_i(|Q_i|-1). Identity templates may be discarded. The full bin-profile operators are supplied by the inherited path-cover kernel or by the positive skeleton compiler in `UNIVERSAL_CONTEXTS.md`. The input T is an intact module at each substitution. No nonuniform exterior attachment or internal arc change is silently included.

A template is **active** if its root lies in a directed triangle. For an n-vertex input T, with h=H(T), an active output contains the actual induced subtournament C3[T,K1,K1]. Hamilton paths of an induced subtournament inject into those of the larger tournament by deterministic vertex insertion. Therefore

    H(C_i[T]) >= H(C3[T,K1,K1]) = h+B_T(2) >= 2^n.              (1)

The last inequality follows from the empty odd-cycle packing term x^n in B_T. This bound concerns the order of the entire actual input, including every padding vertex and internal join.

If the root lies in no directed triangle, it lies in no directed cycle: a shortest directed cycle through it of length greater than three can be shortened using a chord from the root. Thus the root is a singleton strong component. Substitution inserts the strong components of T at that position in the ordered component list, giving

    H(C_i[T])=c_i H(T),       c_i=H(Q_i-r_i)=H(Q_i).             (2)

Call such a template passive. The graph Q_i itself is constructible as C_i[K1].

## 2. Strip only a finished outer shell

**Theorem 1 (logarithmic witness-order bound).** Let P>1 be attainable in this fixed grammar, and let Omega(P) count its prime factors with multiplicity. There is a graph in the same grammar with count P and order at most

    floor(log_2 P)+M Omega(P).                                 (3)

For prime P the bound is floor(log_2 P)+M. For a target 7p, with p an odd prime, the unconditional bound is floor(log_2(7p))+2M; it does not need a prior proof that seven is unattainable.

**Proof.** Start at the root of an expression whose final count is P. Peel downward only through its **finished outer shell** of joins and passive substitutions. At a join, keep both scalar factors. At a passive substitution, record c_i if c_i>1 and continue into its input. Stop immediately on reaching any active substitution, retaining its entire subexpression unchanged. Discard count-one factors and leaves.

No padding inside a retained active input is deleted or moved. In particular this argument never transports a simplified profile across a later active wrapper. All later operations in the peeled shell act only by multiplication on H, by (2) and the ordered-join identity.

The remaining scalar factors h_1,...,h_t are greater than one and multiply to P. Each retained active subexpression has count h_j and, by (1), its original input has order n_j<=floor(log_2 h_j). Its output order is at most floor(log_2 h_j)+M. Each recorded passive factor c_i can be represented by the constructible seed C_i[K1]=Q_i, whose order is at most M+1<=M+floor(log_2 c_i). Ordered-join these representative graphs. This preserves the product P without any subsequent active wrapper.

Since t<=Omega(P), summing their orders gives (3). For prime P there is only one nonunit factor. QED.

The theorem is existential normalization of a **completed count witness**. It does not authorize removing neutral vertices during a live recursive construction, where their higher profile coordinates can change the next output.

## 3. The complete inverse decision algorithm

Given P, set N to (3), or use the sharper prime bound. Enumerate all reachable full bin profiles of degrees 1,...,N. Begin with B_K1=x. At each degree apply every bank operator whose output has that degree, and multiply all pairs of smaller-degree profiles whose degrees sum to it. Deduplicate equal polynomials, retaining one physical expression for each. Discard profiles whose H-readout exceeds P.

This enumeration is complete for the chosen order bound by induction on expression trees. Equal full profiles are interchangeable for these continuations because every bank operation is linear on the complete profile and joins multiply it. Deduplicating only H would be wrong. Count pruning is safe: every operation retains its input as an induced subtournament, and joins multiply positive counts, so no later output can fall below an ancestor's H.

If a retained profile has readout P, replay its stored expression into explicit adjacency rows. If none does, Theorem 1 proves P is omitted by this grammar, at **every** order, rather than merely beyond the searched cutoff.

For a fixed bank this is polynomial in the numeric target P, not necessarily in log P. After identity operations are removed, unary nodes add at least one vertex, and a binary-join expression with n vertices has O(n) nodes. There are at most C^n expression trees for a bank-dependent constant C. Here N=O(log P), since Omega(P)<=log_2 P. Profile degrees and coefficient bit lengths are polynomial in N; literal path-cover bounds or the odd-cycle expansion give O(N log N) coefficient bits. The number of states and pair combinations is therefore P^O(1). No small practical exponent or bit-polynomial complexity is asserted.

This is a genuine inverse decision procedure for a fixed physical language. It does not show that a desired target passes the procedure, and allowing arbitrary new Q_i without constructing them would simply reintroduce the original inverse problem.

## 4. Exact D/E/F/G experiment, with physical replay

Use the four proved port operators D,E,F,G from the companion packets, adding respectively 2,3,4,5 vertices. Thus M=5. The prime interval P<=8191 has complete witness-order cutoff

    N=floor(log_2 8191)+5=17.

The checker enumerates 12,553 distinct full profiles of degree at most 17 and H at most 8191. Of the 1,027 odd primes at most 8191, exactly **77** occur in this grammar. This is exact coverage of the specified prime interval for this bank, not merely an order-17 census: Theorem 1 supplies the global cutoff justification.

The output lists all 77 primes and replays a concrete graph for every one, verifying its Hamilton count by independent vertex-subset endpoint DP. The first primes omitted after the known seven are 29,31,41,43,53,59,61. No claim is made that any of these values is absent from tournaments.

Two explicit constructive outputs beyond the initial seeds are

    H(F(D(Tr_2))) = 197,
    H(D(X E(Tr_2))) = 211,

where X appends a universal source or sink. For the first, D(Tr_2) has even-state coordinates u_0=5,u_1=1, and the F readout is 25u_0+72u_1=197. For the second, E(Tr_2) has polynomial x^5+8x^3+2x; multiplying by x gives even coordinates u_0=11,u_1=13, and D gives 5*11+12*13=211. Both actual eight-vertex adjacency matrices are in the verifier output.

### An all-order obstruction to a returning ternary rule in this bank

The bank contains its H13 seed F(K1), but contains neither 41=3*13+2 nor 43=3*13+4. This can also be proved without the prime-interval search.

At any last nontrivial generator the eight possible scalar forms, classified by input parity, are

    D: 3h+6u,     5h+12u;
    E: 5h+18u,    11h+36u;
    F: 13h+36u,   25h+72u;
    G: 19h+60u,   39h+120u,

with h positive odd and u>=0. Since the target is prime, final count-one join factors can be discarded; eventually one of these forms must produce it.

For target 41 the only nonnegative integer possibilities are D with even input (h,u)=(1,3), and E with odd input (1,2). Every H1 tournament is transitive: a nontransitive tournament contains a triangle, and induced insertion gives H>=3. But these candidates would require B_T(2)=40 or 14 respectively, whereas a transitive input has B_T(2)=2^n. Neither is possible.

For target 43 the sole integer possibility is E with odd input (h,u)=(5,1). It would have B_T(2)=2h+6u=16, so (1)'s empty-packing bound gives n<=4. Odd parity forces n<=3, where tournament Hamilton counts are only one or three, not five. This is impossible. Hence both omissions hold over arbitrary expression size.

### The inherited whole-module barrier still applies

Even enlarging the bank with other whole-module contexts cannot send a fixed nontransitive input T literally to 3H(T)+2 or 3H(T)+4. The companion positive coordinates recover the inherited S194 obstruction transparently. For odd |T|>=3, the empty-packing monomial ensures u_1>=1, so D(T) has count at least 3h+6. For even nontransitive T, D(T) has count at least 5h>=3h+6. An active context dominates D(T) by (1), while a passive context gives an integer multiple of odd h>=3 and cannot equal 3h+2 or 3h+4.

This does not exclude the same numerical targets obtained from unrelated inputs, a change of root followed by single-vertex substitution, or an operation that splits the retained input's external neighborhood types. It identifies precisely why positive whole-module operators and a target-bounded search do not alone close the desired ternary induction.

## 5. Reproduction

Run `python3 verify_inverse_bank.py` and `python3 -O verify_inverse_bank.py`, comparing their stdout files. The standard-library checker uses explicit exceptions. It constructs the dynamic-programming profiles, verifies the two small all-order omission reductions by exact integer enumeration, replays all 77 prime witnesses as actual graphs, and checks their counts independently. The displayed 197/211 witnesses are separately constructed from the stated recipes, not selected only by matching a stored integer.

The order theorem, full-profile equivalence, and monotonicity justify what the finite computation certifies. They must accompany the table; an unqualified finite cutoff would not establish an all-order grammar omission.
