# Parity-context divisors and a prime two-orbit 7p construction

Session197, 2026-11-05. **All-order exact context theorem and an explicitly certified physical counterexample.** The full Hamilton-path spectrum remains OPEN. No literature-priority claim is made.

## 1. Sources and the two directions of the question

Write H(T) for directed Hamiltonian vertex words. The [S196 symmetry proof](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/symmetry/PROOF.md) proves that Aut(T) acts freely on Hamilton words, so its order divides H(T), and gives the exact path-pointed mass identity

    sum_[T], |T|=n H(T)/|Aut(T)| = 2^binom(n-1,2).

It also proves prime rigidity and the sharp first nonrigid7p value4151. The [S196 orbit-spectrum proof](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/orbit-spectrum/PROOF.md) classifies that equality stratum. Its unique moving seven-orbit argument is explicitly bounded in H; it does not claim that every symmetric7p graph has a seven-vertex module.

We use its inherited odd-cycle collection formula (OCF):

    H(T) = sum_Gamma 2^|Gamma|,

where Gamma runs over collections of vertex-disjoint directed odd cycles of length at least3, counted modulo cyclic rotation, including the empty collection. The labelled-bin extension below is reproved directly from OCF, not asserted to be a newly discovered polynomial.

The [S193 rooted-clock profile proof](https://github.com/eliottcassidy2000/math/blob/169792fb5eeb483e351112b5a09706e6d7bd89f0/05-knowledge/results/session193_triangulated_proof_20261105/profile_substitution/PROOF.md) defines the complete labelled-bin polynomial and proves an input-content obstruction: fix an input S, then its profile gcd divides every whole-module continuation. The [S193 block-divisibility proof](https://github.com/eliottcassidy2000/math/blob/169792fb5eeb483e351112b5a09706e6d7bd89f0/05-knowledge/results/session193_triangulated_proof_20261105/block_divisibility/PROOF.md) additionally treats exposed old blocks and apex budgets. The [S195 path-cover algebra](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/pathcover-algebra/PROOF.md) supplies general typed composition. The [S196 transitive-cloning proof](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/spectral-moves/PROOF.md) supplies the nonnegative maximal-run coefficients and their exponential expansion.

Our different question fixes the OUTER context and allows every nonempty tournament as its input, subject only to order parity. We identify its exact common divisor by finitely many physical transitive substitutions. This is distinct from fixing the input and varying the contexts. A separate construction shows that the symmetric7p sector can contain module-indecomposable graphs with two moving orbits, even when both orbit tournaments are Paley type.

## 2. Positive resolution into transitive profiles

Let f_k(T) count unordered spanning covers by k nonempty directed paths. Put

    B_T(d) = sum_k (d)_k f_k(T),

where (d)_k is the falling factorial. For nonnegative integer d, this counts d labelled path bins, empty bins permitted. The transitive j-tournament has B_Tr_j(d)=d^j.

**Lemma1.** For every nonempty n-vertex tournament,

    B_T(d) = sum_Gamma 2^|Gamma| d^j(Gamma),
    j(Gamma) = n - sum_(C in Gamma)(|C|-1).                 (1)

Every j is positive and congruent to n modulo2. The empty collection gives d^n with coefficient1; every other coefficient is even and nonnegative.

**Proof.** For each assignment of vertices to d bins, the number of path choices is the product of the H-values on the bins. Expand each by OCF. A fixed global disjoint odd-cycle collection Gamma can occur precisely when all vertices of each of its cycles occupy the same bin. Treat each such cycle as one object and every uncovered vertex as one object. There are j(Gamma) objects, independently assigned to d bins, giving d^j choices. Its weight is2^|Gamma|. This proves (1) at all nonnegative integers and hence as a polynomial. Disjoint cycles have odd length, so every subtraction |C|-1 is even. At least one object remains when n>0. QED.

Equivalently, the untyped path-cover polynomial sum_k f_k z^k is the corresponding positive combination of Touchard polynomials sum_k S(j,k)z^k. The equivalence follows by the invertible falling-factorial basis change. This is an exact counting identity, not a claim that T is isomorphic to a union or mixture of transitive graphs.

**A small modular consequence.** If n is odd and ell is an odd prime, then ell divides f_(ell-1)(T). Indeed, for every positive odd j,

    (ell-1)! S(j,ell-1)
      = sum_(i=0)^(ell-1) (-1)^(ell-1-i) binom(ell-1,i)i^j
      = sum_(i=1)^(ell-1)i^j = 0 mod ell.

The last equality follows by pairing i with ell-i, since j is odd. Also binom(ell-1,i)=(-1)^i mod ell. The factorial is invertible modulo ell. Apply (1) and the falling-factorial expansion. No analogous assertion for even n is made.

## 3. Exact parity-context divisor

Fix a nonempty rooted tournament (Q,r). Replace r by any nonempty tournament T, uniformly orienting all cross arcs as r did. Write C[T] for this physical substitution. Let a_k be the number of legal words containing each exterior vertex once and k indistinguishable r-symbols, with no consecutive r-symbols. Thus a_1=H(Q)>0 and only finitely many a_k are nonzero. Put D=max{k:a_k>0}, so D>=1.

Cutting at the boundaries of maximal T-runs gives the inherited exact formula

    H(C[T]) = sum_(k=1)^D a_k k! f_k(T).                   (2)

Consequently the map from B_T to H(C[T]) is a fixed integer linear functional in its falling-factorial coefficients. Applying it to (1) gives

    H(C[T]) = sum_Gamma 2^|Gamma| h_j(Gamma),
    h_j = H(C[Tr_j]), j>=1.                               (3)

All terms are positive physical counts. We never use j=0, and formal clone length zero is NOT interpreted as deletion.

For e in {1,2}, interpret e=1 as odd input order and e=2 as positive even input order. Define

    G_e(C) = gcd { H(C[T]) : |T|>0, |T|=e mod2 }.

**Theorem2 (finite exact divisor test).**

    G_e(C) = gcd(h_e,h_(e+2),...,h_(e+2(D-1))).            (4)

Thus the right side divides every count from an arbitrary input module of the prescribed parity, and is the largest universal such divisor.

**Proof.** By (3), the gcd of the h_j of that parity divides every H(C[T]). Conversely every Tr_j is a permitted physical input, so the gcd over all inputs divides every such h_j. These two infinite gcds are equal.

The Stirling expansion gives, for every positive j,

    h_j = sum_(k=1)^D a_k sum_(b=1)^k (-1)^(k-b)binom(k,b)b^j
        = sum_(b=1)^D beta_b b^j,                         (5)

with integer beta_b. The positive-parity subsequence u_t=h_(e+2t), t>=0, is a linear combination of the D sequences (b^2)^t. It satisfies the monic integer recurrence with characteristic polynomial

    product_(b=1)^D (x-b^2).

Every term after the first D is therefore an integer linear combination of earlier terms. The gcd of the first D divides the entire sequence, while the reverse divisibility is immediate. QED.

This is an exact arithmetic decision procedure for the universal divisor, not an inverse construction of every multiple of it. It does not require all positive combinations in (3) to be realizable. The included transitive inputs and the valid expansion in one direction suffice for the gcd equality.

A general prime-modulus implementation can work directly with (5): group bases by their squares modulo ell. The odd subsequence vanishes identically modulo ell exactly when each nonzero square class has sum beta_b b=0; for the positive even subsequence the criterion is sum beta_b b^2=0. This follows by a Vandermonde matrix on distinct nonzero squares. Bases divisible by ell vanish at every positive input length. This optional criterion gives no additional realizability claim.

## 4. Concrete divisors, exact boundaries, and atom restrictions

Use the cyclic carousel C_q on Z/q with differences1,...,(q-1)/2, and root0. Direct legal-run enumeration gives:

- C3: a=(3,1), h_j=2^j+1, G_odd=3, G_even=1.
- C5: a=(15,13,2), h_j=2*3^j+7*2^j-5, G_odd=15, G_even=1.
- C7: a=(175,227,98,6), h_j=6*4^j+74*3^j-31*2^j-9, G_odd=5, G_even=1.

These are finite computations of complete template word sets followed by Theorem2, not numerical extrapolations. The verifier also checks the formulas on actual cloned graphs and arbitrary input graphs.

In particular:

1. Wrapping any odd-order T in the cyclic triangle C3[T,1,1] gives a multiple of3. Divisibility alone excludes prime targets other than3 and7p targets with p other than3. Here the H3 boundary is exactly singleton T. Also H21 is absent in this wrapper without invoking the global H21 theorem: input order1 gives3; input order3 gives9 or15; input order at least5 gives H>=h_n=2^n+1>=33 by (3).
2. Substituting any odd-order T at a carousel-five vertex gives a multiple of15 and at least15, so no numerical prime or7p target occurs.
3. Substituting any odd-order T at a carousel-seven vertex gives a multiple of5 and at least175, so no numerical prime or7p target occurs: the only prime divisible by5 is5, and the only7p with p prime divisible by5 is35.

The even-input hypothesis changes the conclusion: Tr2 substitutions give H=5,41,629 respectively, and the exact common divisors are1. For C5 the even-input H41 is prime. A common divisor of1 is not a coverage theorem.

There is also an exact extremal gap valid for every context C. Formula (3) gives H(C[T])>=h_n. If T is nontransitive, one directed triangle gives

    H(C[T]) >= h_n+2h_(n-2).

Equality holds exactly when T has a unique odd cycle, necessarily a triangle: all terms are positive, and every nontransitive tournament contains a triangle. The one-triangle SCC characterization is inherited from S196. Thus the lower bound is attained by a directed triangle and any n-3 ordered singleton padding, and the minimum is attained only by transitive T. This is a context-level comparison, not a graph embedding.

## 5. A universal profile switch on seven-orbit modules

For the carousel A=C7 and Paley B with differences{1,2,4}, exact cover vectors (f_1,...,f_7) are

    A: (175,819,1113,630,168,21,1),
    B: (189,819,1113,630,168,21,1).

Their labelled-bin polynomials are respectively

    B_A(d)=d^7+28d^5+112d^3+34d,
    B_B(d)=d^7+28d^5+112d^3+48d.

Therefore, in EVERY rooted whole-module context C,

    H(C[B])-H(C[A]) = 14 H(Q).                            (6)

This is an exact all-order physical switch: all internal paths using two or more module runs cancel, while a single run contributes the difference14 times H(Q). The broader positive resolution explains it, but the displayed finite cover vectors already certify (6). They are checked by independent cycle enumeration and subset path covers.

Paley's ordered-profile content is21, whereas the carousel's is7. Thus retaining an intact Paley module can never give7p with prime p other than3. The next section shows why this statement must not be extended from modules to arbitrary automorphism orbits.

## 6. Two synchronized Paley orbits give a module-indecomposable 7p graph

Let A_i,B_i, i in Z/7, be14 named vertices. Each of A and B induces the Paley tournament with difference set{1,2,4}. Set

    A_i -> B_j iff j-i is1 or3 modulo7.

Its exact outgoing rows, ordered A_0,...,A_6,B_0,...,B_6, are

    [1302,2604,5208,10289,4322,8517,651,
     2863,5726,11325,6394,12661,8939,1495].

**Theorem3 (explicit counterexample to an unrestricted module reduction).** This tournament T is strong, has no nontrivial proper module, and satisfies

    H(T)=15144269=7*2163467,
    2163467 is prime,
    Aut(T)=C7.

Its two seven-element automorphism orbits are precisely A and B, and neither is a module.

**Strongness.** Each Paley block contains the directed cyclic successor1 cycle. Cross arcs occur in both directions, so the two internally strong blocks reach each other.

**Symmetry and full group.** The permutation

    (0 1 2 3 4 5 6)(7 8 9 10 11 12 13)

preserves every difference defining an arc, hence is an automorphism of order7. The Hamilton count and primality are exact finite certificates detailed below. Free action gives |Aut(T)| dividing7*2163467. The group also embeds in S14, and the prime2163467 is greater than14, so it cannot divide |Aut(T)|. Hence |Aut(T)|=7 and the displayed rotation generates it.

**Module-indecomposability.** First every regular seven-tournament is module-indecomposable. If M is a module, equal total degrees force its induced tournament to be regular, so a nontrivial proper M has odd size3 or5. Size5 is impossible: an exterior vertex beating M would have outdegree at least5, whereas one losing to M would have outdegree at most1, both incompatible with degree3. For size3, its four exterior vertices are uniform on M. Each beating M must have no outgoing arc among those four, so at most one can do so. Each losing to M must beat the other three exterior vertices, so at most one can do so. Four exterior vertices cannot fit these two classes. Thus neither size occurs.

Now a module M of T intersects each Paley block in a module of that induced block. If it contains all of A but is proper, an omitted B vertex has mixed incidence with A (five wins and two losses), a contradiction; similarly for B. Thus any nontrivial proper M could contain at most one vertex from each block, hence would be a cross-pair. But every A vertex has outdegree5 and every B vertex outdegree8. Two vertices forming a module have identical incidences with every exterior vertex, so their total outdegrees differ by exactly1, not3. This excludes the last case.

The graph count is independently checked by forward endpoint subset DP and backward startpoint subset DP, with different state recurrences. A recursive Lucas primality certificate is included: complete factorizations of p-1, certified smaller primes, and a witness a satisfying a^(p-1)=1 modp and gcd(a^((p-1)/q)-1,p)=1 for every prime q dividing p-1. Independent trial division through floor(sqrt(p)) also checks the final prime. The graph's complete module universe is additionally exhausted as a control, though the argument above is sufficient.

The normalized path-orbit mass for this one ordinary class is exactly2163467, a prime, and its labelled mass is14!/7. These are the inherited free-action consequences applied to this new explicit construction. The example does not contradict the sharp4151 theorem. It refutes only the stronger, unsupported extrapolation that every symmetric7p witness has a single carousel module, or that a Paley seven-orbit always carries a factor3 into the ambient H.

## 7. Reproduction and scope

Run `python3 verify.py` and `python3 -O verify.py`, comparing the JSON receipts byte for byte. All gates are explicit exceptions. The verifier independently checks all labelled tournaments through order5, selected larger profiles, all complete run words for the named templates, parity-gcd recurrences, actual clone and nontransitive substitutions, the profile switch, and the14-vertex graph and primality certificate.

The all-order divisor theorem follows from (1)-(5), not from finite scans. The finite theorem in Section6 is fully explicit and independently certified; no claim of minimal14-vertex order or minimum H among two-orbit examples is made. The construction does not prove infinitely many primes in any parametric sequence. Neither absence of an obstruction nor the existence of a large prime-orbit quotient solves the unrestricted spectrum.
