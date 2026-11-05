# An entropy budget for tournament constructors

Session 197, 2026-11-05. **All-order construction-class restrictions, with an exact bounded-image certificate. The full Hamiltonian-path spectrum remains open.**

The new-in-session consequence is a bridge from the growth rate of a regular language to the number of **distinct integer Hamiltonian-path counts** its physical tournaments can generate, even after arbitrary ordered joins. A genuinely exponentially branching language can still be too sparse. The concrete Fibonacci controller has exponent

    rho = 3 log(phi)/log(5) = 0.8969811534981604..., phi=(1+sqrt(5))/2,

and therefore has density-zero scalar support after arbitrary multiplication. It misses H=11 at every order, although a single forbidden transition produces an actual H=11 tournament. No literature-priority claim is made.

## 1. Source audit and exact advance

The starting guide is [S196](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/00-navigation/S196.md). The following mechanisms are inherited, not presented as this session's discoveries:

1. [S191 fixed-construction sparsity](https://github.com/eliottcassidy2000/math/blob/8d04f8557ea50c985521f039d5f58aa335617c9a/05-knowledge/results/session191_tournament_growth_20261104/two_apex/SPARSITY.md) proves a Dirichlet-product sparsity theorem for polylogarithmically many generator values, then applies it to fixed feedback-arc budgets, fixed cloning templates, and finitely many constant-run patterns.
2. [S192 fixed-apex tensor flags](https://github.com/eliottcassidy2000/math/blob/8d04f8557ea50c985521f039d5f58aa335617c9a/05-knowledge/results/session192_surprising_approaches_20261105/tensor/PROOF.md) already treats arbitrarily many fixed apexes, all mixed-period eigenvalues, sharp Jordan degree bounds, and the physical periodic families `(12)^m` and `(120)^m`. Independently rediscovering those facts during this session did not make them new.
3. [S194 periodic counters](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/two-apex/PROOF.md) supplies stronger short-period, all-prefix nonvanishing obstructions. The present argument does not replace those exact affine-return results.
4. [S196 low-cycle transfer](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/lowcycle/PROOF.md) already uses nonnegative saturation for bounded scalar images at arbitrary word length. General-k saturation is an immediate extension, not the main advance here.
5. Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, [DOI](https://doi.org/10.37236/1141), [paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), Theorem 1, gives the sharp residue-class lower bounds. We need only their consequence

       H(C) >= 5^((|C|-1)/3)                         (1)

   for a nontrivial strong tournament C. The official article and its theorem were checked. This is an external theorem, not proved anew here.

The delta is: allow **positive exponential language growth**, obtain a quantitative scalar-support exponent, and close the weak-output/singleton-padding loophole by factoring a regular language through actual strong components. The Fibonacci example lies outside all fixed-run classes: it has unbounded alternation and exponentially many words.

## 2. Physical fixed-apex presentations

Fix a tournament A on k named apexes. A word w over the alphabet `{0,1}^k` defines T_A(w): its additional vertices x_1,...,x_m form a transitive chain, with x_i->x_j for i<j, and bit a of w_i is one exactly when apex a points to x_i. All original apex arcs remain those of A. The order is k+|w|.

Different words or different choices of apexes need not describe different unlabelled tournaments. Our upper bounds count presentations, so these duplications can only help. In particular word inventories must not be read as coefficients of A000568.

A language L is regular if accepted by a finite automaton. Assume its number l_m of length-m words satisfies

    l_m <= C (1+m)^d lambda^m, lambda>=1.             (2)

This is the usual polynomial-times-exponential bound for a regular language; for a trimmed deterministic automaton, lambda may be taken as the largest spectral radius of its reachable and co-reachable transition blocks. The proof below only needs (2), plus regularity for the factor-closure step.

Let G(A,L) consist of all T_A(w), w in L, with no strongness assumption. Let J(A,L) also permit arbitrary ordered joins of these tournaments and arbitrary transitive singleton padding. Every operation here is on actual tournaments.

## 3. Strong components are projected contiguous factors

**Lemma 1.** Every nontrivial strongly connected component C of T_A(w) contains at least one apex. Its chain vertices, if any, form one contiguous interval of the chain. Consequently C has a presentation

    C = T_(A[S])(pi_S(u)),                            (3)

where S is a nonempty subset of the apexes, u is a contiguous factor of w, and pi_S deletes all other apex bits from each letter. Empty u is permitted.

**Proof.** A tournament induced only on chain vertices is transitive and has only singleton strong components. Thus nontrivial C contains an apex.

If x_i and x_j are in C with i<j, then for i<l<j we have x_i->x_l->x_j. There is a directed path from x_j to x_i inside C, so x_l is mutually reachable with C and must lie in it. Hence the chain portion is an interval. Restricting the actual apex arcs and incidence bits gives exactly (3). This is an induced-subgraph identity, not merely equality of counts. QED.

**Lemma 2.** If L is regular and satisfies (2), its contiguous-factor language Fac(L) also satisfies a bound of the same form, with the same lambda and d and possibly a larger constant. Each letterwise projection pi_S(Fac(L)) has no more words of any given length.

**Proof.** Trim an accepting automaton to the states lying on accepted paths. For every state p fix a path alpha_p from an initial state to p; for every state q fix a path beta_q from q to an accepting state. Their lengths have a common finite bound D.

Any factor u occurring on an accepted path labels a path from some p to some q. For fixed (p,q), the word alpha_p u beta_q is accepted. This map is injective on u of a fixed length, since its prefix and suffix are fixed. There are finitely many state pairs, so

    #Fac(L)_m <= sum_(p,q) l_(m+|alpha_p|+|beta_q|)
                <= C' (1+m)^d lambda^m.

This also handles factors of length zero whenever L is nonempty. Projection cannot increase the number of distinct words. QED.

Regularity matters here. A sparse arbitrary language can contain very long words carrying extremely many short factors, so an accepted-word count alone would not justify this lemma for an arbitrary nonregular language.

## 4. The entropy-to-spectrum theorem

For any set A of distinct integers at least two write <A> for its finite multiplicative closure, including one.

**Lemma 3 (positive-exponent Dirichlet extension).** Suppose

    # {a in A : a<=X} = O((1+log X)^e X^rho), rho>=0. (4)

Then, for every epsilon>0,

    # {n in <A> : n<=X} = O_epsilon(X^(rho+epsilon)). (5)

**Proof.** For every real s>rho, the series sum_(a in A) a^(-s) converges. To see this directly, split A into dyadic intervals; the jth interval contributes at most a constant times `(1+j)^e 2^(-(s-rho)j)`.

The Euler product

    Z_A(s) = product_(a in A) (1-a^(-s))^(-1)

therefore converges, since `-log(1-a^(-s)) <= a^(-s)/(1-2^(-s))`. Its expansion counts finite exponent vectors indexed by the **distinct integer generator values**. Different exponent vectors can yield the same product; this only overcounts. Monotone limits of finite products give

    sum_(n in <A>) n^(-s) <= Z_A(s).

Thus `#(<A> intersect [1,X]) <= Z_A(s) X^s`; set s=rho+epsilon. We do not use the geometric series of ordered construction words, which would count the wrong objects. QED.

**Theorem 4 (regular-controller entropy barrier).** Under (2), put

    rho = 3 log(lambda)/log(5).                       (6)

The number of distinct Hamiltonian-path counts at most X among J(A,L) is

    O_epsilon(X^(rho+epsilon))                       (7)

for every epsilon>0. The same result holds for a finite bank of fixed-apex presentations and regular languages, using the largest rho in the bank. In particular, if

    lambda < 5^(1/3),                                (8)

the scalar support has natural density zero and omits infinitely many odd primes, even after all the allowed ordered joins and padding.

**Proof.** By Lemmas 1 and 2, every nontrivial strong component of every input graph belongs to a finite bank of actual strong tournaments

    T_(A[S])(u), u in pi_S(Fac(L)), S nonempty.

This includes apex-only components with empty word. Let B be the set of distinct H-values of that strong bank. They are at least three. If h in B is at most X, Busch's bound (1) gives

    |C| <= 1 + 3 log(X)/log(5).                       (9)

For fixed S, chain length differs from total order by the fixed constant |S|. Counting words with (2), and summing over lengths through (9), yields

    # {h in B:h<=X} = O((1+log X)^(d+1) X^rho).       (10)

The extra logarithmic power is harmless and covers lambda=1 as well. There are finitely many apex subsets.

The strong components of a tournament are uniquely ordered. Every Hamiltonian path traverses them in that order, so H is the product of their individual H-values; singleton components contribute one. The same is true after arbitrary ordered joins. Hence every count in J(A,L) lies in <B>. Lemma 3 proves (7).

If rho<1, choose epsilon with rho+epsilon<1 to obtain density zero. Also the proof gives `sum_(n in <B>) 1/n < infinity`, since 1>rho. Were all but finitely many primes in <B>, this would contradict Euler's divergence of the prime reciprocal sum. All counts here are odd, so infinitely many omitted primes are odd. No prime number theorem is needed. QED.

**Scope.** Condition (8) is a sufficient obstruction, not a sharp necessary-and-sufficient threshold. Larger language entropy does not prove coverage. The one-apex full binary language has lambda=2 but still cannot realize H=13, as proved in the [inherited one-apex packet](https://github.com/eliottcassidy2000/math/blob/8d04f8557ea50c985521f039d5f58aa335617c9a/05-knowledge/results/session191_tournament_growth_20261104/one_apex/PROOF.md). Neither a modular image nor an exponentially large set of presentations certifies an integer support theorem.

## 5. A Fibonacci controller with a physical escape

Take two apexes a->b and only the two letters

    1: a->x->b,           2: b->x->a.

Let L consist of all words over {1,2} with no consecutive 2s. Its usual two-state controller records whether the previous letter was 2. The transition matrix is

    [[1,1],[1,0]],

so its length-m inventory is F_(m+2), and its growth constant is phi=(1+sqrt(5))/2. This is a genuinely exponentially branching language. It includes arbitrarily many alternations, such as `(12)^m`, so no fixed constant-run bound covers it.

Since `phi^3=2+sqrt(5)<5`, Theorem 4 applies with rho in (6) equal to 0.8969811534981604... . In particular its entire ordered-join closure has O(X^(9/10)) distinct counts through X. This exponent can be certified without decimal logarithms: phi^10=55phi+34=(123+55sqrt(5))/2<125, since 5*55^2<127^2. Therefore rho<9/10. The exponent 9/10 is a convenient consequence, not an optimized estimate.

### A quantitative partition-function certificate

For any real s>rho set z=5^(-s/3), so z<1/phi. The strong-word generating function is

    sum_(m>=0) (F_(m+2)-1) z^m
      = z / ((1-z)(1-z-z^2)).                       (11a)

If B is the set of distinct strong generator values, choose one presentation for each value and then overcount by all strong words. Busch gives

    sum_(b in B) b^(-s)
      <= z^2 / ((1-z)(1-z-z^2)).                    (11b)

Consequently the actual product support M satisfies the explicit bound

    # (M intersect [1,X]) <= exp(C(s)) X^s,
    C(s) = z^2 / ((1-3^(-s))(1-z)(1-z-z^2)).         (11c)

This is the same critical inequality as convergence of the weighted Fibonacci transition matrix z[[1,1],[1,0]]: its spectral radius is z*phi<1. It connects a finite controller's transfer matrix to convergence of an arithmetic Euler product. The constant is not claimed numerically sharp, especially near s=rho; distinct-word collisions only improve it.

### Strongness without a hidden assumption

A word consisting solely of 1s gives the transitive order a,chain,b, so H=1. Every word containing 2 is strong. Indeed b reaches a through any 2-vertex, while a->b is fixed. The apex component reaches the first chain vertex: a does so for letter 1 and b for letter 2. The last chain vertex reaches the apex component: it points to b for letter 1 and to a for letter 2. The increasing chain then places all chain vertices in the same strong component.

### The extra all-order omission H=11

If such a strong word had H=11, (9), or the exact inequality `5^(n-1)<=11^3`, forces n<=5, hence word length at most three. The complete valid word inventory at these lengths gives:

    epsilon, 1, 11, 111: H=1;
    2: H=3;
    12, 21: H=5;
    112, 121, 211: H=9;
    212: H=13.

There is no 11. Arbitrary products cannot create the prime 11 from other positive integer factors. Thus H=11 is absent from the entire controller and ordered-product closure at every order.

This omission is purely architectural. Permitting the single forbidden transition in word `122` produces an actual strong tournament with

    outgoing rows = [6,24,26,17,1], H=11.             (11)

Rows use vertex order a,b,x_1,x_2,x_3, with bit j of row i meaning i->j. The verifier checks the actual adjacency matrix and counts its Hamiltonian paths independently by subset dynamic programming.

The allowed word `212` already gives H=13. Thus this controller genuinely crosses the inherited one-apex H=13 barrier while still failing at the smaller target 11. A construction class need not acquire targets in numerical order as its structural power changes.

### Exact arbitrary-order image through 5,001

Busch gives the integer cutoff

    5^(n-1)<=5001^3  => n<=16.

All nontransitive Fibonacci outputs below that H-bound therefore have chain length at most 14. There are exactly `F_18-2=2582` valid words through that length, of which 2567 contain a 2 and are strong. Exhaustive exact transfer gives:

- 157 distinct strong-generator H-values at most 5001;
- 158 values after also including the transitive count one;
- 414 distinct values after arbitrary ordered products, including one.

These are **all-order results inside this stated language/product class**. The finite word cutoff is deduced from an all-order theorem; it is not guessed from stabilization. `bounded_image.json` records a literal word for every generator value and a factorization into generator values for every product value. Its first omissions include 7,11,19,21,23,31,35,37,41,43. They are not global tournament exclusions.

## 6. Three useful boundary tests

### Singleton padding invalidates a naive order cutoff

With the same two apexes, the word `0^m 2` has H=3 for every m>=0. All 0-vertices are singleton source components, followed by one directed triangle. Its order is unbounded. Applying (9) directly to arbitrary weak words would therefore be false. Theorem 4 deliberately applies Busch only after the actual SCC factor/projection step.

### Bounded mixed-period templates are still sparse

A fixed finite collection of templates

    w_0 v_1^(m_1) w_1 ... v_r^(m_r) w_r, m_i>=0,

where all displayed words are fixed and each v_i is nonempty, is regular and has only polynomially many length-n presentations: every exponent is at most n, so at most `(n+1)^r` exponent vectors suffice per template. The periods may themselves use multiple letters, and their repetitions may create unbounded alternation. Theorem 4 with lambda=1 gives X^(o(1)) scalar support even after arbitrary products, including weak outputs.

This extends the inherited constant-run-template conclusion to fixed bounded expressions of mixed blocks. It still says nothing comparable about unrestricted regular languages with sufficiently large growth.

### Any finite scalar census can be hidden in zero entropy

Let W be any finite bank of actual two-apex witness words over {0,1,2,3}. The infinite regular language

    union_(w in W) 0* w 3*

has only O_W(m+1) words of length m. Prefix 0-vertices are universal source components and suffix 3-vertices universal sink components, so every output has exactly H(T(w)) for some w in W. Thus an arbitrarily impressive **finite** coverage record can be preserved in an infinite zero-entropy constructor that adds no new scalar values at all. This is a physical padding construction, not a formal language analogy. Entropy and eventual coverage require their own all-order arguments.

### Exponentially many A000568 objects can share one scalar value

Let C be the directed triangle, H(C)=3, and D the strong four-tournament `[8,1,3,6]`, H(D)=5. Choose any ordered sequence containing exactly m copies of C and m copies of D, and orient all arcs forward between components. The unique ordered SCC decomposition distinguishes all

    binom(2m,m)

isomorphism classes. Every graph has order 7m and H=15^m. Thus exponentially many ordinary classes can sit in one arithmetic fibre.

Their automorphism groups have order 3^m: each triangle rotates, D is rigid, and the ordered components cannot be permuted. The inherited path-orbit formula therefore gives 5^m Hamilton-path orbits per class. These are explicit examples explaining why A000568 growth, presentation growth, and distinct scalar growth are different quantities. They are deductions from inherited SCC/orbit facts, not a new general multiplicity theorem.

## 7. Reproduction and verification scope

Run in this directory:

    python3 verify.py > normal.stdout
    cp receipt.json receipt.normal.json
    cp bounded_image.json bounded_image.normal.json
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    cmp receipt.normal.json receipt.json
    cmp bounded_image.normal.json bounded_image.json

No third-party packages or repository data are needed. All tests use explicit exceptions and remain active under optimization. The checker independently generates the physical fixed-apex graphs and the inherited nonnegative slot transfer.

It performs 1530 complete word/graph comparisons over all labelled apex tournaments of orders 0 through 3 in the specified finite word universes. These also verify the actual SCC interval/projection identity and saturation compatibility. The two-apex Fibonacci census covers all 2582 required words, with 160 further graph-side checks, including deterministic longer words. Exact integer powers, not floating logarithms, establish all order cutoffs. The 11 escape and the weak-padding hostile are checked separately.

The general entropy, factor-closure, SCC, and Dirichlet statements rest on their proofs above. Finite checks do not establish them by extrapolation. The sharp Busch lower bound is the only imported quantitative theorem. No dirty checkout or unnamed repository path was read or changed.
