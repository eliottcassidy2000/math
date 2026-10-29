# Exact native weights, strip transport, and a smallest class-pair obstruction

Session105, 2026-10-29. **PROVED elementary identities and obstruction; FINITE-EXACT n2..6 census.** No priority claim for the mass formula or reversible quotient construction. This is an ordinary-isomorphism-class operator, not the repository's fundamental paired-black converse graph.

## 1. Provenance and the carrier correction

Recovered inputs:

* `01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md`: fixed-path words in ordinary class A have cardinality h_A=H_A/a_A, where a_A=|Aut(A)|. The automorphism action on ordered Hamilton paths is free. This is inherited, not rediscovered here.
* Session102 `metagraph/SIGN_GRAMMAR.md`: actual native operation, colour guard, SCC bridge-port defect; `CAYLEY_QUOTIENT.md`: switching aggregation is not equitable, and its connectedness does not lift.
* Historical `07-reflections/metagraph-as-transfer-chain.md` proposes a product of class extension matrices as the metagraph. Its missing path/port dependence cannot simply be filled by a class-pair transfer: §5 supplies an explicit obstruction. In particular its displayed “transitive to six neighbours” is not the native K_P support convention here: a transitive source has exactly one Hamilton path and hence one native neighbour.

For ordered directed Hamilton path P, K_P(T) keeps its path arcs and reverses every other arc. Its inverse is itself with the *same* directed P. A path is blue when positional reversal is an anti-automorphism, and black otherwise. Both directions of a native incidence have the same colour.

Let F_n be tournaments on positions 0,...,n−1 with consecutive arcs forward. Its beta=binom(n−1,2) free bits are the nonconsecutive pairs. Let c(x) be ordinary tournament class, and Cx its nonpath-bit complement. Then C is a fixed-point-free involution for n>=3. It is **not** tournament converse. All matrices below count directed incidences; a diagonal M_AA includes two words for each C-pair, not a single simple loop.

## 2. Three weights, not one

Define

    M_AB = #{x in F_n : c(x)=A, c(Cx)=B},
    W_AB = #{P in HP(T_A) : K_P(T_A) belongs to B}.

The normalization orbit map immediately gives

    W_AB = a_A M_AB,
    M_AB = M_BA,
    sum_B M_AB = h_A,
    sum_B W_AB = H_A.                                  (1)

Each identity separately holds with the blue or black filter, replacing the total degree by that colour's degree. The involution C proves symmetry, and freeness of Aut(T_A) on paths proves the multiplicity a_A. Equivalently the stabilizer of the *directed pair* (T,K_P(T)) under simultaneous relabelling is trivial: their agreement arcs recover the spanning directed path P, whose every vertex is distinguished by its path position.

Thus the support graph remembers only whether M_AB>0; M is the symmetric path-orbit incidence matrix; W is the generally nonsymmetric labelled-witness adjacency from a specified source representative. The actual aggregate number of labelled incidences between all representatives of A and B is

    (n!/a_A) W_AB = n! M_AB.                            (2)

In particular W_AB/a_A=W_BA/a_B. This is exact detailed balance with orbit-size measure, not permission to discard a_A. At n3, W_(transitive,cycle)=1 and W_(cycle,transitive)=3, while both M entries are1. The labelled graph remains the two K1,3 components from Session102, even though its ordinary class quotient is one support edge. Ordinary class quotient connectedness does not imply labelled connectedness either.

## 3. An authentic Markov operator and safe spectral bound

If a labelled source chooses a uniformly random directed Hamilton path and performs its native move, ordinary isomorphism classes form an equitable Markov quotient:

    P_AB = W_AB/H_A = M_AB/h_A,
    pi_A = h_A / 2^beta.                               (3)

Relabelling transports guards and all path choices, which is the reason this quotient is valid. It is not true for switching classes. Let D=diag(h_A), and let J map a class vector to the word vector constant 1/sqrt(h_A) on its fibre. Then J*J=I and

    S = D^(-1/2) M D^(-1/2) = J* C J.                  (4)

Here C also denotes the permutation operator of native complement on F_n. Consequently S is symmetric and -I <= S <= I. Equivalently, for every real f on classes,

    |sum_AB M_AB f_A f_B| <= sum_A h_A f_A^2.           (5)

Equality conditions are literal: a lifted fibre-constant function must be C-even or C-odd respectively. Powers of P describe actual successive random native moves with a *fresh uniform Hamilton path at each visited source*. They do not describe iterating one named path, which simply returns after two moves. The projection/resampling J J* between complements must not be suppressed.

Converse does preserve this Markov quotient: it transports P to reverse(P). If classes are merged into converse orbits U, the aggregate symmetric mass is Mbar_UV=sum_(A in U,B in V) M_AB and the degree is hbar_U=sum_(A in U)h_A. Then Mbar/hbar is the correct quotient Markov kernel. This does **not** identify it with the fundamental paired-black graph, which has a different line/loop convention. Switching has no analogous guard-preserving equivariance, as the n4 and n>=8 Session102 nonlift examples demonstrate.

## 4. Exact source-preserving recursion that does work

Append a final position n to x in F_n. For b in {0,1}^{n−1}, let E_n(x,b) preserve x, impose n−1→n, and use b for pairs (i,n), 0<=i<n−1. Then

    F_(n+1) = disjoint_union_(x in F_n,b) {E_n(x,b)},
    C_(n+1) E_n(x,b) = E_n(C_n x, complement(b)).       (6)

The second equation checks *every arc*: old nonpath arcs complement, the old path stays forward, the new consecutive arc stays forward, and the other new arcs complement. Thus it is an authenticated native intertwining, not an arbitrary product of matrices. It yields the exact adjacency recurrence

    M_(n+1)[A,B]
      = sum_(x in F_n) sum_b
          1{c(E_n(x,b))=A, c(E_n(C_n x,complement(b)))=B}. (7)

All source bits and the final port survive until the last class projection. Equation (7) builds the weighted operator and enumerates ordinary classes, but costs 2^beta at order n: no improved A000568 complexity is asserted. The number of nonempty source fibres is the class count, while their masses are the inherited H/Aut values. Reweighting source words by 1/h_(c(x)) would count one per class, but needs those global fibre masses and is not a closed scalar recurrence.

An exact SCC composition specialization is also available. For the canonical chain T=T_1→...→T_k, choose independently one normalized path word x_i from each child's fibre, concatenate the words, reverse all cross-child arcs except each exit-to-next-entry bridge, and complement child nonpath arcs. Call the resulting target G(x_1,...,x_k). Then

    M_[T],B = sum_(x_1,...,x_k) 1{c(G(x_1,...,x_k))=B},
    W_[T],B = (product_i a_(T_i)) M_[T],B.              (8)

Proof: SCCs occur contiguously in every Hamilton path, their automorphism group is the direct product, and the normalized global path words are exactly these tuples. The native arc rule gives G. This refines the inherited product of *row masses* to actual target weights, retaining the boundary ports and child words. It is not a product of unported child adjacency matrices: the old nonsingleton modules cease to be modules at bridge exceptions. For arbitrary substitutions with interleaved runs, (8) does not apply without retaining those runs and their joint ports.

## 5. Smallest class-pair obstruction to strip closure

Bits are backward arcs in lexicographic unordered-pair order. At n4 the pair order is (01,02,03,12,13,23), and C toggles mask22. Take x=4 and y=18=Cx. Their source and native-target classes all have canonical minimum code4. Each has eight strip extensions b in {0,1}^3. Their class-pair histograms differ. A single separating event is sufficient:

    #{b: (c(E_4(4,b)),c(C E_4(4,b)))=(10,40)} = 1,
    #{b: (c(E_4(18,b)),c(C E_4(18,b)))=(10,40)} = 0.    (9)

The first event occurs at b=4 (bits b_i occupy position i), giving the literal n5 source mask260 and native target mask106. All endpoint class codes in (9) use lexicographic n5 pairs and minimize over all vertex permutations. `checks.json` contains the full eight-entry histograms. At n2 and n3 every same-pair-class group has identical extension histogram; hence this is the smallest obstruction in this convention.

The problem is not only directionality or ordinary-vs-converse bookkeeping. Take n4 words x=2 and y=6. After merging each endpoint class with its converse and forgetting source/target order, both have the same pair key [2,4]. Their eight extension histograms still differ: folded endpoint pair [2,8] occurs once for x and never for y. Full receipts are included. Thus even that coarser metagraph cannot support a universal strip transition law depending only on its current endpoint-pair key.

Formally, the strip-extension counting kernel is not strongly lumpable for these partitions. This rules out a word-independent local transition assigned to those keys. It does **not** rule out every possible identity between total matrices at consecutive orders, nor a larger sufficient state, nor computation of averages after recovering the omitted fibre disintegration.

The pair-rigidity statement in §2 explains a genuine boundary: retaining the entire simultaneous-isomorphism type of the directed native pair retains the whole normalized word. Consequently there are exactly 2^beta simultaneous-isomorphism types of directed native pairs: any finite state that reconstructs such a pair must distinguish at least that many states. If the two tournaments are allowed to exchange, there are 2^(beta−1) types for n>=3, since C has no fixed words. This is a reconstruction-state lower bound, not a lower bound for computing class counts or aggregated adjacency. Independent endpoint isomorphism types forget exactly the compatibility that the obstruction detects. This is consistent with the earlier 81/27 decorated-cycle warning; no claim that scores, cyclic triples, or Hamilton counts recover automorphism transporters is used.

## 6. Exact bank and independent enumeration routes

|n|ordinary classes|total M mass|blue M mass|ordered support entries|diagonal M mass|
|---|---:|---:|---:|---:|---:|
|2|1|1|1|1|1|
|3|2|2|2|2|0|
|4|4|8|4|7|2|
|5|12|64|16|50|8|
|6|56|1024|64|732|16|

The checker builds a canonical dictionary from **every labelled tournament** by permutation orbits. It computes M and B by fixed-path words, and separately computes W and its blue part by all Hamilton permutations of each representative and actual XOR native targets. It checks (1), labelled orbit mass, exact Markov balance, blue invariance, and the full strip recurrence for n2→3 through n5→6. No floating-point arithmetic or external package is used. Runtime checks raise explicit exceptions, including under optimization.

    python3 check_weighted_transport.py > checks.json
    python3 -O check_weighted_transport.py > checks_optimized.json
    cmp checks.json checks_optimized.json

The matrix receipt gives canonical class code, a_A, H_A, h_A, and every nonzero M_AB and blue entry, so both ordinary and converse-merged operators are reconstructible. Native labelled connectivity, ordinary quotient connectivity, switching Cayley connectivity, and fundamental paired-black connectivity remain distinct questions. No all-order native connectivity result follows from this bank or from (7).

## 7. A compatible proposal clock and a gap-controlled mixed operator

The Session105 up/down audit gives a second, different authentic process K: choose a uniformly random vertex, delete it, and reinsert a fresh vertex with independent uniform incident orientations. On ordinary classes its stationary law is

    mu_A = (1/a_A) / sum_B(1/a_B),

not pi_A=h_A/2^beta from (3). Its self-adjoint spectrum is nonnegative and its nontrivial gap is 3/n for n>=3; see `../updown/UPDOWN_HARMONIC_AUDIT.md` for the proof and normalization. In this section K denotes the row-stochastic (backward) kernel, the transpose of the column-stochastic K in the up/down note; its spectrum, stationary law and reversibility are unchanged. Directly mixing K with (3) therefore generally changes the stationary law. This mismatch is already visible at n3: pi=(1/2,1/2) but mu=(3/4,1/4) in transitive,cycle order.

There is a native clock repair. At each tick choose one of **all n! vertex permutations** uniformly. If that word is a directed Hamilton path in the current tournament, perform its native move; otherwise stay. Do not condition on the guard succeeding. The ordinary class operator is

    R_AB = W_AB/n!                  if A != B,
    R_AA = 1 - H_A/n! + W_AA/n!,
    R = I + (W - diag(H))/n!.                         (10)

For each *fixed* word the guarded move, extended by identity on invalid sources, is an involution: it preserves all guard arcs when valid. Thus this proposal chain is doubly stochastic on labelled tournaments. Equivalently, (1) directly proves reversibility of R with mu. Native self-class moves remain in W_AA; it would be incorrect to subtract them a second time.

Reversal pairs the n! words without fixed points for n>=2, and at most one of each pair can be directed in a tournament. Therefore

    H_A <= n!/2,   R_AA >= 1/2.                       (11)

This is built-in laziness, not an added operation. Hence R=(I+Q)/2 for a reversible stochastic Q, so its entire spectrum is nonnegative. For any 0<alpha<=1, the actual mixed process

    L_alpha = alpha K + (1-alpha) R                   (12)

chooses the refresh process with probability alpha and the native proposal process otherwise. Both share mu and are positive semidefinite. On mean-zero functions, their Dirichlet forms satisfy

    <f,(I-L_alpha)f>_mu
      >= alpha <f,(I-K)f>_mu
      >= (3 alpha/n) <f,f>_mu.                        (13)

Thus **both the right and absolute spectral gaps of L_alpha are at least 3 alpha/n**, n>=3. There is an important clock-cost caveat: under mu, the native proposal acceptance probability is exactly 2^{−(n−1)}, since a fixed word requires n−1 independent forward arcs in a uniform labelled tournament. Thus this construction does not promise rapid native-only exploration; its proved gap comes from the refresh fraction. The refresh component makes the process irreducible (all labelled orientations can be installed by a finite sequence of vertex refreshes); laziness makes it aperiodic. This is a genuine gap-controlled metagraph exploration operator using native moves, not a proof that the native graph alone is connected and not a uniform sampler of isomorphism types. Its stationary class probabilities are proportional to labelled orbit sizes.

At n3 the native proposal matrix is exactly

    R = [[5/6,1/6], [1/2,1/2]],

whereas the conditioned-path operator (3) is the deterministic class swap. The checker verifies rational row sums, detailed balance, and half-laziness of (10) through n6 and this explicit n3 hostile. The new mixed-gap statement combines these elementary clock facts with the separately audited up/down theorem; no unverified spectral computation is a dependency.
