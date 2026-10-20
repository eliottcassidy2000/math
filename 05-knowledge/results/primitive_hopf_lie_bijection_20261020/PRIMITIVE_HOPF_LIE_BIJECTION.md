# Connected sign-even graphs are the primitive Lie objects of tournaments

**Status: PROVED over Q**, relative only to the classical Fourier and free-Lie/PBW facts explicitly identified below: the normalized Walsh coalgebra isomorphism, its connected/primitive restriction, triangular primitive generators indexed by strong tournament types, and an all-order deterministic primitive-support type bijection. **PROVED:** the converse/edge-parity refinement and its integer recurrences; the paired-component correction to inverse-Euler enumeration with an involution. **FINITE-EXACT:** primitive matrices through order five, connected graph edge-parity census through order six, and integer enumeration through order sixty. No local graph construction, relabelling-equivariant set bijection, odd-prime preservation for the primitive matching, or major-conjecture consequence is claimed.

## 1. Inheritance, novelty, and the two even predicates

A graph X is **sign-even** when every automorphism reverses an even number of its edges, relative to any fixed orientation. This is not the condition that every vertex has even degree. The path P3 is sign-even; the triangle is not.

The important prior results are:

- [THM-4531, twisted Mallows–Sloane and natural-map obstructions](../../../01-canon/theorems/THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md).
- [THM-4721, restricted twisted Burnside, components and Walsh basis](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md), especially its components and Walsh sections, and its [full note](../even_graphs_20261019_sylow_refinement_components_walsh_duality.md).
- The earlier [October 11 compositional codec](../tournament_even_bijection_20261011c.md) already proves the component lemma, the weighted-Lyndon identity, and a complete type bijection using rank matching of connected atoms. Its [converse continuation](../tournament_converse_context_20261011d.md) already proves the SCC-palindrome formula. Those facts are inherited, not new here. In particular, the component lemma predates October 19.

The new bridge is that connected sign-even graphs are not merely equinumerous with Lyndon words: their signed orientation lines map **explicitly and linearly onto primitive tournament combinations**, and these primitives are the free Lie algebra on explicit corrections of strong tournaments. The resulting primitive coefficient matrix supplies a nonzero-support atom matching instead of arbitrary equal-rank pairing.

The live objects are: SCC words; graph components; the induced-subtournament coproduct; Walsh edge characters; primitive Lie brackets; and the converse involution. Their shared predicate is primitiveness. The main hostile is that a raw strong tournament is usually not primitive. Another is that normalized Walsh is not multiplicative for ordered join. Both are checked below.

The old endpoint-transfer packet [THM-266](../../../01-canon/theorems/THM-266-endpoint-transfer-boundary.md) concerns **degree-even** graphs, as its counts 1,2,3,7,16 show. Its row/column method is useful as a reminder to retain operation data, but it is not a sign-even theorem and is not imported into this proof. Likewise the June Glaisher/Witt reflection's rooted degree-parity completion is a labelled, different-vertex-count construction, not the present unlabelled correspondence.

## 2. Two coalgebras and one explicit Fourier map

Let H be the Q-vector space with one basis vector [T] per tournament isomorphism type, including the empty tournament 1. Grade it by vertex count. Define

- [T] [U] = [T → U], the ordered join, with all arcs from T to U;
- Δ[T] = sum over A subset V(T) of [T[A]] tensor [T[V(T) minus A]].

The product is associative. The coproduct is coassociative and cocommutative. Splitting the vertices of an ordered join independently in its two factors proves that Δ is multiplicative. Thus H is a connected graded cocommutative Hopf algebra. No automorphism denominators occur in its definition: all vertex subsets are counted, with their multiplicities.

For graphs, retain the orientation sign. On labelled graphs on [n], impose the relation

    [gX] = epsilon(g,X) [X],

where epsilon(g,X) is minus one to the number of graph edges whose upward orientation is reversed by g. An odd graph becomes zero; each sign-even graph type contributes a one-dimensional **orientation line**. A representative with upward-oriented edges chooses a basis of that line. Another orientation-parity choice can negate it, so there is no unmarked preferred positive generator.

Let C be this graded signed graph space. Define

    δ[X] = sum over A with no X-edge between A and its complement
             [X[A]] tensor [X[V minus A]].                         (2.1)

Induced orientations, rather than independently recanonicalized signs, are retained in this formula. This is exactly splitting the connected components. A graph is sign-even iff all its components are sign-even: orient all identical copies using the same model; permutations of copies reverse no edges, while internal component automorphisms have their own signs. Consequently C is the usual symmetric component coalgebra, and its primitive subspace is precisely the span of connected sign-even graph orientation lines.

For a graph X on [n], m = binom(n,2), put

    F(X) = 2^(-m) sum over all labelled tournaments T on [n]
                      phi_X(T) [T],
    phi_X(T) = product over ij in E(X), i<j, of t_ij,
    t_ij = +1 if i→j, and -1 otherwise.                            (2.2)

**Theorem 2.1.** F descends to a graded coalgebra isomorphism C → H. In particular,

    F(span of connected sign-even graph lines) = Prim(H).         (2.3)

**Proof.** On labelled objects, the characters phi_X form the Walsh basis. Relabelling acts by a signed permutation on the characters and by an ordinary permutation on tournaments; the Fourier transform intertwines these actions. It therefore descends to an isomorphism on their coinvariants. The graph coinvariants are exactly C, because a graph stabilizer with sign -1 kills its orientation line over Q; the tournament coinvariants are H.

For the coalgebra identity, fix a vertex split A,B. In ΔF(X), sum independently over orientations of edges internal to A, internal to B, and crossing A,B. If X has a crossing edge, its character is odd in that independent orientation bit, so the whole term cancels. If it has none, the crossing orientations contribute 2^(|A||B|), cancelling the corresponding part of 2^(-m). What remains is exactly F(X[A]) tensor F(X[B]). Summing over A proves ΔF = (F tensor F)δ. A coalgebra isomorphism preserves primitives. In characteristic zero the primitives of the symmetric component coalgebra are just its connected generators. ∎

This proof explains the mechanism: **connectedness is detected by cancellation across every nontrivial cut**, not by an equality of two counting sequences.

F is not an algebra homomorphism from graph disjoint union to tournament ordered join. For example, if a is the singleton and b is the directed triangle,

    F(three isolated vertices) = (3 a^3 + b)/4,

whereas F(singleton)^3 = a^3. This failure is essential, not a normalization bug.

## 3. Strong tournaments give triangular primitive generators

Every tournament has a unique ordered SCC word; its condensation is a transitive tournament. Thus, as an algebra, H is freely generated by the strong tournament types, counting the singleton as strong. Write S for their graded span.

The Eulerian idempotent is the finite-in-each-degree convolution logarithm

    e = log-star(id) = sum_{k>=1} (-1)^(k-1) (id-u epsilon)^star k / k.

Explicitly, for a nonempty tournament T,

    e(T) = sum over ordered partitions (A1,...,Ak) of V(T)
             (-1)^(k-1)/k · [T[A1] → ... → T[Ak]].                (3.1)

Every block is nonempty. This is a finite, effective rational formula.

**Lemma 3.1.** Every e(T) is primitive. For a strong type s,

    e(s) = s + a linear combination of ordered joins of smaller types. (3.2)

**Proof.** In the convolution algebra Hom(H,H tensor H), let alpha(x)=x tensor 1 and beta(x)=1 tensor x. Then alpha star beta = Δ and beta star alpha = Δ by cocommutativity. Since Δ is multiplicative, applying Δ to the convolution logarithm gives log-star(alpha star beta). Commutativity of alpha and beta splits the logarithm, giving e(x) tensor 1 + 1 tensor e(x). Connected grading makes all series finite on each homogeneous input. Formula (3.2) follows directly from (3.1): its k=1 term is s, and every k>=2 term is a nontrivial ordered join, whose factors have smaller vertex counts. ∎

**Theorem 3.2.** H is the free associative algebra on the primitive elements e(s), s strong. Its primitive Lie algebra, with bracket [x,y]=xy-yx, is the free Lie algebra on those elements.

The generator change is invertible recursively by vertex weight using (3.2). This identifies H, as a Hopf algebra, with the tensor algebra on S equipped with primitive generators. The assertion that its primitives over Q are the free Lie algebra is the classical free-Lie/PBW theorem. For a primary account, see [Reutenauer, Free Lie Algebras (1993)](https://doi.org/10.1093/oso/9780198536796.001.0001), or his [Handbook of Algebra chapter (2003)](https://doi.org/10.1016/S1570-7954(03)80075-X). We use this classical theorem, not a new PBW claim.

Hence (2.3) upgrades the inherited scalar equality to the explicit linear correspondence

    connected sign-even graph orientation lines
       --normalized Walsh--> free Lie algebra on e(strong types). (3.3)

A raw strong tournament does not itself satisfy the primitive condition. Already

    e(b) = b - a^3,
    F(P3) = -e(b)/4,                                             (3.4)

for the canonical-mask orientation of P3. The subtraction is the missing coordinate in the naive strong-versus-connected identification.

## 4. A primitive-support atom bijection, then a full type bijection

Order strong tournament types by (vertex count, least adjacency mask). Choose the usual Lyndon bracketing: for a Lyndon word L=UV, take V to be its longest proper Lyndon suffix and set P_L=[P_U,P_V]; a one-letter word s has P_s=e(s). The P_L form a homogeneous free-Lie basis.

For every weight n, expand

    F(X) = sum_{L of weight n} K_n(L,X) P_L                       (4.1)

over connected sign-even graph types X. Theorem 2.1 and Theorem 3.2 make K_n square and invertible over Q.

**Theorem 4.1.** The bipartite graph with vertices L and X and edge condition K_n(L,X) != 0 has a perfect matching at every weight n. Choosing the lexicographically first perfect matching gives an explicit all-order bijection rho_n from weighted Lyndon words to connected sign-even graph types, with a nonzero Lie/Fourier coefficient certifying each match.

**Proof.** The determinant of K_n is nonzero, so some permutation term is nonzero. Its entries form a perfect matching. Every finite matrix can be computed by enumerating graph and tournament masks, their relabellings, the ordered partitions in (3.1), and rational linear elimination in (4.1). An ordered search with ordinary bipartite-matching feasibility tests chooses the declared first matching. ∎

Now apply the inherited Chen–Fox–Lyndon factorization to the ordered SCC word of a tournament, giving its unique nonincreasing sequence of Lyndon factors. Replace each factor by rho_n(L) and take disjoint union. Conversely, decode each graph component by rho_n inverse, sort the factors nonincreasingly, concatenate and take ordered tournament join. The two maps are inverse and preserve vertex count. Graph component sizes are exactly Lyndon-factor weights.

This improves the old atom codec by a structural support condition. It still uses a globally selected matching and a declared alphabet order. It is not a local graph operation, is not relabelling-equivariant on labelled objects, and does not preserve automorphism orders or graph edge counts. The precise remaining primitive problem is to replace the support matching by a direct combinatorial rule producing one connected graph from one primitive word and explaining its nonzero coefficient without enumerating the whole degree. Equal dimensions do not provide that rule.

### Small certificates and the denominator wall

With pair order 01,02,... and least masks under relabelling:

- weight 3: the letter b is matched to P3, coefficient -1/4;
- weight 4: the word ab, bracket [a,e(b)], is matched to K1,3, coefficient 1/8; the unique strong letter c is matched to C4, coefficient 1/8;
- weight 5: aab is matched to K1,4, coefficient -11/384; ac is matched to connected graph mask 62, coefficient -1/64. The six strong letters are matched to masks 29,63,58,126,221,495 in strong-mask order 8,10,12,40,41,76.

The full matrices, generator expansions and matching are in primitive_walsh.json. In particular,

    e(strong mask 12 at n=5)
      = [12] - 2[4] - 2[9] - [0] + (2/3)[2] + (8/3)[5] + (2/3)[11].

All masks in this display are tournament masks on five vertices. Thus odd denominators occur already in the chosen primitive generator basis. The all-type Walsh Z[1/2] theorem does **not** imply that K_n is invertible over Z[1/2], nor any odd-prime automorphism-valuation refinement of this primitive matching. Rational nonsingularity suffices for Theorem 4.1. This is an explicit boundary, not an unresolved numerical issue.

There is nevertheless a stronger **determinant-only survivor**.

**Proposition 4.2.** With the stated unscaled graph lines and standard Lyndon brackets, det K_n is a signed integral power of 2 at every weight n, although individual entries need not lie in Z[1/2].

**Proof.** Full normalized Walsh is an isomorphism over R=Z[1/2]: its labelled inverse is again a Walsh transform, and passing to coinvariants kills the graph-orbit sign-2 torsion. Hence its degree-n determinant in graph-type and tournament-type bases is an R-unit. Replacing raw strong generators s by e(s) is unitriangular in SCC-word length, since each correction in (3.2) has at least two strong factors. Its determinant on each homogeneous free-associative piece is one. Passing from words in e(s) to ordered PBW products of standard Lyndon brackets is also unitriangular, by the leading-word form of Lyndon bracketing and unique Lyndon factorization. Thus the determinant of full F in the PBW basis is still an R-unit, even though the change of basis has rational entries.

The component-count filtration on C is its coradical filtration; PBW degree is the corresponding filtration on H. A coalgebra isomorphism preserves it. The associated graded map is the symmetric algebra map induced by F on primitives: apply the iterated reduced coproduct to see that its top k-fold term is the symmetric tensor of the k primitive images. Therefore the degree-n determinant factors as det K_n times a product of determinants of tensor and symmetric powers of K_j with j<n. Such determinants are integer powers of det K_j. Induction on n gives the assertion. ∎

The exact determinants at weights 1 through 5 are 1, 1 (empty matrix), -2^-2, 2^-6, and 2^-41. The proposition does not restore an integral primitive matrix: the 3-denominators above remain, so a product term of p-adic valuation zero may mix positive and negative entry valuations. No entrywise or automorphism p-local matching statement follows from this determinant alone.

## 5. Converse becomes edge parity on primitive Lie objects

Let theta send each strong type to its converse, and extend theta as a Lie automorphism of the free Lie algebra on e(s). Tournament converse reverses ordered joins. Formula (3.1), with block order reversed, shows that it commutes with e. Consequently, on Lie words of letter length r, actual tournament converse acts as

    eta = (-1)^(r-1) theta.                                      (5.1)

Meanwhile, reversing every tournament arc in (2.2) multiplies F(X) by (-1)^|E(X)|.

**Theorem 5.1.** The number a_n of connected sign-even graph types with an even number of edges is the +1 eigenspace dimension of eta on weight-n free Lie primitives. The number b_n with an odd number of edges is its -1 eigenspace dimension. In particular a_n+b_n=c_n, the usual connected sign-even count.

This is a linear eigenspace statement. It does not identify each graph with an individual Lyndon bracket eigenvector. From order six, theta nontrivially exchanges strong types, so letter-length parity alone is insufficient. At weight six there are 35 strong types but only 9 self-converse strong types; ignoring the other 13 converse pairs predicts the wrong parity split.

Let T(z) count all tournament types, R(z) count self-converse types, S(z) count strong types, and Q(z) count self-converse strong types, with constant terms T_0=R_0=1 and S_0=Q_0=0. The inherited SCC and palindrome laws are

    T(z)=1/(1-S(z)),      R(z)=(1+Q(z))/(1-S(z^2)).                (5.2)

Component decomposition and Theorem 5.1 give

    T(z)=product_n (1-z^n)^(-c_n),
    R(z)=product_n (1-z^n)^(-a_n) (1+z^n)^(-b_n).                (5.3)

The second identity also follows directly from the Walsh converse trace and component edge parity. It is a signed Euler product, not the ordinary Euler transform of connected self-converse objects.

For any integer series A(z)=1+sum A_n z^n, define its ghost sequence

    g_n(A)=n A_n - sum_{k=1}^{n-1} g_k(A) A_{n-k}.

Then g_n(A)=[z^n] z A'(z)/A(z), requiring only integer arithmetic. Put d_n=a_n-b_n. The efficient recurrence is

    n c_n = g_n(T) - sum_{k|n, k<n} k c_k,
    n d_n = g_n(R) - sum_{k|n, k<n}
                         k (d_k if n/k is odd, else c_k),
    a_n=(c_n+d_n)/2,      b_n=(c_n-d_n)/2.                       (5.4)

A separate Lie-character calculation provides a cross-check with extra structure. For vertex weight n and letter length r,

    dim L_(n,r) = (1/r) sum_{j|gcd(n,r)} mu(j) [z^(n/j)] S(z)^(r/j),
    tr(theta|L_(n,r)) = (1/r) sum_{j|gcd(n,r)} mu(j)
                           [z^(n/j)] A_j(z)^(r/j),              (5.5)

where A_j=Q for odd j and A_j=S for even j. Multiply the second line by (-1)^(r-1) and sum over r to obtain d_n. These are the character versions of Witt inversion: PBW gives the symmetric-power product for the tensor algebra, and logarithm plus Möbius inversion separates primitive powers. They are classical transform formulas applied to this specific alphabet and involution.

Initial new parity columns:

| n | c_n | a_n: even edges | b_n: odd edges | d_n |
|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 0 | 1 |
| 2 | 0 | 0 | 0 | 0 |
| 3 | 1 | 1 | 0 | 1 |
| 4 | 2 | 1 | 1 | 0 |
| 5 | 8 | 7 | 1 | 6 |
| 6 | 43 | 23 | 20 | 3 |
| 7 | 398 | 237 | 161 | 76 |
| 8 | 6413 | 3247 | 3166 | 81 |
| 9 | 184596 | 93584 | 91012 | 2572 |
| 10 | 9540998 | 4773466 | 4767532 | 5934 |

## 6. Recovered enumeration themes and a concrete historical correction

### 6.1 Integer Burnside buckets plus ghost inversion

The old [Burnside extension program](../burnside_enum_extensions.md) used factorial/LCD scaling, buckets of equal orbit exponents, and inverse Euler transforms. We reuse the mechanism rather than just its saved output. For a partition lambda let z_lambda=product d^(m_d) m_d! and c2(lambda)=sum floor(lambda_i/2)+sum_{i<j} gcd(lambda_i,lambda_j). Then

    T_n = sum over odd-part lambda partitions of n of 2^c2(lambda)/z_lambda.

For R_n, sum instead over partitions of n having all nontrivial parts congruent to 2 mod 4, and either no fixed point (n even) or one fixed point (n odd). These are exactly the permutations for which reversal followed by relabelling can fix a tournament. Every pair orbit either has even length without endpoint exchange, or odd length with endpoint exchange. Each allowed pair orbit gives two choices.

The script multiplies each term by n!, groups equal c2, performs integer bit shifts, and divides exactly once. It then applies (5.4). There is no floating-point rounding and no graph enumeration in the order-sixty run. Ordinary and signed Euler products are reconstructed independently as a check, while (5.5) recomputes the primitive dimensions and traces by letter length.

The older [closed-orbit/Witt reflection](../../../07-reflections/the-zeta-function-and-the-ocf-read-complementary-halves.md) supplies the second underused theme: log-derivative ghosts separate repetitions from primitives. Here the primitives are strong-type words and Lie elements, not closed walks in one tournament. No spectral or Hamiltonian-path claim is transferred between those different carriers.

### 6.2 Inverse Euler does not commute with taking involution-fixed multisets

This corrects a concrete historical enumeration label. Suppose connected types carry an involution, with c_n total connected types and q_n fixed connected types. A fixed multiset may contain a fixed component individually, or contain an equal number of each member of a nonfixed pair. Therefore its generating function is

    F(z)=product_n (1-z^n)^(-q_n)
                     (1-z^(2n))^(-(c_n-q_n)/2).                 (6.1)

If f_n=InvEuler(F)_n, then

    f_n = q_n + [2 divides n] (c_(n/2)-q_(n/2))/2.                (6.2)

Thus InvEuler of all self-converse oriented graphs counts **fixed connected components plus paired non-self-converse component atoms**. It does not count only connected self-converse oriented graphs.

The saved [connected_selfconverse_oriented.out](../connected_selfconverse_oriented.out) labels its inverse-Euler term at six as 727 connected self-converse oriented graphs. The correct connected value is **726**. The extra atom is the disjoint union of the out-star and in-star on three vertices. These are the unique nonfixed pair among the five connected oriented graph types on three vertices; three of those five are self-converse. This is the first possible paired-component defect.

The enumeration script independently recomputes ordinary and self-converse oriented-graph counts through six by summing over **every permutation**, with pair states absent/forward/backward. It obtains InvEuler fixed count 727 and corrected connected fixed count 726. This is not just subtraction from an unverified old table. The original inverse-Euler numbers remain valid for their repaired paired-atom interpretation; only their connected-object label is wrong.

## 7. Reproduction, checks and boundaries

Files in this packet:

- primitive_walsh.py: standard-library exhaustive graph/tournament atlas, rational Eulerian primitives, connected Walsh images, full coalgebra identities, Lyndon basis elimination, nonzero-support matching and parity controls;
- primitive_walsh.json: all matrices, matches and primitive generators through weight five;
- primitive_walsh.out: complete compact finite report;
- primitive_enumeration.py: integer Burnside buckets, unsigned/signed primitive recurrences, independent length-refined Lie-character transform, Euler-product reconstruction, and the oriented-graph involution correction;
- primitive_enumeration.csv: all nine enumeration columns through order sixty;
- primitive_enumeration_length.json: dimension and converse trace at each nonempty (vertex weight, letter length);
- primitive_enumeration.out: summary and validation count.

Run from this directory:

    python3 -B primitive_walsh.py
    python3 -B primitive_enumeration.py --max-n 60

Both scripts retain validation under Python -O because they use explicit check functions, not removable assert statements. The finite atlas bounds are declared; the all-order theorems do not rely on extrapolating these checks. The integer enumerator's documented default resource boundary is sixty.

The connection contract is now explicit. Source: connected sign-even graph orientation lines. Target: primitive Lie combinations of strong tournament types. Map: (2.2), followed by the triangular generator change (3.1). Preserved: vertex weight, primitiveness, converse/edge-parity eigenvalue. Lost by a set matching: linear coefficients and orientation-line signs. Necessary sidecars: the declared strong alphabet order, sign gauges if coefficients rather than support are retained, and the selected matching. Cheapest decisive tests: the nonprimitive raw triangle, the nonmultiplicative three-isolate image, the weight-five denominator 3, and the nontrivial strong converse pairs at weight six.

The principal open improvement is a direct local primitive atom rule. The present construction provides a rigorous linear bridge and a support-certified effective type bijection; it does not turn equal counts into a natural labelled correspondence, and it gives no LRC, plane-colouring, Collatz or other unrelated prize closure.
