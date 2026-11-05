# Four physical binary/ternary unit-clock families

Session 197. Constructive all-order families with a finite exact tensor certificate. This extends the previously studied primitive-amplitude-60 bridge by two bounded internal gaps, giving primitive amplitudes 60,100,156,260. It does not classify every unit clock or every cross-neighborhood word.

## 1. Explicit tournaments

Choose i,k in {0,1}, j>=0 and m>=0. Take three vertices a,u,b with transitive order

    a -> u -> b,    a -> b.

Take an ordered transitive chain I whose cross-neighborhood word is

    1 0^i 1 0^j 2 0^k 2.                             (1)

A chain vertex of type 1 loses to a and beats b,u. A vertex of type 2 loses to b and beats a,u. A vertex of type 0 beats all three. Earlier chain vertices beat later ones. Finally add a transitive module R of order m, with

    I -> R -> {a,u,b}.

There are no other unspecified arcs. For m=0 the module R is absent. Denote this actual tournament by T_(i,k)(j,m). Its order is 7+i+k+j+m.

Put

    B_(i,k)=4(3+2i)(5+8k) in {60,100,156,260},
    beta_j=-1 if j=0, and beta_j=1 if j>=1.

**Theorem.** For every permitted parameter,

    H(T_(i,k)(j,m))=B_(i,k) 2^j 3^m+beta_j.          (2)

The positive clone range m>=1 consists of strong tournaments. Indeed every chain vertex reaches R, R reaches all three outside vertices, a reaches the first type-1 chain vertex, b reaches a type-2 chain vertex, and u reaches b. Every vertex is therefore reachable from, and can reach, the R module.

The m=0 values in (2) are separately verified deletion values, not an unsupported evaluation of a positive-clone formula at zero.

## 2. Uniform run bounds and exact interpolation

For a rooted tournament, collapse each maximal run of transitive root clones in a Hamilton word. If there are t runs, the exterior fragments between them give t-1 vertex-disjoint directed paths from the root's out-neighbors O to its in-neighbors I. Each fragment contains an O-to-I edge. Hence t-1 is at most the maximum matching size of that bipartite crossing graph.

If a_t counts collapsed words, the exact clone count at positive multiplicity n is

    sum_t a_t t! S(n,t).

Indeed, the distinct clone labels are assigned to the t ordered nonempty runs, and each run has its unique transitive order. Expanding the Stirling terms expresses the sequence as a linear combination of q^n for q=1,...,t_max. Therefore a matching bound of two gives a sequence in the span of 1,2^n,3^n, determined by its first three positive values. A matching bound of one gives a sequence in the span of 1,2^n, determined by its first two positive values. These assertions concern positive multiplicities only.

For the R module, keep one representative r. Its out-neighbors are exactly a,u,b. Its crossing edges all meet {a,b}: a points to the two type-1 chain vertices, and b to the two type-2 vertices. These include a disjoint pair, so the matching rank is exactly two for every j. Thus, for fixed j, the positive-m sequence is in the span of 1,2^m,3^m.

For the middle 0^j module, keep one representative z, with any fixed m>=0. Its in-neighbors are exactly the chain segment before z. Among its out-neighbors, only a points into that earlier segment, and it points to both type-1 vertices. All later chain vertices lose to the earlier segment, every vertex of R loses to it, and b,u do not point to it. The crossing graph is consequently a nonempty star centered at a, uniformly in m. Thus, for each fixed m>=0, the positive-j sequence is in the span of 1,2^j.

For each of the four pairs (i,k), the verifier supplies the exact 3-by-2 tensor grid m=1,2,3 and j=1,2. First interpolate in m on each of the two j rows; then for each fixed m interpolate in j. This proves (2) for all j,m>=1. No joint-degree assumption or stabilization from a finite search is used.

For j=0, the same rank-two bound and three positive-m counts prove the negative-unit formula separately. For m=0, the same rank-one middle-module bound and the j=1,2 deletion counts prove the positive-j deletion formula separately. The single corner j=m=0 is counted directly. Thus all boundary cases in (2) are covered.

## 3. Exact finite certificate

The complete grid is m=0,1,2,3 and j=0,1,2 for each of four gap pairs: 48 actual graphs, of orders at most 14. Every graph is counted by both endpoint subset dynamic programming and literal enumeration of Hamilton vertex words. These methods agree with (2).

For reproducibility, the four primitive coefficients are B_00=60, B_10=100, B_01=156 and B_11=260. The grid entries are therefore B times 3^m minus 1 when j=0, B times 2*3^m plus 1 when j=1, and B times 4*3^m plus 1 when j=2. The receipt lists all 48 independently counted integers and all graph parameters, rather than trusting those descriptions as input counts.

The endpoint DP starts with D[{v},v]=1 and appends an unused out-neighbor. Summing the full-set endpoint counts gives H. The literal enumerator starts at each vertex and recursively chooses unused out-neighbors, incrementing only for complete words. All computations are integer exact. At order at most 14 even the crude bound 14! is below 2^64, so no count overflows the unsigned 64-bit implementation. Rank checks at every applicable grid source also verify the claimed two and one matching ranks.

Run the standalone `verify.cpp` under C++17. Optimized and unoptimized builds reproduce the same receipt. No graph bank, external source code, or symbolic package is required. The finite certificate plus the uniform run bounds is an all-order proof; the separate exploratory cross-word census is not a premise of this theorem.

## 4. Physical updates, padding, and exact root boundary

For j>=1 and m>=1, two genuine transitive-module cloning operations commute:

- Add one vertex to R: H -> 3H-2.
- Add one vertex to the middle zero module: H -> 2H-1.

For j=0, enlarging R instead gives H -> 3H+2. Inserting the first middle zero vertex changes the sign regime by H -> 2H+3; only subsequent middle clones use 2H-1. This distinction prevents a false claim that the binary update is unchanged at the empty-module boundary.

A leading block 0^ell in the cross word consists of universal-source padding and does not alter H. A trailing block 0^d consists of vertices that, together with R, form one transitive module of order m+d. Hence it multiplies the displayed ternary amplitude by 3^d. These are physical graph statements, not scalar substitutions without a graph.

For the source with R a singleton, let A=B_(i,k)2^j and beta=beta_j. The separate deletion case proves that its full root-run polynomial is

    R(z)=(z+1)[A(z+1)^2+beta].                        (3)

Indeed the positive clone sequence gives (a_1,a_2,a_3)=(3A+beta,3A,A), and the deletion count gives a_0=A+beta. Thus every source in this family has zero alternating boundary R(-1)=0. This is a theorem for this explicit family, not for arbitrary unit ternary sources.

Its negative amplitudes are B*3^d, always 4 modulo 8. Its positive amplitudes are B*2^j*3^d with j>=1, always 0 modulo 8. The apparent corresponding congruence for all physical unit ternary clocks remains unproved; these families supply examples and a mechanism, not that general theorem.

## 5. Actual prime and 7p witnesses

The receipt contains full adjacency rows, and both graph counters check these examples:

- (i,k,j,m)=(1,0,1,1): order 10 and H=601, prime.
- (0,1,1,1): order 10 and H=937, prime.
- (0,0,2,1): order 10 and H=721=7*103, with 103 prime.
- (0,0,1,3): order 11 and H=3241=7*463, with 463 prime.

Primality is certified by exact trial division through the square root. The maximal needed trial bounds are 24,30,10,21 for 601,937,103,463 respectively. The verifier checks every integer divisor up to those bounds, so it does not rely on a probable-prime test or an external table.

These are actual strong tournament witnesses with reusable module interfaces. No assertion of infinitely many primes or 7p values in any of these exponential families is made. The theorem also does not prove that (1), with its two bounded internal gaps, exhausts all cross words giving unit clocks. The broader finite search suggested the construction but has not established such an all-length classification.
