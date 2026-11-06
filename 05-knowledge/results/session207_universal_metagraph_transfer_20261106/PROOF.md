# Universal tribonacci transfer of degree-labelled metagraph incidence

Status: PROVED using the inherited odd-cycle collection formula. This generalizes the session's H=3 incidence calculation to every fixed core and every Hamilton path. It is not a new global count-spectrum completeness theorem.

## 1. Arbitrary pointed-core tail theorem

Let B be a nonempty tournament with a distinguished vertex q. Let R_a be the tournament on 1,...,a whose consecutive arcs point forward and whose longer arcs point backward. Define L_a(B,q) by placing R_a before B, orienting every cross arc from B toward R_a, except the single arc a->q when a>=1. Put f_a=H(L_a(B,q)), f_0=H(B), and use H(empty)=1 in residual counts.

For each directed simple path Q in B starting at q, including the singleton path, assign weight H(B-V(Q)). Let W_even and W_odd be the sums of these weights over paths with respectively even and odd numbers of vertices. Then

    sum_(a>=0) f_a x^a
      =[H(B)(1-x^2)+2x W_even+2x^2 W_odd]/D(x),
    D(x)=1-x-x^2-x^3.

In particular

    f_0=H(B),
    f_1=H(B)+2W_even,
    f_2=H(B)+2W_even+2W_odd,
    f_a=f_(a-1)+f_(a-2)+f_(a-3) for every a>=3.

### Proof by complete crossing-cycle decomposition

Use the inherited odd-cycle collection formula H(T)=sum_C 2^{|C|}, over collections of vertex-disjoint directed odd cycles of length at least three. Let h_i=H(R_i), h_0=1. As in the known tribonacci partner calculation, every cycle of R_i is a consecutive interval, so

    sum_(i>=0) h_i x^i=(1-x^2)/D(x).

If an odd cycle meets both the tail and B, its smallest tail vertex i forces the entire suffix i,i+1,...,a and then the exceptional arc a->q: no other arc exits the tail into B. Its remaining vertices form a directed simple path Q in B starting at q, and its final vertex returns to i. Conversely this construction closes whenever t+|Q| is odd, where t=a-i+1>=1. The only exceptional failure of the closing arc would be t=1 and Q=(q), but this has even total size2 and is excluded automatically. Thus every allowed odd total is at least3 and closes legally.

Any two crossing cycles share a and q, so at most one belongs to a disjoint collection. With no crossing cycle the contribution is h_a H(B). With one crossing cycle using t tail vertices and Q, the remaining cycles lie in the prefix R_(a-t) and the induced core B-V(Q), independently; its weight is2h_(a-t)H(B-V(Q)). There is no remaining eligible crossing cycle because a and q were removed. Consequently

    f_a=h_a H(B)
       +2 sum_(1<=t<=a,t odd) h_(a-t)W_even
       +2 sum_(2<=t<=a,t even) h_(a-t)W_odd.

Summing over a gives the theorem. No Hamilton-path insertion assumption or scalar closure approximation is used.

Since D is an irreducible cubic over Q (neither1 nor-1 is a root), its nonzero numerator of degree at most2 cannot cancel any factor. Therefore every fixed pointed core has minimal constant-coefficient linear recurrence order exactly3 over Q, not merely at most3.

## 2. Two-sided transport for every actual metagraph edge

Let T have a directed Hamilton path P=(v1,...,vm), m>=1. Put U=K_P(T), where K preserves the path arcs and reverses all others. For a,b>=0 define

    T_ab=Tr_a -> T -> Tr_b,
    P_ab=(left forward chain,P,right forward chain),
    U_ab=K_(P_ab)(T_ab),
    g_P(a,b)=H(U_ab).

The source always satisfies H(T_ab)=H(T). The induced blocks of U_ab are R_a,U,R_b. All cross arcs point from right to left except last(left)->v1 and vm->first(right). This is the same actual path-complement operation used in the repository, not a graph whose edges add vertices.

Fix b and regard the middle-plus-right block as an arbitrary pointed core with marked v1. The tail theorem applies to a. Fix a and apply its converse version to the right tail with marked vm. Hence

    g_P(a,b)=g_P(a-1,b)+g_P(a-2,b)+g_P(a-3,b), a>=3,
    g_P(a,b)=g_P(a,b-1)+g_P(a,b-2)+g_P(a,b-3), b>=3.

The entire infinite array is determined by the nine specified corner values g_P(i,j), 0<=i,j<=2; this is a sufficient signature, not a claim that nine is minimal information for every family. Its bivariate generating function is

    sum_(a,b>=0)g_P(a,b)x^a y^b=N_P(x,y)/(D(x)D(y)),

where N_P has degree at most2 in each variable. For 0<=i,j<=2 its coefficient is

    [x^i y^j]N_P
      =sum_(0<=r<=i,0<=s<=j) d_r d_s g_P(i-r,j-s),
    d_0=1, d_1=d_2=-1.

The numerator need not have nonnegative coefficients. The formula is a proof-based finite reduction: compute nine ordinary Hamilton counts of tournaments of orders at most m+4, then generate every padded neighbour count at every larger order. It does not enumerate all cores or identify distinct neighbour isomorphism classes.

In particular the H=3 packet is one explicit specialization; its nine counts are not an accidental tribonacci fit.

## 3. Actual fibre degrees stabilize after three vertices on both sides

The source automorphism group is always Aut(T_ab) isomorphic to Aut(T), by its unique ordered SCC decomposition. Its self-converse status is unchanged under equal padding a=b, but not asserted unchanged for unequal padding.

For a,b>=3 there is an exact target automorphism identity

    Aut(U_ab) isomorphic to Aut(U)_{v1,vm},

where the right side fixes the two endpoints pointwise (one vertex if they coincide). Indeed the left tail scores are 1,1,2,...,a-2,a-1. A middle vertex v has score

    d_U(v)+a-1[v=v1]+1[v=vm],

which is at least a-1>=2. Right-tail vertices have score at least a+m>=4. Thus the only outdegree-one pair is the first two left-tail vertices; its directed internal edge fixes them individually under every automorphism. The second vertex's unique outneighbour pins the third. Successively, the i-th tail vertex has exactly one not-yet-fixed outneighbour, the next vertex; the final vertex pins v1. This fixes the entire left tail. The dual indegree-one argument fixes the right tail and vm. Restriction gives an automorphism of U fixing the endpoints; conversely every such automorphism extends by fixing the tails. This proves the group identity.

Therefore each padded neighbour's ordinary fibre degree, for a,b>=3, is exactly

    g_P(a,b)/|Aut(U)_{v1,vm}|.

The denominator is independent of a,b and can be computed once from the original marked neighbour. The same tribonacci recurrences hold for these fibre degrees within the region where all terms have stabilized; for the a recurrence this means a>=6,b>=3, and dually for b. They need not hold across the small-padding boundary where automorphism orders can change.

## 4. Exact incidence measures along every fixed-core padding family

Every Hamilton path of T_ab uniquely extends a Hamilton path P of T. Thus for any integer k and every nondegenerate padded order m+a+b>=3, the number of labelled neighbours of a fixed labelled T_ab with count k is exactly

    #{P in HP(T):g_P(a,b)=k}.

Different paths give distinct labelled neighbours because the common arcs recover the preserved path. If m+a+b=1 or2, the complement operation is the identity, so this distinct-neighbour statement is not applied.

For original class-line incidence use one representative of each Aut(T)-orbit of Hamilton paths. Every orbit is free, and each determines one incidence flag at the source. Then the same count over orbit representatives gives the exact incidence multiplicity toward count k; a loop contributes two source flags, as usual. This is a pushforward incidence measure, not identification of the target class by its count. It remains valid when multiple marked paths lead to isomorphic targets.

For a,b>=3, replacing each g_P by g_P/|Aut(K_P(T))_{v1,vm}| gives the exact degree-labelled incidence measure of ordinary class fibres. The source fibre degree H(T)/|Aut(T)| is fixed. Thus the finite set of path-orbit signatures (nine counts and one endpoint-stabilizer order per orbit) determines every degree-labelled incidence row in this infinite two-sided padding family.

The aggregate sum over path-orbit neighbours of their ordinary counts also obeys both tribonacci recurrences. No such linear recurrence is claimed for the number of distinct target counts, which involves collisions of these sequences.

## 5. Consecutive-size adjacency rays, not only subsequences of orders

The same degree-pinning argument with no right tail gives, for a>=3,

    Aut(U_a0) isomorphic to Aut(U)_{v1}.

There is no need to fix the old last endpoint in this one-sided statement. All core vertices still have score at least a-1, and the unique left degree-one pair forces the tail and v1 pointwise.

Therefore every attained source count h has a realizing core T of order m and a selected path P such that at every order m+a, a>=3, there is an actual class line with source fibre H(T)/|Aut(T)| and target fibre

    f_a/|Aut(K_P(T))_{v1}|,

where (f_a) is a single exact tribonacci sequence. The source's ordinary count remains h at every such order. Its class-fibre degree remains fixed, even though unequal padding need not preserve converse status.

The tail formula gives f_0<=f_1<=f_2; the recurrence preserves monotonicity and yields f_(a+3)>=3f_a for every a>=0. Hence

    f_a>=3^floor(a/3) H(U).

Thus these neighbour degrees grow without bound along all sufficiently large consecutive sizes, rather than only an arithmetic subsequence. More precisely, the primitive tribonacci companion matrix and its Perron root lambda (the real root greater than1 of lambda^3=lambda^2+lambda+1) give f_a~C lambda^a with C>0. This is a standard finite-dimensional Perron-Frobenius consequence of the positive seed vector, not a heuristic estimate. It concerns the constructed family only, not all neighbours or all tournaments at that size.

## 6. Scope, dependencies and controls

This proves a universal exact transport for padded fixed cores, not an all-core classification of the metagraph. In particular it does not prove the global odd-except7-and21 count conjecture, general connectivity, or equitable lumping by count; the H3 packet already gives a fixed-size obstruction to the latter.

Inherited dependencies: THM-002-ocf.md; THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md; the repository's actual path-complement and free Hamilton-path orbit conventions. This packet reuses the previously audited score-pinning idea from session207_metagraph_fibre_amplification_20261106, but its general tail recurrence and complete marked-core transport are additional deductions. No external novelty claim is made.
