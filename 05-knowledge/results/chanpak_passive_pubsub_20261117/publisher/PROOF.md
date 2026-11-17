# A publisher cannot independently program its two-cover count at fixed H

**Status: PROVED for every finite tournament.** The ingredients are the classical SCC decomposition and the elementary Hamiltonian-cycle property of strong tournaments. The application is to the proposed variable publisher interface (H,P₂); this does not assert a new proof of the classical ingredients.

Write H(T) for Hamiltonian paths, P₂(T) for ordered two-path covers with empty paths allowed. Use H(empty)=P₂(empty)=1.

## 1. A finite core for every fixed path count

Every tournament is the unique ordinal sequence of its strongly connected components. Its H and P₂ are products of the respective component quantities: a path traverses the components in order, and each of two cover paths does likewise, possibly skipping components. In particular a singleton contributes (H,P₂)=(1,2).

Delete all singleton SCCs, retaining the order of the others. Call the result C the **singleton-free publisher core**. If s vertices were deleted,

    H(T)=H(C), P₂(T)=2^s P₂(C).

Every nontrivial strong tournament on r vertices has a Hamiltonian cycle, so H>=r>=3, by taking its r distinct cuts. For completeness, take a longest directed cycle. Any outside vertex with both in- and out-neighbors on it can be inserted. Otherwise all outside vertices divide into those dominating the entire cycle and those dominated by it. Strongness forces an edge from the latter set to the former, which inserts two vertices and contradicts maximality. The cases where either set is empty contradict strongness. Thus the cycle spans.

Let h_i be the H values of the nontrivial components and r_i their orders. Then r_i<=h_i and h_i>=3. If there are at least two components, sum h_i<=product h_i; the inequality follows inductively from a+b<=ab for a,b>=3. With one component it is immediate. Consequently, for h=H(T)>1,

    |C|=sum r_i <= sum h_i <= product h_i = h.

For h=1 the core is empty and T is transitive. This bound is independent of the number of singleton SCCs and hence of the total order of T.

## 2. Exact classification of the P₂ fiber at fixed h

Let C_h be the finite set of isomorphism classes of tournaments C having no singleton SCC, at most h vertices, and H(C)=h. Then, allowing all tournament orders,

    {P₂(T): H(T)=h} = union_(C in C_h) {2^s P₂(C): s>=0}.

The converse is constructive: append s singleton SCCs to C. Thus this is equality, not merely a bound. For h=1 use the empty core and restrict to s>=1 if empty tournaments are excluded.

Equivalently, there is a finite dictionary mapping each attainable odd part p of P₂ to a smallest exponent e_h(p), and

    P₂ is attainable at H=h iff oddpart(P₂)=p in the dictionary
    and v₂(P₂)>=e_h(p).

So **unbounded independent programming of (H,P₂) is mathematically impossible**: at any fixed h, only finitely many odd parts of P₂ are allowed. The number of distinct two-cover values up to X is O_h(log X). This is not limited to a recursive construction, a finite set of seeds, or positive glues.

At a prescribed order n>=h, the normalized interface stabilizes exactly:

    {P₂(T)/2^n : |T|=n, H(T)=h}
        ={P₂(C)/2^|C| : C in C_h}.

For n<h retain only cores with |C|<=n. The positivity condition is precisely the available singleton padding. As a simple consequence, v₂(P₂(T))>=max(0,n-h) whenever h>1.

## 3. Target-first publisher inversion and the necessary bypass

For a fixed exterior reader N=alpha*P₂+beta*H with alpha>0, fixing H=h confines N to finitely many geometric rays

    N=alpha*P₂(C)*2^s+beta*h.

Given a candidate core C, a prescribed N selects s exactly: (N-beta*h)/(alpha*P₂(C)) must be a nonnegative integer power of2. At prescribed order n it must be exactly2^(n-|C|). No long search over padding or singleton positions is needed.

A constructive bypass must therefore let the **nontrivial core** vary with the target, and normally let h vary as well. Adding singleton states cannot independently supply missing P₂ odd parts. Two candidate strategies remain meaningful:

1. Choose h from the target's residual congruences first, then construct a core with the required odd part and valuation. This is a target-owned core problem, not independent coordinate programming.
2. Keep a common unmarked publisher core but expose rooted or marked response coordinates. Such coordinates are not determined by (H,P₂); earlier neutral rooted pairs are concrete evidence. An exterior capable of reading these coordinates can avoid the fixed two-coordinate fiber obstruction. No universal realization theorem for this richer interface is claimed.

The theorem concerns only the unmarked (H,P₂) interface. It does not rule out a sufficiently rich marked publisher with a target-dependent H.

## 4. Exact eventual multiplicity, not only a finite census

There is also a complete finite-core formula for how often a fixed h appears by tournament order. For each core C in C_h, let r(C) be its number of nontrivial SCCs and R(C)=|C|. Its s singleton vertices can be placed in the r(C)+1 gaps before, between and after the nontrivial components. The SCC sequence is canonical, so distinct gap lengths give nonisomorphic tournaments, and every tournament arises exactly once. Thus

    M_h(n)=sum_(C in C_h, R(C)<=n) binom(n-R(C)+r(C), r(C)),
    sum_(n>=0) M_h(n) x^n
        =sum_(C in C_h) x^R(C)/(1-x)^(r(C)+1).

For n>=h, M_h(n) is exactly a polynomial in n, of degree max_C r(C), at most floor(log_3 h). The coefficients depend on the finite but potentially expensive core census. This is not an efficient enumeration algorithm or a proof that C_h is nonempty for every odd h except7,21.

For an attained prime h, every core is strong, so the polynomial is linear:

    M_h(n)=A_h(n+1)-B_h,  n>=h,
    A_h=#C_h, B_h=sum_C |C|.

For h=1 there is exactly one transitive tournament at every positive order. One may refine the same generating function by the normalized P₂ value P₂(C)/2^R(C); the finite interface dictionary then carries exact multiplicity data as well.

The mechanism explains a useful split: fixed-h large-order growth eventually consists entirely of singleton placement, while genuinely new publisher response types must already occur in bounded nontrivial cores. It does not classify those cores uniformly as h grows.

## 5. The entire unmarked slot response stabilizes simultaneously

The same argument is not limited to two paths. For k labeled path slots, including empties, P_k multiplies across ordinal SCCs and a singleton contributes k. Hence

    P_k(T)=k^s P_k(C) for every positive integer k,
    P_k(T)/k^n=P_k(C)/k^R.

The function P_z is a polynomial in z: partition the vertices into j nonempty directed paths, then inject these paths into z labeled slots, giving a multiple of the falling factorial(z)_j. Thus the displayed equality holds as a polynomial identity before normalization, and as a Laurent-polynomial identity after it.

For fixed h, **the full normalized unmarked slot polynomial P_z(T)/z^n belongs to a finite menu**, independent of n once n>=h. Measuring more unmarked cover orders does not manufacture new normalized publisher types by increasing singleton padding. Conversely this does not classify rooted response polynomials: singleton positions and marked endpoints can remain visible to rooted subscribers even when all the unmarked P_k agree. That is precisely why marked context, rather than another scalar unmarked statistic alone, is a plausible bypass.

### Replay

`python pubsub/publisher/verify.py` exhaustively checks all1099 labeled tournaments of orders1 through5. Canonicalization gives the expected1,1,2,4,12 isomorphism classes. It checks core deletion and H/P₂ factorization for every graph, and the exact joint(H,P₂) multiplicity formula for every class at these orders. The all-order statements above are proved algebraically, not inferred from this finite replay.
