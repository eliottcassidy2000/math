# Fixed feedback strata: asymptotic class counts and forced wiring

Fix k >= 1, put n=m+k, and let m tend to infinity. All implicit constants below may depend on k. A presentation consists of an ordered transitive base B of size m and k individually labelled core vertices K, with arbitrary cross and core edges. Thus L=2^{km+binom(k,2)} presentations. Classes are ordinary isomorphism classes, not converse-merged classes. For finite identities assume 0<=k<=n; endpoint independence needs m>=2. The asymptotic claims only require fixed positive k and sufficiently large m.

## Uniqueness, rigidity, and ordinary classes

If a second transitive-deleting k-set exchanges r vertices of K with r vertices of B, its complement retains r core vertices and m-r base vertices. Each retained core vertex must have one of the m-r+1 insertion-cut profiles on that remaining ordered base. These are independent cross-edge conditions, so the probability of this necessary event is ((m-r+1)/2^{m-r})^r. Consequently

p(nonunique K) <= sum_{r=1}^{min(k,m)} binom(k,r) binom(m,r) ((m-r+1)/2^{m-r})^r = O_k(m^2 2^{-m}).

The conditions between retained core vertices can only decrease this probability. The r=1 contribution has order m^2 2^{-m}; fixed r>=2 terms are smaller asymptotically.

The chance that any two core attachment profiles coincide is at most binom(k,2)2^{-m}. Call presentations generic when K is unique and the profiles distinct. Any automorphism of a generic tournament preserves K, fixes the uniquely ordered transitive base pointwise, and fixes all core vertices by their profiles. Thus it is the identity. Each generic ordinary class has exactly k! presentations. Each nongeneric class has at least one nongeneric presentation, giving

N_{n,k} = L/k! (1+O_k(m^2 2^{-m})).

Here tau<=k, not tau=k: a larger transitive induced set always contains one of size m. In fact tau<k is automatically nongeneric for m>=1 because deleting one additional base vertex can replace a superfluous core deletion.

## Strong classes

For distinct core vertices u,v, each base vertex independently supplies u->b->v with probability 1/4. Therefore the probability that some ordered pair lacks such a path is at most k(k-1)(3/4)^m. Outside this exceptional event all core vertices lie in one SCC. The same conclusion suffices even when not all these two-step paths exist.

Conditional on the core lying in one SCC, the whole tournament is strong exactly when (i) the first base vertex fails to dominate all core vertices and (ii) the last base vertex fails to lose to all core vertices. Necessity follows because otherwise the first is a global source or the last a global sink. For sufficiency, the core SCC reaches the first base vertex, hence every base vertex; every base vertex reaches the last base vertex, which reaches the core SCC. Unconditionally, for m>=2 the two endpoint conditions are independent, each of probability 1-2^{-k}. Their intersection therefore has probability (1-2^{-k})^2. The strong event differs from that intersection only when the core is not in one SCC, so

|Pr(strong)-(1-2^{-k})^2| <= Pr(core not in one SCC) <= k(k-1)(3/4)^m.

No conditional independence is asserted. It follows that

N^strong_{n,k} = L/k! [(1-2^{-k})^2 + O_k((3/4)^m + m^2 2^{-m})].

For k=1 the core SCC condition is automatic; the endpoint probability is exact for presentations when m>=2.

## Self-converse and merged classes

For a generic self-converse tournament, every anti-isomorphism preserves K and reverses the base order. Its square is an automorphism, hence the identity. Its core permutation pi is an involution. There are only O_k(1) possible pi and core orientations. On cross edges the anti-isomorphism pairs positions (u,b_i) with (pi(u),b_{m+1-i}), forcing opposite edge indicators. Fixed positions are impossible; otherwise all orbits have size two and there are 2^{km/2} choices. Hence generic self-converse classes number O_k(2^{km/2}). Nongeneric classes number O_k(m^2 2^{(k-1)m}). Both are o(N_{n,k}); the same upper bound is o(N^strong_{n,k}). The converse-merged count equals (N+SC)/2, and analogously in the strong subfamily.

## Hamilton-path mass

Let f(T) be the number of k-subsets whose removal is transitive. On labelled choices of such a subset and its k! core labellings, Aut(T) acts freely: an automorphism fixing the core labelling also fixes the ordered base. Therefore a class T has k! f(T)/|Aut(T)| presentations, including nongeneric classes. For any isomorphism-invariant H,

(1/k!) sum_presentations H = sum_classes f(T) H(T)/|Aut(T)|.

For H equal to the number of directed Hamilton paths, H<=k!(k+1)^m: choose the core order and assign each base vertex to one of the k+1 intervening slots; its order within a slot is forced. Thus if V=sum_classes H/|Aut|, then

0 <= (1/k!)sum_presentations H - V <= (#nonunique presentations)(k+1)^m.

The [exact random-attachment transform](../session26-attachments/proof.md) proves the uniform-core expectation

E[H | core] = (k!/4^k)(k+1)^m + O_k(k^m),

we obtain

V = 2^{binom(k,2)-2k}[2^k(k+1)]^m [1+O_k((k/(k+1))^m + m^2 2^{-m})].

This is a weighted class mass, not an ordinary class count. No replacement of an individual H by its expectation is used: the expectation computes the complete presentation sum, followed by the controlled multiplicity correction.

Any blue-incidence mass bounded by SC times max H is o(V), with relative bound O_k(2^{-km/2}+m^2 2^{-m}). Combining the [proved all-edge lower transfer inequality](../session26-blue/RESULT.md) tau(T)+tau(kappa_P T)>=floor(n/3), the tau<=k stratum is independent for n>=6k+3. Under the established incidence convention boundary black edges therefore equal (V-B)/2~V/2.

## Upper transfer bound

Let S be a minimum feedback set in T, of size tau, and let B=V(T)\S. The induced T[B] is transitive. Under kappa_P, all edges reverse except consecutive arcs of the directed Hamilton path P. The exceptional edges with both endpoints in B form a disjoint union of paths. Delete a minimum vertex cover C of that path forest. Each path on ell vertices needs floor(ell/2) vertices, hence |C|<=floor((n-tau)/2). On B\C there are no exceptions, so the transformed tournament is the reverse of a transitive tournament. Therefore

tau(kappa_P T) <= tau + floor((n-tau)/2) = floor((n+tau)/2).

For tau(T)<=k, its neighbors lie in [floor(n/3)-k, floor((n+k)/2)]. This is a valid universal transfer bound, not a claim of sharpness.

## Exhaustive controls and scope

check.py exhaustively checks every labelled tournament on n=1,...,5, every directed Hamilton path, both transfer inequalities, every presentation for 1<=k<n, the presentation-multiplicity Hamilton mass identity, the uniqueness/profile/SCC probability bounds, and the conditional endpoint characterization whenever m>=2. Results are in results.json. All assertions passed. The expectation theorem is a separate dependency, not independently proved by this file. No classical-novelty or literature-search claim is made. These results concern fixed-low-feedback strata, not all tournament classes.
