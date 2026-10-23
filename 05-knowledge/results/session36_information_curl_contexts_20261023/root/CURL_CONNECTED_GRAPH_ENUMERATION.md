# Curl classes, connected-graph numbers, and Hamiltonian-weighted logarithms

Status: PROVED, conditional only on the explicitly proved companion curl-fibre classification. This is an enumeration identity and a structural SCC quotient, not a graph-isomorphism-preserving identification.

## 1. Source, quotient and forgotten datum

For a fixed labelled vertex set, identify tournaments when their full signed triangle circulations agree. The companion lattice theorem proves that a fibre consists of the cyclic rotations of the ordered strong-component list. If there are k nonempty strong components, the fibre has exactly k elements. The labels remain fixed; rotations of disjoint nonempty labelled blocks are free.

This quotient forgets which gap in the circular list of strong components is the initial gap. It does not forget internal strong tournaments or their cyclic order. In particular every strong tournament is a singleton fibre. Reordering blocks arbitrarily need not preserve the curl.

## 2. Exact generating function

Let a_n=2^(n(n−1)/2), a_0=1, and let s_n be the number of strongly connected labelled tournaments (s_1=1). Use formal exponential generating functions

T(x)=Σ_(n≥0) a_n x^n/n!,   S(x)=Σ_(n≥1) s_n x^n/n!.

Unique ordered SCC decomposition gives T=1/(1−S). Let c_n count distinct labelled curl tensors. Passing from an ordered k-tuple of SCC blocks to its cyclic rotations divides its contribution by k, hence

C(x)=Σ_(n≥1)c_n x^n/n! = Σ_(k≥1) S(x)^k/k = −log(1−S(x)) = log T(x).

No positive analytic radius is claimed: T has radius zero. Every identity is coefficientwise in formal power series.

All labelled simple graphs also have total counts a_n. Their unique decomposition into a SET of connected components gives the same logarithm for their connected counts. Therefore c_n is EXACTLY the number of connected labelled simple graphs, for every n≥1, not merely a matching numerical prefix. The familiar connected-graph recurrence, also obtained by differentiating C=log T, becomes

c_n = a_n − Σ_(j=1)^(n−1) binom(n−1,j−1)c_j a_(n−j).

The initial values are 1,1,4,38,728,26704,1866256,251548592.

Primary exposition of the standard connected-graph generating function: University of Toronto MAT344 lecture notes, section 4, https://www.math.toronto.edu/balazse/2019_Fall_MAT344/Lec_10.pdf . The contribution here is its interpretation as full tournament-curl fibres, not a new discovery of connected-graph enumeration.

## 3. Weighted refinement preserves Hamiltonian counts

Let w(T) be any weight multiplicative across ordered strong components and unchanged by cyclic rotation. The same argument over a commutative Q-algebra gives

Σ_(curl classes [T]) w(T) x^|T|/|T|! = log(Σ_(labelled tournaments T) w(T) x^|T|/|T|!).

Hamiltonian path counts are multiplicative across SCCs: every directed Hamiltonian path visits the ordered SCC blocks consecutively, and the arcs between successive blocks allow every independent internal path choice. Consequently w(T)=H(T)^s is multiplicative for any fixed real s, as are the formal monoid weights [H] with [a][b]=[ab]. The latter retains the entire H-frequency distribution rather than only a moment.

For the first moment, double-count a tournament together with a Hamiltonian ordering. Each of the n! orders fixes n−1 arcs and leaves choose(n,2)−(n−1) arcs free. Thus for n≥1,

B_n=Σ_T H(T)=n!·2^((n−1)(n−2)/2), B_0=1.

If d_n is the sum of H over all distinct labelled curl classes, then

D(x)=Σ_(n≥1)d_n x^n/n! = log(1+Σ_(n≥1)2^((n−1)(n−2)/2)x^n),

d_n = B_n − Σ_(j=1)^(n−1) binom(n−1,j−1)d_j B_(n−j).

Equivalently d_n=Σ_T H(T)/k(T), where k(T) is the SCC count. The divisibility is explained by free rotation, not by numerical rounding. This D is a locally defined new series, NOT the user's earlier D(z) whose separate definition must be retained.

## 4. What this does not show

- Labelled counts do not determine unlabelled counts: automorphisms and periodic component necklaces require their own Burnside analysis.
- The equality with connected graphs is not an adjacency-preserving bijection, nor a proof that graph complementation represents path reversal.
- Although H is constant on each curl fibre, evaluating H from a reconstructed tournament may still require exponential work.
- This supplies neither the missing metagraph wiring nor the conjectured full Hamiltonian-count spectrum.

The structural message is precise: global source-to-sink placement is partly a gauge for Hamiltonian counting, while internal directed-cycle incidence is not.

## 5. Sharp naturality obstruction: equal counts, unequal actions

For every n≥3 there is no S_n-equivariant bijection between labelled tournament curl classes and connected labelled simple graphs. The connected complete graph K_n is fixed by every permutation. No tournament curl tensor is: fix a triple of distinct vertices, then its transposition changes the sign of that alternating triangle coordinate, which is always ±1 or ±3, never zero. An S_n-fixed tensor would require this coordinate to equal its negative. Thus the fixed-point counts already disagree for the full symmetric group.

At n=1,2 both sets are singletons and no obstruction exists. At n=3 the curl classes form two orbits of size2, whereas connected graphs form orbits of size1 and3. This gives an elementary exact boundary for the tempting 'trijection': enumerative agreement is proved, but a relabelling-natural bijection is impossible from n=3 onward. It does not forbid arbitrary bijections after choosing external label orderings.
