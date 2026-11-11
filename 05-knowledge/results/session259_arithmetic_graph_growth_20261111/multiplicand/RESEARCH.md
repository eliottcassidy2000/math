# Multiplicand growth: a squarefree persistence clock

Status: PROVED here (elementary arguments), with FINITE-EXACT experimental receipt. These are session research results, not independently audited canon promotions or claims of literature priority.

## Inheritance and concept board

Inherited: THM-362 (`01-canon/theorems/THM-362-natural-operation-graph-shadows.md`) identifies the ordinary multiplicand shadow with proper divisibility, and THM-2422 (`01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md`) separates fibres from shadows and gives strict/weak fibre counts. HYP-3003 (`05-knowledge/hypotheses/HYP-3003-summand-multiplicand-farey-basis-merge.md`) motivates retaining operation packets. No inherited theorem is re-audited here.

Live board: targetwise factor involution; unit-deleted divisibility order; prime-exponent box; prime-product simplicial complex; persistence under integer births; strict diagonal deletion. Hostiles: 4 and 9 (strict composite atoms), 12 (nonsquare yet nonsquarefree), 15 (first positive-dimensional topological birth), 30 (its death). The least-used sidecar is the order complex, distinct from the undirected graph's one-dimensional cycle space.

## 1. Exact one-node graph and order growth

Write D_N for the directed graph on 2,...,N with x→z iff x properly divides z. Its new vertex n receives τ(n)-2 edges. The corresponding Hasse graph receives ω(n) edges, except a prime n receives zero because its sole predecessor 1 was removed. In prime-exponent coordinates n↦(v_p(n)), divisibility is coordinatewise order; Hasse steps add one basis vector. Numeric cutoff is the weighted halfspace Σv_p log p≤log N, not an unweighted degree cutoff.

With the unit present the order complex is always a cone. Removing the unit is therefore structurally essential, not cosmetic.

For N≥4 the undirected D_N has exactly

    c(N)=1+π(N)-π(floor(N/2))

components: one core and the isolated primes in (N/2,N]. Every composite joins a prime divisor ≤N/2, and every prime p≤N/2 joins 2 via 2p (p=2 is already in the core). Hence its ordinary graph cycle rank is

    b1_graph(N)=Σ_(n=2)^N (τ(n)-2) -(N-1)+c(N).

This graph cycle rank is NOT the H_1 of the clique/order complex: comparable triangles are filled there.

## 2. Natural radical retraction

Let X_N be the order complex of the divisibility poset {2,...,N}. Let K_N be the simplicial complex whose vertices are primes ≤N and whose faces are finite prime sets S with product(S)≤N, including the empty face.

PROVED: X_N is naturally homotopy equivalent to K_N, compatibly with all inclusions N≤M.

Proof. The map r(n)=rad(n) is order preserving, idempotent, and r(n)|n. It retracts the poset onto its squarefree subposet. The pointwise inequality r≤id gives a simplicial homotopy (via the standard order homotopy/prism construction); r is independent of N. The squarefree subposet is exactly the nonempty face poset of K_N, whose order complex is its barycentric subdivision. Thus this is a natural filtered homotopy equivalence, not just a coincidence of Betti numbers. Exponents are destroyed; support and weighted birth times remain.

## 3. Wedge decomposition and the complete persistence barcode

For every N≥2 and integer k≥0,

    reduced b_k(X_N)=#{m: N/2<m≤N, m odd squarefree, ω(m)=k+1}.

More strongly, over any field (indeed integrally), the reduced persistence module has exactly one interval

    [m,2m) in degree ω(m)-1

for each odd squarefree integer m>1, and no others. Ordinary H_0 additionally has the essential interval [2,infinity).

Proof. Let A_N=star_{K_N}(2). This is a nonempty contractible subcomplex, nested with N. A face outside A_N corresponds precisely to an odd squarefree m satisfying N/2<m≤N. Every proper face of such a face has product at most m/3≤N/3<N/2, so lies in A_N. Consequently K_N/A_N is a wedge of spheres, one sphere of dimension ω(m)-1 per such m. This includes S^0 summands for isolated vertices. The relative chain complex C_*(K_N,A_N) has zero boundary and has those outside faces as its basis. Under inclusion each basis face remains itself until 2m enters, when it belongs to A_N and maps to zero. Relative homology agrees naturally with reduced homology because A_N is nonempty and contractible. This proves the interval decomposition with explicit maps, not merely dimensions.

An explicit integral cycle is the boundary of the simplex with vertices {2}∪{p:p|m}. The facet excluding 2 is born at m; every other facet has product 2m/p<m. The missing full simplex enters at 2m and kills this cycle. In degree zero this cycle is [p]-[2].

This produces a remarkable distinction: multiplication is diagonal-sensitive at squares, but its unit-deleted order topology is support-sensitive at squarefree numbers. A prime square changes the factor diagonal while changing no homotopy type. The nonsquare 12 is also topologically silent. The squarefree composite 15 births a loop; its double 30 kills it. All finite bars have the same multiplicative lifetime 2, or the same additive length log 2 in logarithmic time. No corresponding nontrivial barcode exists for the additive order shadow, whose order complex is always a simplex (even after deleting 1).

Generalization (PROVED by the same proof): restrict to any nonempty set P of allowed primes with least prime q. Replace 2 by q, odd by q-free, and intervals [m,2m) by [m,qm). Because every other allowed prime exceeds q, all proper facets of a q-free m are present before m. Thus the least allowed atom is a uniform multiplicative death clock.

## 4. Local birth links, Mobius signs, and squarefree events

At the birth of n in X_N its lower link is the order complex L_n of proper nonunit divisors of n. If n is squarefree with r=ω(n), this is the barycentric subdivision of the boundary of an (r-1)-simplex: L_n≃S^(r-2). For r=1 this uses the empty-complex S^(-1) convention.

If n is not squarefree, radical retracts its proper nonunit divisor poset to the squarefree divisors, which have maximum rad(n)<n; so L_n is contractible. Therefore every nonsquarefree birth is a homotopy equivalence. A squarefree birth contributes one relative homology generator in degree r-1. It births a class when odd, kills one in degree r-2 when even (except n=2, the initial component). The global barcode resolves what the local link alone cannot: whether the relative cell creates or kills.

The reduced Euler characteristic is

    reduced χ(X_N)=-Σ_(n=1)^N μ(n),  N≥2.

Proof: the nonempty faces of K_N correspond to squarefree n≥2, with sign (-1)^(ω(n)-1)=-μ(n). Subtracting 1 gives the formula. Equivalently the signed count of surviving odd squarefree m in (N/2,N] equals -M(N). This is an exact identity; no estimate for M(N), zero-free region, or new number-theoretic cancellation theorem follows from the topology alone.

## 5. Strict diagonals are a different graph, not the same poset

Deleting equal-parent factor witnesses removes precisely the edges a→a² from D_N (a≥2); hence exactly floor(sqrt(N))-1 edges are removed. This strict graph is generally not transitively closed. Its transitive closure is divisibility with precisely the relations p→p² removed for primes p. For composite a, choose a proper factor 1<d<a: the path a→ad→a² survives, since d<a and a/d<ad. All other divisibility edges were retained from the start. For a prime p there is no intermediate divisor between p and p².

Thus the radical retraction and barcode above apply to the weak ordinary divisibility shadow, NOT automatically to the strict graph or its reachability poset. In particular the strict reachability poset has minimal vertices primes and prime squares, whereas the weak one has only primes. This is the exact order-theoretic version of the strict-fibre atom defect.

A prime square p² is an isolated vertex of the strict undirected graph precisely while p²≤N<2p²: its only possible incoming proper divisor p has been removed, and its first multiple is 2p². At 2p² it joins p through the two edges p²→2p² and p→2p². This is a second doubling clock, distinct from squarefree homological death. Section 8 below classifies the complete strict-reachability order-complex persistence; applying the weak barcode alone to it already fails at N=4.

## 6. Connection contract and useful next questions

Source: labelled multiplicand filtration → weak divisibility order → radical support → prime-product complex. Map: forget second parent, then retain radical; subdivision connects posets and simplicial complexes. Preserved: divisibility homotopy and the entire filtered homology module. Lost: exponent multiplicity, factor counts, diagonal-square marker, and detailed ancestry. Needed sidecars: exponent vector and strict/weak witness flag. Cheapest decisive tests: 4 vs 6, 12 vs 15, and the interval [15,30).

Addition comparison: n adds n-1 shadow edges and its lower link is a simplex; multiplication adds τ(n)-1 with unit, or τ(n)-2 without, and its lower link sees squarefree support. The common swap involution detects doubles/squares in target fibres; the topology instead converts doubling into a death operator. These are two genuinely different uses of the same integer 2.

Niche continuation: strict reachability poset persistence is solved in Section 8, including higher prime-square corrections. Wildcard continuation remains open: operation-fibre enrichment that restores arbitrary exponent multiplicities to the radical barcode, rather than trying to infer them from topology that intentionally forgot those multiplicities.

## 7. Reproduction and scope

Run `python3 05-knowledge/results/session259_arithmetic_graph_growth_20261111/multiplicand/experiment.py`. The program independently reduces simplicial boundary matrices over F_2 through cutoff 300 and compares all finite and still-live intervals to the formula. It also checks graph and strict-reachability descriptions through 120. Those bounded checks are FINITE-EXACT corroboration only; the proofs above establish the all-order laws.

## 8. Extension: complete strict-reachability persistence (PROVED)

This section resolves the strict-reachability question left open in Sections 5–6. It concerns the ORDER COMPLEX OF THE TRANSITIVE CLOSURE of the strict graph, not the flag complex of its unrepaired edge set.

Let Y_N be that order complex. Its reduced persistent homology is the direct sum of:

1. The weak intervals [m,2m), degree ω(m)-1, for odd squarefree m>1.
2. One interval [p²,2p²), degree 0, for each prime p.
3. One interval [p²m,2p²m), degree ω(m), for every prime p and odd squarefree m>1.

These are indexed intervals with multiplicity. It is permissible that p divides m: p³ contributes a degree-one interval [p³,2p³). Thus diagonal deletion has higher interactions, not just extra isolated vertices. The first higher correction is [12,24) in degree one, from p=2,m=3. The next include [20,40) and [27,54). At p=2,m=15 it produces [60,120) in degree two.

### Proof by a common retract and relative chains

Remove all prime-square vertices from X_N and Y_N to obtain the same subcomplex Z_N: the two posets differ only in the relations p<p². Mapping every vertex p² to p and fixing all other vertices defines an order-preserving retraction onto this common subposet. For the weak poset this retraction is pointwise ≤identity, so Z_N≃X_N naturally. For the strict poset it remains an order-preserving retraction, although it is NOT pointwise comparable to identity at p², and therefore cannot be called a homotopy equivalence.

Distinct prime-square vertices are incomparable. In the strict poset, p² has no lower vertices. Its upper link L_(p,N) consists of p²k≤N with k≥2. Division by p² identifies this link with X_floor(N/p²). No prime-square vertex occurs in this upper link. Therefore Y_N is obtained from Z_N by attaching one cone p²*L_(p,N) for each p²≤N, and every cone intersection lies in Z_N.

The strict retraction splits the homology sequence of (Y_N,Z_N), naturally in N. At the relative-chain level the summands are the cone complexes modulo their bases. Consequently the extra homology is precisely the one-degree suspension of the REDUCED homology of each L_(p,N). When p²≤N<2p², the link is empty and the relative group is one copy of H_0, providing the additional isolated-component interval. Equivalently use augmented reduced H_(-1) of the empty complex, with lifetime [1,2) in the k=N/p² clock. Once the link becomes nonempty its H_(-1) vanishes. The previously proved weak barcode for its nonnegative reduced homology gives every interval in item 3 after scaling time by p² and shifting degree by one. Since the splitting and quotient maps commute with inclusions, this proves the persistent-module decomposition, not merely individual Betti numbers. The argument is integral and hence holds over every field.

### Consequences and boundary

Strict diagonal deletion retains the universal logarithmic lifetime log 2, while creating a tower of additional bars one degree higher at prime-square rescalings of weak births. This is an exact bridge between the square diagonal and squarefree topology. It does not retain every prime exponent: each correction remembers a distinguished squared prime and a squarefree odd support, not arbitrary multiplicity.

The simple proposal “weak barcode plus H0 intervals [p²,2p²)” is REFUTED by the first extra loop [12,24). Its repaired form is items 1–3 above. The formerly open strict-reachability persistence question is now resolved here; persistence of the raw strict graph's clique complex remains a different question.

Run `python3 05-knowledge/results/session259_arithmetic_graph_growth_20261111/multiplicand/strict_experiment.py`. Independent full order-complex boundary reduction over F2 at cutoffs 30,100,300 matches the entire predicted barcode, including all higher correction intervals. See `strict_receipt.json`. No weak barcode algorithm is used in constructing or reducing those strict chain matrices.


## 9. Strict Euler law and the single-event arithmetic word (PROVED)

Let M(t)=Σ_(n≤t) μ(n). For N≥2, the strict reduced Euler characteristic is

    reduced χ(Y_N) = -M(N) + Σ_(p prime, p²≤N) M(floor(N/p²)).

Indeed each prime-square relative summand contributes the negative reduced Euler characteristic of its weak upper link. For t=floor(N/p²)≥2 this equals M(t) by Section 4; for t=1 the empty link has reduced Euler characteristic -1, so the same contribution M(1)=1 holds. This is an exact identity, not a new cancellation bound for M.

The entire reduced barcode is a single-event word in integer time. Its birth labels are exactly:

- Odd squarefree m>1, in degree ω(m)-1.
- p²m with p prime and m odd squarefree, including m=1, in degree ω(m).

Every death label is twice its birth label. No two births coincide, no two deaths coincide, and no birth coincides with a death. Hence at each integer there is at most one reduced homology generator born OR killed.

For uniqueness, a correction birth has exactly one prime with exponent ≥2: that prime is p, with exponent 2 or 3. This reconstructs p and m and excludes ordinary squarefree births. If p=2 its 2-adic valuation is 2; if p is odd the birth is odd. Thus births have v_2 equal to 0 or 2, whereas deaths have v_2 equal to 1 or 3, proving the birth/death separation. The essential ordinary H_0 birth at 2 is excluded from this reduced event word.

Consequently every integer with v_2≥4 is homologically silent. An odd integer with two or more distinct repeated primes, or with any prime exponent ≥4, is also silent. Here “silent” means the induced map on every reduced homology group is an isomorphism, as follows from the complete barcode. It does NOT assert that the strict-complex inclusion is a homotopy equivalence.

The companion strict experiment checks this event uniqueness through cutoff 300 and the Euler identity at all integers 2≤N≤300 directly from independently computed strict bars and divisor-recursive Mobius values.
