# Switching turns path-forest edge degree into boundary degree

Status: PROVED-IN-NOTE (elementary identities and recurrence); FINITE-EXACT checks in companion script. Chan–Pak remains OPEN. No new theorem ID or literature priority claim.

## Inheritance and object distinctions

Closest inherited mechanism: THM-076 (`01-canon/theorems/THM-076-walsh-ocf-factorization.md`) and THM-259 (`01-canon/theorems/THM-259-walsh-degree-formula.md`) retain even-length path forests in the edge-Walsh expansion of H. We reprove the cancellation needed here and **push edge supports through their vertex boundary**, rather than claiming a new edge-Walsh formula.

The historical `07-reflections/the-cut-side-is-classical-clebsch-and-the-permutohedron.md` identifies K5 cut switching with the Clebsch graph: 16 switch states, five star generators. This is provenance, not a claim that its other cut/cycle prose is needed. A 16-state switch orbit is not the 16 entries of gl(4), nor automatically the requested 1-5-4-5-1 partition. Exact provenance of that particular partition is not established in this note.

Hostile: graph switching need not preserve H or tournament isomorphism class. A regular cyclic five-tournament has switching H histogram {1:5,9:5,11:5,15:1}. Its five vertex-star neighbours have H=9 although it has H=15. Least-used sidecar: the vertex boundary of an edge-Walsh forest.

Concept board: Hamiltonian paths / edge-Walsh forests / cut characters / folded-cube walks / scalar support versus labelled census. Anchor: a genuine recurrence for H observables; niche: a sharper Fourier carrier on switching orbits; wildcard: small eigenfunction cones as local constraints on attainable spectra.

## 1. Exact carrier and cancellation

Let T be a tournament on [n], n>=2, with a_ij=+1 for i->j and a_ji=-a_ij. Let s in {+1,-1}^n modulo global negation act by a_ij -> a_ij s_i s_j; write T^s and f_T(s)=H(T^s). There are 2^(n-1) states. For each permutation p,

f_T(s)=2^{-(n-1)} sum_p product_{i=1}^{n-1}(1+a_{p_i,p_{i+1}}s_{p_i}s_{p_{i+1}}).

Expand by selected edges of the permutation path. A nonempty selected component is a contiguous directed path. If one component has odd edge length, reverse its vertex block, choosing the first such block by its starting position. The selected edge positions and other components remain unchanged; its sign changes by -1. This is a fixed-point-free involution of expansion terms. Unselected boundary edges make no contribution, so they create no exception.

Survivors consist of r vertex-disjoint paths of positive even edge length, each using at least three vertices. The switching monomial cancels all internal s-coordinates and retains exactly the 2r endpoints. Thus:

**Boundary-degree theorem.** The switch Fourier expansion of f_T uses only even vertex sets S with |S|<=2 floor(n/3).

This is not the edge-degree statement 2 floor((n-1)/2): it is its boundary pushforward. Distinct edge forests with the same endpoint set aggregate with signs. The boundary quotient loses their internal vertices, lengths, edge signs and cancellation addresses; the aggregated Fourier coefficient is the necessary sidecar.

The constant term equals n!/2^(n-1). One proof: no nonempty path forest has empty boundary. A second proof: each ordered vertex permutation specifies exactly one switching state making all its n-1 path edges directed. Therefore

sum_s H(T^s)=n!.

## 2. A family of uniform finite-order recurrence engines

Let A be the switching adjacency operator, (Ag)(s)=sum_{v=1}^n g(s with s_v negated). The n=2 case counts its two coincident generators with multiplicity. On the even-set character chi_S(s)=product_{i in S}s_i, A has eigenvalue n-2|S|. Hence

product_{j=0}^{floor(n/3)}(A-(n-4j)I) f_T=0.

For W_k(T)=(A^k f_T)(identity), the same characteristic polynomial annihilates the sequence k>=0. W_k is the sum of H over all ordered length-k star-switch words, including repeated visits, NOT a count of distinct isomorphism classes. Its ordinary generating function has denominator dividing

D_n(z)=product_{j=0}^{floor(n/3)}(1-(n-4j)z).

An eigenvalue zero contributes only at k=0. The stated operator polynomial is **minimal uniformly over all labelled tournaments**, although an individual tournament can omit modes. To prove this, fix j<=floor(n/3), choose j disjoint three-vertex paths, and let S be their 2j endpoints. In the Fourier coefficient of chi_S, their edge monomial has coefficient, up to orientation convention, 2^j (n-2j)!/2^(n-1), which is nonzero: contract the j paths to blocks, arrange all n-2j blocks, and independently reverse the j even paths. Distinct edge monomials are linearly independent on the tournament sign cube. Thus this boundary Fourier coefficient is not identically zero, so the eigenvalue n-4j occurs for some T. Any common polynomial annihilator must vanish at every one of these distinct eigenvalues. This is uniform minimality, not simultaneous nonvanishing in every T. This is a recurrence in switch-word length for each fixed n, not an n-recurrence for A000568.

At n=5 the full Clebsch spectrum would allow 5,1,-3. H has **no -3 component**. Universally:

A f_T=f_T+30,  W_k(T)=H(T)-15/2+(15/2)5^k,

sum_{k>=0}W_k(T)z^k=(H(T)-15/2)/(1-z)+(15/2)/(1-5z).

Consequently each state's five neighbours have H-values summing to its own H plus 30. A local H assignment on Clebsch violating that law cannot arise from any five-tournament. At n=6 the allowed eigenvalues are 6,2,-2 and the recurrence is W_(k+3)=6W_(k+2)+4W_(k+1)-24W_k. At n=7 they are 7,3,-1.

## 3. What transfers toward the larger targets

Source: labelled tournaments and their switch states. Target: low-degree even Fourier functions on the folded cube. Map: T -> (s -> H(T^s)); predicate preserved: exact H at every switch and their walk sums. Lost: tournament identity, automorphisms and number of labelled tournaments representing one function. Sidecar for A000568: switching-class representatives with stabilizer/orbit data; these are not supplied by the scalar recurrence.

For Hamiltonian-spectrum reachability, the positive contribution is an explicit necessary linear constraint on every switch-generated candidate bank, together with exact average n!/2^(n-1). This supports constructing coupled H packets instead of guessing individual values. It does not prove that every positive odd eigenfunction or every locally consistent packet is realizable. Nor do walk moments imply interval support. Those remain OPEN inverse problems, directly relevant to using switching operations in constructive Chan–Pak searches.

Cheapest next substantive test: enumerate switch packets at n=6, retain the size-four boundary coefficients, and compare packet support under a chosen vertex-extension operation. Determine whether an endpoint-coefficient update closes after a bounded marked sidecar; no such closure is asserted here.

## 4. Reproduction and scope

Run `python3 05-knowledge/results/session250_sixteen_20261110/geometry/check_switching.py`.

The script enumerates every labelled tournament at n=2,3,4,5, and 14 explicitly seeded cases each at n=6,7 (transitive controls plus 12 pseudorandom bitmasks). Every switch is included. H uses exact integer subset dynamic programming. It checks the factorial sum, exact Fourier support bound, and the pointwise Clebsch identity for all n=5 cases. The permutation argument supplies an independent exact proof path. Higher-order random cases are finite tests, not exhaustive claims. The companion output records actual completion.
