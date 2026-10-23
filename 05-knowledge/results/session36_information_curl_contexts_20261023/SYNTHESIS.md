# Session36: what the tournament remembers

Status: audited proofs and exact finite certificates, with imported theorems and open boundaries identified below. Session36 revisits older work rather than treating all recovered facts as new. This synthesis is an account of specific maps and their losses, not a universal theory of mathematics.

## 1. The central finding

A tournament can carry several different kinds of information at once:

1. its internal strongly connected pieces;
2. their linear order;
3. the cyclic order left when an initial gap is forgotten;
4. the labelled triples which are cyclic;
5. the run counts visible when the whole tournament is used as one module;
6. the much smaller scalar Hamiltonian count.

These are NOT interchangeable invariants. Some support composition; some only support closed observations. The session produced exact equivalence theorems, ranks, generating functions and counterexamples specifying the distinctions.

The most useful general research move is therefore: choose the allowed interface and operations first, then determine the smallest state which predicts every permitted observation. A rich invariant can still fail this test if the next operation needs a coordinate it discarded.

## 2. Full signed curl: a circular SCC word

[Metric curl recovery](lattice/METRIC_CURL_RECOVERY.md) proves that the standard integral augmentation lattice, its Euclidean metric, and its tournament alternating form together classify tournaments up to cyclic rotation of their ordered strong-component lists. For strong tournaments this is a complete isomorphism invariant. For all tournaments it determines H, and in fact the complete multiplicative path-colour/path-cover profile.

This is full signed triangle circulation, not just the set or number of cyclic triangles. A transitive triangle still has circulation±1; a cyclic triangle has±3.

[The inverse decoder](root/CURL_REALIZABILITY.md) gives necessary and sufficient sign constraints for an arbitrary proposed rooted curl tensor and an O(n³) decoder. There are at most n labelled realizations, exactly one per SCC rotation. Local allowed triangle magnitudes and tetrahedron identities alone fail already at four vertices: the sign potentials must be compatible.

Neither the alternating Smith type alone nor the real metric-operator spectrum alone retains enough information. The integral placement relative to the root lattice is essential.

## 3. Two different Hamiltonian interfaces

[Module contexts](context/RESULT.md): at known order n, the exact linear connection rank for one-hole uniform module substitution is ceil(n/2). The complete path-cover profile is necessary and sufficient for all those responses. Known order and ceil(n/2)−1 actual cyclic-context counts recover it; the unmeasured leading coefficient is fixed. The path-colour polynomial, endpoint jet and full signed linear-forest polynomial are inherited invertible versions of the same data.

[Vertexwise contexts](context/CURL-AND-CLAN-COMPLETION.md): when the tournament instead supplies the individually addressable quotient slots, all nonempty substitution responses are equivalent to the full multivariate walk series, and exactly to its labelled cyclic-triangle hypergraph, or absolute triangle curl. A positive full-support walk coefficient detects each cyclic triple.

The final equivalence uses the CITED theorem of Boussaïri–Ille–Lopez–Thomassé: equal labelled cyclic-triangle structures are related by homogeneous-interval reversals. The response interpretation and explicit high-content detector are proved in the packet. Indecomposable tournaments are recovered up to converse; merely strong tournaments are not. An explicit strong seven-vertex pair is nonisomorphic and nonconverse but has identical complete vertexwise responses.

[Cyclic-triple realization](context/REALIZATION_PHASES.md) also has an exact positive-parity/negative-NAE recognition model. Positive consistency alone first fails at order5; the residual clauses cannot be discarded. This is a sound and complete finite constraint reduction, not a polynomial-time recognition or enumeration claim.

These theorems concern exact interfaces, not every possible attachment and not a general fast H algorithm. One-hole linear rank does not apply if the context duplicates the varying input.

## 4. Enumerative bridges which survive proof

[Labelled curl counts](root/CURL_CONNECTED_GRAPH_ENUMERATION.md) have exponential generating function log T(x), where T counts all labelled tournaments. Consequently they equal the connected-labelled-simple-graph counts for every order. The proof uses free cyclic rotation of disjoint labelled SCC blocks, not a numerical prefix. A weighted logarithm gives the total H and the full multiplicity refinement over curl fibres.

There is nevertheless NO relabelling-equivariant bijection to connected graphs for n≥3: the complete graph is fixed by every permutation; no nonzero alternating tournament curl tensor is. Equal cardinality does not mean equal symmetry.

[Unlabelled curl classes](root/UNLABELLED_CURL_NECKLACES.md) instead obey the necklace transform of A000568. Their first values are1,1,2,3,9,45,399,6416. Primitive necklaces recover the older October11 Lyndon/connected-automorphism-even-graph count correspondence. That correspondence and its rank codec are inherited, not claimed anew. The new interpretation is that these necklaces are actual metric-curl classes.

Hamiltonian-weighted unlabelled counts require an Adams operation: a d-fold repeated word raises its H-weight to the d-th power. Omitting this correction first fails at order6, on two copies of C3. These are exact transforms and enumeration checks, not yet a demonstrated acceleration of A000568.

All generating functions here are formal unless explicitly stated otherwise. The labelled tournament series has radius zero. The locally named series are not silently substituted for older unrelated D(z) definitions.

## 5. The Lie/topology connection has a precise obstruction

[The cocenter interpretation](root/CURL_COCENTER_AND_LIE_BOUNDARY.md) identifies the vector space of unlabelled curl classes with positive-degree HH_0 of the free associative ordered-SCC algebra. This is an exact basis map: the commutator SUBSPACE identifies cyclic words.

It is not an algebra quotient. The equal curl words ab and ba, with a a singleton and b a directed triangle, become distinct after appending the same ab: abab and baab are not rotations. Eight vertices is the first possible order for this join-congruence failure. Recovering an initial gap repairs it.

The free Lie algebra on strong tournament types has the inherited weighted-Lyndon dimensions, hence the connected-even counts. But its natural projection to the cocenter kills every bracket of letter length≥2. A matching basis count does not turn this into a Lie embedding. H and the path-colour polynomial are commutative characters and likewise cannot retain formal commutator information.

This is the session's most concrete answer to 'what do these structures represent?': ordered composition, cyclic trace and commutative response are genuinely different mathematical objects, linked by explicit forgetful maps.

## 6. Finite tensor depth has a real horizon

The metric-curl packet proves that joint Gram/alternating coefficients at root-lattice genus n−1 recover the metric pair. This is a decorated theta observable, not an ordinary scalar modular form.

At fixed order, the entire two-root response contains exactly the information in (c3,c4). In particular genus2 root-shell coefficients are universal throughout every skew-conference tournament family. At order8 a genus3 coefficient separates the two H-values189 and333 among all128 tested conference switches. This is a precise finite-depth separation, not an assertion about all genus2 full-lattice theta coefficients.

[The fixed-genus obstruction](lattice/FINITE_GENUS_OBSTRUCTION.md) goes further: for every fixed root-shell genus cutoff g, at every sufficiently large order there are strong tournaments with identical labelled score vectors and identical coefficients through that cutoff but different H. The mechanism is polynomially many bounded-depth signatures versus arbitrarily many multiplicatively independent path-count choices in Moon interval chains. The proof supplies an explicit, very large threshold, not a small effective census.

The concrete order9 pair IC(2,2,4), IC(2,3,3) has H81 and75 with the same scores and root-pair distributions. No global minimality is asserted. Full lattice theta data and growing genus remain outside this no-go.

## 7. Metagraph arithmetic goes beyond inventories

[Dyadic electrical structure](metagraph/DYADIC-ELECTRICAL-METAGRAPH.md) starts from the actual complement-line metagraph, not a different arc-flip graph. Its reflected black lines occur in twins, so its Laplacian is L_blue+2L_paired-black.

The exact critical-group2-rank is

c(blue)−c(full)+dimension(binary blue cut∩cycle space).

The next layer is a half-Laplacian pairing on blue harmonic potentials; its radical counts factors divisible by4. Here the black wiring becomes visible. Through order7 the exact2-primary groups are computed and independently checked. Degree/color-preserving synthetic rewiring can change these groups, so they provide constraints not reducible to node inventories.

A two-colour spanning-tree polynomial has a complete interval of supported exponents and exact endpoint product formulas. Its black-component contraction links the older bridge criterion to the leading small-blue-conductance term.

Boundaries matter: the order8 extension is BLUE-ONLY, not a full order8 metagraph census. It already refutes a tempting alternating-half-form pattern. Synthetic rewired graphs are not claimed tournament-realizable. No all-order connectivity, black-core connectivity, universal torsion bound or complete wiring classification follows.

## 8. Runner geometry: exact observations and invisible boundaries

[The runner packet](runner/README.md) proves that complete phase-overlap responses at speeds1,b,...,b^(J−1) lose exactly mean-zero1/b^J-periodic differences. A regular-closed interval carrier with fewer than b^J components is recovered almost everywhere against arbitrary [0,1]-valued competitors, and exactly among regular-closed carriers, once the required moments are known. The component threshold is sharp; explicit Hölder noise bounds are proved.

A finite list of scalar phase samples still cannot determine even one unknown interval. Fourier-magnitude spectra, or responses with the phase origin forgotten, lose translation. Haar measurements lose isolated points: the tight13-speed AP has six isolated weak-safe times, all invisible to every overlap integral.

Thus geometric reconstruction is not an LRC proof. It does not supply a uniform component bound, lawful translated tests, owners, common ancestry or arrival clocks. LRC remains open in this repository's audited work.

## 9. Recovery versus new work

The key recovered threads were the old path-cover substitution kernel and endpoint jet (THM3121/3134), ordered-join path-colour law (THM3166), linear forests (THM4757), the October11 SCC/Lyndon codec, automorphism-even graphs (THM4721 and the published Royle et al theorem), Hodge currents (THM4739), bicycle obstructions (THM2467), Moon interval products (THM4744), and the runner clock/endpoint guardrails. Each packet records exact source paths and hostile controls.

The advancement is in the new proved connections and sharp boundaries, not in relabelling these older theorems as discoveries. Literature priority beyond the explicitly cited sources has not been established.

## 10. Where this points next

- For complete wiring, seek a certificate which includes actual path/port compatibility, then test whether the new electrical layers reject all remaining synthetic alternatives. Counts and torsion alone are not sufficient.
- For enumeration, exploit the exact absolute-curl response quotient and modular decomposition, but benchmark against existing canonical enumeration before claiming speedup.
- For higher invariants, grow the necessary incidence/tensor depth rather than extrapolating a fixed low-order pattern.
- For runner and Collatz-style coverage, retain exact boundary and owner/clock data. This session advances the observer theorem, not the outstanding uniform153N coverage argument.

The novel direction is not that every problem is secretly a tournament. It is that one can often locate the exact missing coordinate by asking which composition or observation first distinguishes two apparently identical states.
