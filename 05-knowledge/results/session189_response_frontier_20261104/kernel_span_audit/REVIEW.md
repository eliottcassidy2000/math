# Independent review: selected odd-transitive kernel span

## Verdict and exact novelty

PASS for the mathematical statements of the producer PROOF.md, SHA-256 01c05852cde9c315c199f34f89211ede98a0ab24b8b9efa2f99f8aef53af6e51, subject to the exact dependency and scope recorded here. The all-order proof establishes a particular mixed-quartic battery, not merely finite nonvanishing tests.

Generic detection of a single nontrivial module converse in an arbitrary fixed exterior was already proved by incoming E4. The pinned E4 source was independently retrieved, and its provenance is retained in E4_PROVENANCE.json. The new contribution here is the complete span of selected odd-transitive kernels, an explicit finite-jet inverse, and a maximum-matching certificate of nonvanishing within that restricted battery when its stated physical size bound holds.

The theorem compares a source with reversal of the distinguished transitive module alone. It does not prove a prime-root XOR-family reduction, permit quotient-only changes without restoration, ensure the exterior mark lies within a prescribed core, or replace physical exterior names by unweighted quotient vertices. These limitations are essential, not cosmetic.

## 1. Span and finite jet

The exact two-channel identity from the separately reviewed odd-transitive packet is the sole structural response input. Its variables retain physical names in a squarefree algebra, so repeated names annihilate products. There is no endpoint-marginal substitution.

The finite-difference transform sends binom(s,j) to t^j and is invertible on every bounded-degree polynomial space. Vandermonde's identity gives

    E(P(-1-s))(t)=E(P)(-1-t).

For odd L=m-2, each P_u is fixed by s↦-1-s. The selected polynomials P_1,...,P_(k-1),P_L have different lowest nonzero powers: P_L has constant 1, and P_u has leading low-degree term -s^u for 1≤u≤k-1. Thus they are independent in the k-dimensional reflection-invariant space. Invertibility and equivariance send this entire space to Q[t(1+t)] through degree k-1. This proves full span, rather than inferring it from computed ranks.

The raw moment cutoff also checks exactly. Every tau monomial uses at least two names, the singleton terms of alpha S_i-beta E_i cancel, and alpha beta(ell_i-rho_i) has support degree at least one. Consequently the jth raw moment requires at least 2j-1 names, and vanishes for j>floor((n+1)/2).

Writing w=t+t², the coefficients of w^h start with t^h with coefficient one. Truncation at degree d is unitriangular. The explicit Lagrange-inversion coefficient

    (-1)^(h-j) (j/h) binom(2h-j-1,h-j)

does give t^j modulo t^(d+1). The argument concerns a finite nilpotent jet, so it invokes neither analytic convergence nor a limit of finite experiments. Rational combinations can be multiplied by a common denominator to yield integer certificates.

The independent polynomial controls retain every kernel, its exact w coordinates, and inverse rational certificates for odd m=3,...,51. Explicit inverse-jet formulas are checked through d=30. These are implementation checks of the all-order algebra above.

## 2. Maximum matching and constructive strictness

The directed O-to-I cross-edge graph is formed on physical exterior names. Each O-to-I path contains an O-to-I edge. Disjoint paths yield disjoint selected edges, so tau^(q+1)=0 for maximum matching size q. Grouped quotient counts would not suffice for this step.

For a nonperfect maximum matching, the (q+1)st moment loses its A_i term. If i lies in O, the ell_i channel also vanishes by the same matching obstruction; if i lies in I, the rho_i channel vanishes. The remaining channel has a fixed sign on each side. A maximum matching leaving i unmatched supplies q disjoint edge paths, the singleton i, and Hamilton paths on the leftover I and O vertices. Their supports are disjoint and cover the exterior. Ordered tau factors give the stated lower bound (q+1)! with no extra factorial. This works also for q=0, including an empty incidence side.

In the perfect case n=2q, the A_i term cannot fit. For i in O, a nonsingleton O-to-O rho_i path together with q-1 cross paths would require q+1 distinct O names and is impossible. The two remaining support patterns are exactly: q-1 cross edges plus an extra singleton I name from alpha, or q-1 cross paths covering all names except i. Subtracting the positive cross-edge ell_i channel produces the stated negative formula. The I case is its symmetric counterpart.

For q≥2, a directed edge between two matched I endpoints extends one matching edge into a three-name path; omitting the partner O name yields a positive T_i term. The O-edge construction supplies an I mark. Applied to any two specified marks on one side, one is strict. The bound q! follows from ordering q-1 disjoint tau paths and the prefactor q. For q=1, direct substitution gives constant moments +2 and -2; this exceptional index is handled separately and correctly.

This proves a nonzero moment in every nonempty exterior. Its index lies within the stated jet dimension. All signs needed for grouping are on physical incidence sides; groups must not mix O and I.

## 3. Actual named and grouped coefficients

The selected physical coefficient is [a_0 a_1 c_r c_i], with four distinct names and fixed port edge x→y. For a uniform-incidence exterior group G, identifying its c bits gives exactly the sum over i in G: the coefficient already saturates both port degrees, so no higher-support monomial can alias into it.

A nonzero rational moment combination forces at least one selected coefficient to be nonzero. A four-variable Boolean finite difference uses at most sixteen actual attachment queries. The size condition m≥2 floor((n+1)/2)+3 is sufficient because it supplies every moment through the nilpotent cutoff. The graph-dependent bounds follow by the stronger matching nilpotence. No necessity or sharpness of the general size bound is claimed.

The independent implementation uses a dense squarefree convolution and an endpoint subset recurrence, rather than the producer's sparse multiplication and path-permutation enumeration. It exhausts all 1,098 exterior/root-incidence cases through exterior order four and adds deterministic controls through order eight, for 1,375 cases total. It enumerates all maximum matchings and checks the exact signed formulas, factorial lower bounds, finite support cutoff, and all selected saturated quartics.

For the producer's seven actual examples, a separate literal old-permutation/port-placement engine retains complete source and local-converse named polynomials, including all zeros: fourteen graphs and 26,822 coefficients. All 42 selected physical quartic contrasts agree. This engine expands the actual port-edge indicators directly and does not use the two-channel identity or three-cover formula.

A third calculation builds the actual augmented tournament and counts Hamilton paths by a predecessor subset recurrence. All 960 source/target Boolean query counts agree with the complete literal polynomials and reproduce all 30 grouped finite differences. Thus the grouping statement is checked at actual physical masks, including groups with multiple names, not merely at formal sums.

## 4. Boundary controls and source-relative scope

Six additional complete named banks (996 coefficients) retain two useful boundaries:

- A transitive triple with one exterior name has its selected mixed quartic zero, even though its full response differs under local converse. The restricted battery's size condition is a real issue and cannot be dropped merely by citing generic E4 separation.
- For a transitive triple with the two exterior names forming the directed contracted triangle, the two selected quartics are +4 and -4. Grouping both incidence sides together cancels them. The raw positive moments also vanish, exactly as the producer states.

Same-side grouping preserves a signed witness; arbitrary grouping does not. Further, the surviving exterior mark can lie outside a designated core. Therefore this proof alone does not justify a root-family conclusion where ports must imitate that core throughout an outside ambient. No such conclusion is approved here.

## 5. Reproducibility and corrected receipt issue

The first clean producer replay found one reproducibility defect: a prose novelty field had been added to the frozen receipts but was not emitted by verify.py. It did not affect any mathematical value. The author corrected the verifier to emit that field and regenerated both modes. The corrected verifier has SHA-256 3bad86d48997389c7aade3c6a3939ee3df5bbab24dc18c7e40e19ed11d2606b2, and both corrected producer receipts have SHA-256 9fd2170aac6e2e8107031c4c16d2580071647e5e11d1d654645f653b012213b0.

The final release pins all producer files and the inherited odd-transitive proof and verifier. Independent scripts, retained exact banks, actual queries, source receipts, and normal/optimized outputs are covered by PUBLICATION_FILES.json and SHA256SUMS. Executables and temporary replay directories are excluded. Clean producer replays in both Python modes regenerate the frozen receipt byte-for-byte, and both audit modes agree. Every listed hash is checked before release.
