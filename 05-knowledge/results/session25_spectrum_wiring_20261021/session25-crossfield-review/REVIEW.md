# Independent review: exact clocks, signed actions, and palette interfaces

## Verdict

ACCEPT the crossfield certificates, including the exact hard-inequality-skeleton reduction, as elementary theorems under their stated interface hypotheses. ACCEPT the paper bridge only in the corrected copy `PAPER_CONNECTIONS.md` in this review directory. The original has a misleading arbitrary-orbit/fixed-reflection sentence and an incorrect A.9 visual-source page reference. No repository files or producer originals were edited. No conclusion about BB halting, ordinary Collatz convergence, the native 153N ledger, or the plane's chromatic number is certified here.

The source PDF was independently extracted, and printed/PDF page18 was rendered and visually inspected. All ten old-family rules, seven merged guards, and the initialization agree with the source. Raw-machine correctness and Lean build success remain reported source claims, not independently verified results.

## Paper interface proof audit

1. Quotient: set Q=q(M). Necessity of qG being constant on each q-fibre follows from commutation; conversely define the induced map by that common value. Its definition is unique on the image. Terminal status must be constant on fibres. A quotient preserving the assigned duration labels additionally needs tau constant on fibres; endpoint equality alone does not establish a unique time.
2. Stuttering: adjoin absorbing halt symbols, use total domain-closed G, require halt reflection and strictly positive integer durations, and verify the initial raw prefix. Inducting the local simulation equality gives equality at cumulative times. If the raw trajectory halts, a cumulative time eventually exceeds its halting time, so an encoded macro state is terminal by reflection. Conversely any macro halt encodes a raw halt. Injectivity is unnecessary. Zero durations break this argument: F can reach halt while G stays at a nonterminal encoded point forever with tau=0. No raw-time formula follows from an unweighted reachability rule.
3. Binary encoding: (i,a,b,c) maps bijectively to (2a+i,b,c), inverse parity and integer division. The two c=0,b>0 branches both become A+3; the two b,c>0 branches both become A+6; the b=0 exceptional branches give exactly the printed rules. Hence the conjugacy is an all-parameter branch calculation, not a consequence of testing. Dropping i instead fails already at (a,b,c)=(0,0,1).
4. Signed group: B^3=-I, J^2=I and JBJ=B^-1 imply twelve distinct normal forms B^k J^e; determinants separate rotations/reflections. The commuting subgroup consists precisely of six rotations. The central factor <B^3> and the S3 factor <B^2,J> commute, intersect trivially, and generate. However fixed J does not stabilize arbitrary B-orbits. For v=(1,2), its image (3,-2) is absent from the six listed orbit points. The repaired text fixes v=(1,0) for the displayed action and transports reflection by P_v=[v,Bv] otherwise. Its determinant x²+xy+y² is positive for nonzero real/rational v and P_v commutes with B.
5. Stronger hostile characterization: fixed J stabilizes the B-orbit of v iff Jv=B^k v for some k. The six rank-one equations give respectively y=0, x+2y=0, x+y=0, 2x+y=0, x=0, x−y=0. Thus the fixed-reflection claim fails generically, not just in one sample. This characterization was checked symbolically and is not needed by the corrected theorem.
6. Tournament qualification: a nonidentity order-two automorphism swaps two vertices and contradicts their directed edge. Cauchy's theorem therefore excludes even-order ordinary tournament automorphism groups. Reflections are anti-automorphisms, not commuting directed symmetries.

## Palette proof audit and edge cases

The tree-decomposition proof of tw(H+uv)≤tw(H)+1 is valid: add u along a tree path from a bag containing u to one containing v. Connectivity for u is preserved, all other connected bag sets stay unchanged, and uv is covered. Disconnected graphs can use a connected decomposition tree with empty/intermediate bags. The usual bag-removal argument supplies the coloring bound. Permuting a separating leaf coloring realizes any desired ordered unequal pair; disjoint interiors and absence of extra edges make simultaneous extension valid.

The stronger reduction is exact, not merely one-way. With every leaf admitting a separating coloring, color-name symmetry leaves only two supports: unequal pairs only, or all pairs. Retain just the former as H_D. Restriction and independent extension identify the entire skeleton-vertex coloring relation with proper colorings of H_D. Equality forcing is therefore equivalent to H_D being k-colorable but H_D+uv not k-colorable. If chi(H_D)≤k−1, recolor v with an unused kth color; consequently a nonvacuous equality assembly has chi(H_D)=k and tw(H_D)≥k−1.

Scope controls:
- Distinct marked terminals and distinct leaf ports are necessary. Coincident terminals force equality trivially and lie outside the theorem.
- All skeleton vertices are retained, so isolated terminals and an edgeless skeleton cause no exception. An empty vertex set cannot have two distinct marked terminals.
- The finite simple skeleton excludes loops and parallel edges. A loop would conflict with distinct ports. Parallel leaf constraints can be treated by intersection if explicitly generalized, but no such generalization is needed here.
- Any loop inside a leaf prevents a proper coloring and violates the separating-leaf hypothesis. Duplicate edges do not change ordinary coloring support.
- Uncolorable leaves and uncolorable assemblies are excluded by the hypotheses/nonvacuity condition. The conclusion explicitly constructs a coloring, so it does not rely on a vacuous universal statement.
- Arbitrary precolors, lists, or pinned color names can defeat the permutation argument and are explicitly excluded. A chosen skeleton coloring can be extended, but the theorem does not promise compatibility with extra external precoloring constraints.
- Interior overlap, unintended vertex identification, or cross-fragment edges invalidate independent extension. Plane embeddings must certify this interface separately.
- A unit-distance graph can specify selected unit edges. If the intended graph includes every incidental unit-distance edge, those extra constraints must be included in the assembly analysis. No Euclidean realization is supplied by the abstract sharpness graph K_(k+1)−uv; for k=6 it contains K6 and is not a planar unit-distance realization with distinct vertices.
- Sharpness is algebraic for every k≥2: K_(k+1)−uv has a K_k, admits two size-k bags, and has exactly k! proper k-colorings, all identifying u,v. The k=2 instance is the length-two path and correctly sits above the width-zero threshold.

## Clock and return-map audit

The ordinary Collatz family n_j=(2^(2+12j)−1)/3 is odd and congruent to1 modulo105 because 2^12≡1 modulo315. Its next odd state is1 and its exact ordinary time is3+12j, unbounded while both modular guards and valuation modulo6 remain fixed. j=0 is a return from1; j≥1 is first arrival, as the intermediate powers of two exceed1. This proves the information-loss claim without asserting anything about the repository's different native map.

For A.9, at c≥2, write a=2k+e. The first rule reaches (2+3e,k+3,c−2). Drain the last two coordinates until one vanishes. If c−2<k+3, one final c=0 rule yields first return after c rules; otherwise return occurs after k+4 rules. These are exactly the two printed formulas, including the equality boundary c=k+5. At c=1 an odd a returns in one rule; even a returns after three, with the k=0 special endpoint. c=0 halts after one rule. Every intervening b is positive, proving first return, not just eventual reachability. Initialization (5,2,0) reaches (8,0,1) in one macro rule. L counts only A.9 rule applications. It is neither raw TM time nor classical/native Collatz time.

## Independent executable evidence

`independent_checks.py` imports no producer module and raises explicit exceptions, including under python -O. It checks exact symbolic group identities and reflection transport, the arbitrary-orbit hostile, 44,521 independently simulated A.9 first returns, 123 applicable bounded-width separation cases drawn from all labelled graphs on two to four vertices, 54 explicit expanded-fragment support comparisons, and 101 exact integer clock examples. Normal and optimized receipts agree. The original producer checkers were also replayed. Finite controls corroborate the proofs and do not replace them.

The final manifest freezes the source PDF, producer note/checker/receipts, corrected note, crossfield final note/checker/receipts, and review artifacts. Review acceptance applies only to those byte versions.
