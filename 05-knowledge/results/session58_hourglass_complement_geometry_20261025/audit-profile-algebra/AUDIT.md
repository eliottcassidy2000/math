# Independent audit: order-averaged complement-pair profiles

Verdict: PASS for the algebraic claims and exhaustive n≤6 claims described below. Producer THEOREM.md has now been read in full and agrees with the independently established statements below. No remote actions were performed.

## Definitions and normalization

Use the uniform inner product on the Boolean cube of pair orientations. For graph X let χ_X be its Walsh character, ε(p,X)=(-1) to the number of edges reversed by p, and O_X its ordinary graph orbit, of size m_X. A signed orbit sum S_X with coefficients ±1 is well-defined exactly when every automorphism of X has sign +1 (Royle-even, not merely even edge count). Its coefficient at χ_X is fixed to +1. Write U_X=Σ_{Y∈O_X}χ_Y and α_X=(1/n!)Σ_p ε(p,X).

The invariant subspace for the affine tournament action has one line spanned by S_X for each Royle-even graph orbit, and zero lines for the other orbits. This follows immediately by propagating a coefficient through the signed monomial action; an odd stabilizer forces that coefficient to be zero. The ordinary graph invariant line is spanned by U_X. Both vectors have squared norm m_X, and their inner product is m_X α_X. Consequently, for ordinary and tournament averaging projections P_G and P_T,

P_G S_X=α_X U_X; P_T U_X=α_X S_X; P_T P_G S_X=α_X² S_X.

The first two maps have different source/target invariant spaces. Only their composition is an endomorphism with eigenvalues α_X². It would be imprecise to call the profile operator itself an endomorphism with eigenvalues α_X. The profile linear map and its adjoint have rank equal to the number of surviving nonzero α_X.

Global converse/complement multiplies χ_X by (-1)^|E(X)|. Thus, for the profile indexed by unordered graph-complement classes and tournament-converse classes, retain only even-edge X. Graph orbits failing Royle-evenness still contribute zero under P_T. Ordinary graph-complement indicators span the entire graph/complement invariant space, so no additional rank can escape this diagonal calculation.

## Independent cancellation proof and a parity family

Let G_r have vertices u=0,v=1, common neighbors 2,…,r+1, and pendant vertex r+2 attached only to u; its edges are uv, both spokes to each common neighbor, and the pendant edge. It has 2r+2 edges. Assign independent uniform [0,1] variables to vertices. Their relative order is uniform, so α is the integral of the product of edge comparison signs.

Fix x_u=a,x_v=b. Each common neighbor integrates to 1−2|a−b|; the pendant integrates to 1−2a, and the uv edge contributes sign(b−a). Hence

α(G_r)=∫∫ sign(b−a)(1−2a)(1−2|a−b|)^r da db.

Average this expression with its a↔b substitution. The first factors become |b−a|. Integrating over d=|a−b| gives

α(G_r)=2∫₀¹ d(1−d)(1−2d)^r dd
       =0 for odd r,
       =1/((r+1)(r+3)) for even r.

The odd-r cancellation is the involution d↦1−d. For even r substitute t=1−2d to obtain (1/2)∫₀¹(1−t²)t^r dt. This also verifies the normalization and sign, including α(P3)=1/3 at r=0.

For r≥2 degrees distinguish u, v, and the pendant, and all remaining automorphisms permute common neighbors. Thus Aut(G_r)=S_r. In the chosen labeling those permutations reverse no edges, so G_r is Royle-even. For odd r≥3 it therefore gives a genuinely nonzero S_X killed by P_G. At r=3 it is the six-vertex graph of mask255, eight edges and Aut=S3. At r=1 an automorphism interchanges the two degree-two triangle vertices and reverses their edge, so the graph is Royle-odd; its zero integral does not provide a nonzero signed invariant. This exception is essential.

Adding k isolates to G_3 preserves α: relative order on its six nonisolated vertices remains uniform. Its automorphisms are S3×S_k, and the isolate permutations reverse no edges. Therefore every n≥6 has a nonzero even-edge signed invariant in the kernel. The extension establishes rank deficiency for every such n, not failure to distinguish individual converse classes.

The root's separate eight-vertex example, two disjoint K1,3, also works. Each star is Royle-even and has three edges. Reversing the relative order within one component negates the product, so α=0. Component swaps are signed-even (or follow the inherited component-even lemma). This is an independent disconnected check, not the minimal obstruction.

## Independent finite census

`audit.py` imports none of the producer's code. It constructs permutation actions explicitly, exhausts all ordinary graph orbits and affine tournament/converse orbits, builds the actual integer order histograms, and computes exact rational ranks by elimination. It separately enumerates signed graph orbit coefficients and compares their nonzero count against those ranks.

n: 1,2,3,4,5,6
ordinary graph classes: 1,2,4,11,34,156
tournament/converse classes: 1,1,2,3,10,34
actual profile matrix rank: 1,1,2,3,10,32
distinct individual profiles: 1,1,2,3,10,34

Exactly two Royle-even even-edge modes vanish at n=6: masks255 and2015 in lexicographic pair order. Their automorphism group orders are6 and2. Both edge lists are included in `independent_results.json`. The actual six-vertex matrix has34 rows and78 columns. This proves vertex-minimality of a linear-kernel obstruction at n=6 within this exact profile setup. It does not prove minimal edge count among all obstructions unless that separate question is explicitly investigated.

Crucially all34 rows remain distinct despite rank32. A linear relation among distributions over classes need not be a difference of two point masses. Thus no collision between two tournament converse classes is obtained here, and no all-order completeness or incompleteness claim is justified. A nonzero zero-mean kernel function can be split into positive and negative parts to yield distinct mixtures with equal profiles, using uniform-labelled orbit weights, but that is a mixture statement.

## Checks and provenance

`python audit.py` and `python -O audit.py` produced identical outputs. All checks use explicit RuntimeError exceptions, not assertions. `family_check.py` independently evaluates the finite binomial integral expansion for r=0,…,7 and full signed permutation sums for r=0,…,5 (through eight vertices); normal and optimized outputs agree. Two-star signed order sum is independently enumerated over all8! orders.

The signed Walsh basis is inherited, not a new general basis theorem: repository THM-4721, Section5/TheoremW, already proves it, and Section4 proves the component-even lemma. The present specialization calculates the ordinary-averaging coefficient, its explicit zero family, and the finite profile rank. The earlier Session58 complementary-graph note provides the actual order-profile definition. `REVIEWED_SOURCE_HASHES.json` pins all reviewed inputs at audit time. Nothing here claims a literature-wide novelty search.

Final producer review: `check_producer.py` verifies every entry of the producer’s34×78 literal matrix against the independent census, and both displayed mixture laws against all78 independent columns. Positive and negative probabilities each sum to1 and have disjoint support. Normal and optimized checks agree. `negative_controls.py` confirms explicit rejection of false rank34, individual-collision, and nonzero-integral claims under normal and optimized execution. These negative controls test the check mechanism, not additional mathematical cases.

## Quadratic completion supplement

The producer added `quadratic_certificate.json` after the initial result. `quadratic_audit.py` reconstructs its entire34×34 selected minor from the independently generated literal profile matrix, using the32 selected original profile columns and squares of columns whose graph-pair representatives are masks3 and15. Exact rational elimination, independent of the producer's Bareiss determinant implementation, yields

211991592542375584355844096000000 ≠ 0.

The determinant is967372386 modulo1,000,000,007. Normal and optimized runs agree. Hence degree≤2 functions of the profile span all functions on the34 converse classes at n=6. This repairs linear expressivity on the finite class set, not single-class distinguishability, which already held. No general all-n degree-two result follows.

For mixtures this has an operational interpretation only with shared latent provenance: sample one tournament T from the mixture, then independently sample two orders conditional on that SAME T. The probability of observing graph-pair type q twice is E_T[(A_T(q)/6!)²]. This is not (E_T[A_T(q)/6!])². If the tournament is independently resampled for each observation, the resulting joint law is the product of the first-order mixture law and cannot distinguish the two displayed mixtures. Recording two independent observer orders conditional on a single retained tournament does supply the needed second moments. Mixed products similarly arise from two possibly different type outcomes.

More generally, on a finite set of N objects with distinct profile vectors, polynomials in their profile coordinates span all class functions: for each target x and each other y choose a coordinate differing between x and y, normalize that coordinate's affine linear function to equal1 at x and0 at y, and multiply these N−1 factors. This constructs its indicator. A profile collision prevents separation at every polynomial degree. The n=6 computation improves the general degree≤33 interpolation bound to2; it does not assert2 for larger n.

Final quadratic prose review: Section6 of the producer theorem was read and passes, including its finite Lagrange interpolation argument and same-latent interpretation. Independent direct-matrix evaluation confirms squared-count contrasts−9 and9/4 and normalized probability contrasts−1/57600 and1/230400 for masks3 and15, respectively.
