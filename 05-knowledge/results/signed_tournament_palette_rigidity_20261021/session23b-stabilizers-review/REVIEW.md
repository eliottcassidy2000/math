# Independent review of Session 23b stabilizer/palette packet

Verdict: **PASS for the exact REPORT.md SHA-256 below.** No blocking mathematical or computational defect found. This is an independent-agent review, not the producer's second implementation. All review outputs were written outside the repository; no Git operations or repository modifications were performed.

Audited report:

`d197f61808e4cd65d8f26d27ad5663d793c4caf1f8b6133d948914b2dc5405d7`

Uploaded primary PDF:

`b9b484ab2e28417b3f6f2c696862e10cd73ebe36569c2037eda58a4320d3f01c`

## What passed

1. **Imported theorem scope.** The actual uploaded Gluck paper's Theorem 1 is about primitive solvable permutation groups, not arbitrary primitive or arbitrary solvable groups. The complete exceptional list in the report agrees with printed page 65, including the degree-nine order-eight dihedral case. Printed pages 65–66 were checked visually. Theorem 2's maximal chain, permutation images, and exclusion of index two are retained correctly. Corollary 1 does apply to all faithful odd-order permutation groups, including intransitive groups by taking their orbit images and intersecting their kernels. The paper explicitly invokes solvability of odd-order groups; the report does not mislabel this dependency as elementary.

2. **Paper details.** The regular subset is setwise, not pointwise; complements preserve its stabilizer. Proposition 1 has m=2^n and group order 2^(m−1), not a linear exponent typo. Corollary 2's monomial-module/odd-quotient/faithfulness hypotheses and merely abelian, not necessarily trivial, vector stabilizer are stated correctly. The publisher confirms original publication in 1983, with 2018 only the online date. [Publisher](https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/trivial-setstabilizers-in-finite-permutation-groups/9667B85D9DDD0FE911C8382287B86167).

3. **Existing literature.** The ordinary tournament two-color statement and general minority bound are not new. The cited Meslem–Sopena abstract explicitly derives its cyclic-tournament result from Gluck and states the canonical Paley result; the cited Lozano abstract explicitly supplies the bound and extremizers. Lozano's journal/year/pages also agree with the journal issue. [Meslem–Sopena](https://arxiv.org/abs/1608.04866), [Lozano](https://arxiv.org/abs/1707.05549), [journal issue](https://www.iam.fmph.uniba.sk/amuc/ojs/index.php/amuc/issue/view/34). This audit does not certify absence from all literature or any priority claim.

4. **Residual path transfer.** The ordinary path stabilizer fixes every marked vertex, so restriction to the omitted carrier is faithful. Its odd order follows from the absence of tournament involutions. Gluck gives a regular subset, and complementation gives the minority bound. Empty and singleton omitted carriers cause no ordinary-action defect. The regular-color fiber over one uncolored isomorphism class is R/|H|, not 2^m/|H|. The same proof applies to individually marked scaffolds but does not erase an unfaithful kernel.

5. **Signed three-color theorem.** A regular binary coloring for the odd ordinary subgroup leaves a signed stabilizer injecting into C2. Any remaining negative element is an involution. Recoloring one of its moved vertices with a new named color refines the original partition and kills it. For n≥2 the vertex representation determines the sign faithfully. For n=0 or 1 a formal negative identity is invisible and cannot be killed by any vertex colors. For reversed-converse marked paths the omitted action is faithful when m≥2 because an anti-map cannot fix two vertices. The m=0,1 warnings and the odd-k/m=1 exclusion are correct.

6. **Disjoint-product transfer.** The inherited Corollary 5C really is a disjoint-support permutation-action factorization, not merely an abstract direct product. A C3 triangle factor has exactly three regular singletons and three regular doubletons. Thus r_T(z)=r_H(z)(3z+3z²)^t(1+z)^u, R_T=R_H6^t2^u, rho(T)=rho(H)+t, and the stated free-orbit quotient follow. The size bound uses |B|≤2delta only after applying the minority bound on B. Ordered joins of arbitrarily many directed triangles correctly show unbounded minority size at delta=0.

7. **Fano/Paley actions.** The actual permutation groups have orders 21, 42, 168 and distinguishing numbers 2, 3, 4. The two cyclic Fano systems are different and disjoint; the anti-coset swaps them, so the displayed groups are not an inclusion chain 21<42<168 in those actions. The vector identification (1,2,4,3,6,7,5) is valid. The elementary Fano three-color obstruction covers all four nonempty class-size partitions correctly. The degree-eight solvable affine-semilinear group of order 168 is a different group/action and does not import Fano conclusions.

8. **Equivariance and color semantics.** An isomorphism-equivariant selection of a single subset is fixed by the entire original automorphism group; it cannot be regular when that group is nontrivial. Keeping the regular-color family or its orbit is a different lawful output. The report correctly separates distinguishing colors from acyclic/dichromatic colors. No properness condition is imposed here: a proper coloring of the underlying complete graph would require n colors. These predicates must not be interchanged just because some numerical values coincide.

## Reproduction of the packet

Both supplied scripts were run under ordinary Python and `python3 -O`, with all outputs redirected to this review directory. The four generated JSON receipts match their four supplied counterparts byte-for-byte. All seven hashes in the supplied manifest match their files. The scripts use explicit exceptions, not removable assertions. The original 1,283 partial-path-object count was reproduced.

The producer's independent-representation script necessarily reads its producer receipt as comparison data, but does not import the producer's algorithm. Its independence claim is appropriately limited in the report. The reviewer code below adds an independent agent and an algorithm not used in either supplied script.

## Additional reviewer computation

`reviewer_check.py` imports no producer module or receipt. It generates the two affine groups from generators and GL(3,2) from coordinate swaps, cycling and a shear. Rather than enumerating named colorings or set partitions, it uses inclusion–exclusion on the fixed-color sets of prime-order elements. A nontrivial stabilizer contains a prime-order element. Each fixed-color condition identifies coordinates along cycles; intersections join those equivalence relations. This yields the full polynomial counting regular named q-colorings:

- Ordinary Paley: q^7−7q^3+6q.
- Signed Paley: q^7−7q^4−7q^3+7q^2+6q.
- Full Fano: q^7−21q^5+28q^4+56q^3−112q^2+48q.

At q=1,2,3,4 the three sequences are respectively (0,84,2016,15960), (0,0,1512,14280), and (0,0,0,4032). Thus the report's sharp 84/1512/4032 counts and corresponding 4/36/24 regular-orbit counts are confirmed by a different counting method.

The same reviewer script checks every ordered nonempty directed path in every labeled tournament of orders 1 through 5: 1+6+60+1216+53760 = **55,043 objects**. It verifies ordinary restriction faithfulness, odd residual order, binary regular-subset existence and minority bound; signed faithfulness for at least two omitted vertices; at most one surviving anti-map for every ordinary-regular binary coloring; its involution property; and successful third-color refinement. Signed kernels appear only at zero or one omitted vertex. The odd marked-length/singleton-omission exclusion is also checked. This is finite corroboration, not a replacement for the general proofs.

Ordinary and optimized reviewer runs produce byte-identical receipts. Reproduction:

    python3 reviewer_check.py --output reviewer-receipt.json
    python3 -O reviewer_check.py --output reviewer-receipt-O.json
    cmp reviewer-receipt.json reviewer-receipt-O.json

## Scope boundaries

The repository's full fixed-slack skeleton theorem was not reproved from scratch; its stated action factorization and proof were inspected, and the report explicitly makes its transfer conditional on that inherited theorem. The current research does not settle a flagship conjecture, establish global literature priority, produce a canonical equivariant symmetry-breaking representative, or alter uncolored Burnside counts without retaining the marked fibers. No changes to the audited report are required for PASS.
