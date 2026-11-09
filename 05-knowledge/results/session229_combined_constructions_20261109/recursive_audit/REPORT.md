# Independent recursive/template audit

**PASS: mathematical proofs and stated finite claims, with a small total-algorithm boundary convention.** Chan–Pak completeness remains OPEN. No Git operations were performed.

Audited source: sibling `recursive/PROOF.md` and `recursive/verify.py` in the session229 deliverable (original scratch folder `session229-recursive`). Dependencies: THM-3121-path-cover-walk-content-substitution-kernel; THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform; THM-3166-falling-factorial-order-join-path-colour-transform; session228 recursive compiler.

## Proof review

1. **Sharp arbitrary-C3 order bound: PROVED.** The eight positive run-count cells with each coordinate in {1,2} give exactly 3+Σa_i+Σa_i a_j+3a_1a_2a_3, for a_i=2^{n_i}−2. Subtracting 2^{n−2}+1 gives (Σa_i a_j)/2+11a_1a_2a_3/4. Transitive blocks of sizes (n−2,1,1) realize equality for every n≥3. Thus the prior cap remains valid without a singleton requirement.
2. **General cycle-certified finite-template bound: PROVED.** The q one-run-per-block words and q choices of the one repeated block have distinct maximal-block-run signatures. Decorating repeated block i contributes F_i(2)∏_{j≠i}h_j. The resulting sum is Σ_i(h_i+F_i(2))∏_{j≠i}h_j≥Σ_i(2^{n_i}−1). Pairwise balancing proves the exact integer lower envelope L_q(n). No connectivity-to-Hamilton-cycle theorem is needed because the cycle is an explicit certificate.
3. **Target-bounded membership/minimum-depth compiler: PROVED as an algorithmic theorem, not implemented as a full expanded-language census here.** Every template output is bounded in order, descendants cannot have larger H, and full path-cover profiles determine future substitutions. Order induction captures all bounded profiles. Flattening only exterior joins retains arbitrarily large padded outputs; finite seeds above the cap are retained as exterior factors. The min/max depth recurrence is exact after deduplication by full profile. General-cap boundary: if X<q, that template has no admissible node and must be discarded; if no template remains, use cap zero and only multiplicative closure of eligible seed scalars. This convention is included in the frozen producer checkpoint1.
4. **Three-moment family: PROVED.** Substitution with two TT2 blocks gives 11H+17F(2)+4F(3), equivalently 4Q(3)+5Q(2)−11H. Hence TT_m gives 4·3^m+5·2^m−11, including 137 for m=3. The negative coefficient is merely a basis change, not subtraction of physical objects.
5. **19 remains excluded: PROVED relative to the specified seed grammar.** A C3 node with all three blocks of order at least two has H≥45. Therefore any expression of scalar19 in the enlarged language would use only singleton-C3 nodes, joins, and the original seeds. The complete old compiler excludes19. This is not tournament nonexistence.

## Independent exact computation

`audit.py` imports no producer or repository module and reads no input files. Its seed engine differs from the producer: it counts induced Hamilton paths by endpoint subset DP, then forms unordered path covers by selecting the component containing the least-labelled remaining vertex. Profiles are retained in the falling-factorial basis, not the producer's ordinary-power basis. Newton differences reconstruct exact joined/cyclic profiles.

- All 1,099 labelled seed tournaments through order5 give15 distinct profiles.
- Independently reproduces old complete membership through4095:577 scalars; depth counts79,390,107,1 at depths0,1,2,3. Scalar/depth SHA256 is `7c895b5f3a40b78ab702a2d63b0cd43daccc6fe39f4d353006f98c331b44fa50`, matching the session228 receipt.
- Every one of the1,331 ordered triples of labelled tournaments of orders1..3 passes literal Hamilton DP against the inherited C3 kernel.
- The family identity passes all15 seed-profile representatives with literal graphs.
- Exactly17 scalar escapes occur among single C3 nodes with all children old profiles of orders≤5:137,159,483,509,1121,1669,1685,1779,1959,2097,2233,2239,2283,2415,2721,2823,2961. Each is independently reconstructed and counted on a literal graph; outgoing rows are exported in the receipt.
- All26 cycle-certified labelled templates of orders3,4, with each block TT1 or TT2, give400 literal checks of the general bound and balancing envelope.
- All512 transitive C3 size triples in1..8 satisfy the sharp lower bound; all applicable singleton-pair cases attain it.

Normal and optimized executions use explicit exception checks, not assertions. Their receipts are byte-identical. These are finite controls alongside the proofs, not extrapolations to all templates or all odd targets. The17 additions are not a complete census of the enlarged recursive language, nor claims of new tournament attainability in the literature.

## Reproduction

From this directory:

    python3 -I -B audit.py > replay.json
    cmp replay.json CHECK.json
    python3 -I -B -O audit.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json
    sha256sum -c SHA256SUMS
