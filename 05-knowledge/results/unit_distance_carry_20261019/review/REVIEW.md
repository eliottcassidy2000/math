# Independent review: unit-distance lifting checkpoint

Reviewed 2026-10-19. Outcome: **PASS**, within the hypotheses and scope stated below. No claim about a new chromatic number of the Euclidean plane is established by this checkpoint.

## Mathematical findings

1. **Joint finite quotients.** Theorem 1 correctly proves both ingredients that are needed: every tuple of local norm-one steps lifts to an exact global norm-one element, and the additive component has full vertex image. The local construction `w=t+u c(t)`, followed by global CRT and `w/c(w)`, establishes the first ingredient without assuming that arbitrary vertex CRT solutions have norm one. The difference-span argument and `H+pA=A` establish the second. Projection to a deepest factor and the compatible diagonal have opposite graph-homomorphism directions and give the asserted equality of chromatic numbers. This concerns colorings factoring through the finite observations, not a lower bound on the global graph.

2. **Carry obstruction and loss bound.** The trace contradiction is correct for every odd residue characteristic and every prime-power residue field, with an unramified base. It also applies to the straight `p^r`-cycle between depths `r` and `r+1`. The partition into straight cycles proves the failed-edge bound at every depth. At depth one it gives 435 failures among 12,615 edges for q=29. The q=3 example correctly demonstrates that all successive sections can fail while the chromatic number remains 3 at all depths. Thus failure of a section is not a chromatic decrease or failure of an unrelated embedding.

3. **Stress and local lifting.** The carry/rigidity equation is exact modulo p². Its solvability is equivalent to annihilation by every left-kernel stress, and its cokernel class is independent of the initial representatives. Full row rank allows multivariable Hensel lifting. This yields an exact local framework; it does not by itself give coordinates in the prescribed global number field or a Euclidean chromatic consequence.

4. **Hidden poles.** The CRT construction in Theorem 3 has the claimed valuations and norm. Its congruence is semilocal: `v_Pi(u−1)≥n_i`, not necessarily `u−1∈I` in the global integer ring. The theorem permits other poles. The explicit fixed-place integrality sentence and semilocal interpretation requested during review were added to the checkpoint.

5. **Application at 29.** The cited field census is consistent with the direct factorization `(X−3)(X²+5X+11)` modulo 29; the remaining quadratic has nonsquare discriminant 10. There are eight conjugation-fixed degree-two completions and sixteen swapped places. Only the eight fixed completions carry the one-completion quadratic involution used in Theorem 1. The swapped pairs supply the pole examples for Theorem 3.

6. **CM baseline.** The Kronecker argument, valuation-signature equivalence modulo roots of unity, finite denominator-budget bound, and Gaussian sharpness example are correct with the CM assumption. The coloring classification is properly cited to [Davidson, Theorem 3.3 and Corollary 3.4](https://arxiv.org/html/2610.09474v1#S3.SS2), whose statements and elementary proof were source-checked. The non-CM counterexample explains why that assumption cannot be omitted.

## Reproduction and independent controls

- The original `reproduce.py` was copied to this separate review directory and run there. Its JSON and stdout matched the checkpoint's saved files byte-for-byte. The reviewed script SHA-256 was `35ced69c62daa41e59c79619885e8c05cbf23d39ae52bdd9f7602ad3da0800ca`.
- `independent_checks.py` was written separately and imports no checkpoint implementation. Run `python3 independent_checks.py`; it uses only the standard library and the included `section-witnesses.json` and `rigidity-witness.json`.
- It enumerates local norm-one quotients at p=2,3,5,7, including depth three at p=2; checks full additive generation at mixed depths (1,2) for p=2 and p=3; and checks a genuine q=9 residue-field carry example.
- It verifies all candidate edges using exact integer norms, independently certifying minimum section defects 6 at q=3 and 15 at q=5, with one failed edge in every straight cycle. It also checks the higher-depth obstruction for (p,r)=(3,2),(3,3),(5,2), and the q=3 three-coloring at the tested depths.
- Results are saved in `independent-results.json` and `independent.out`. These finite controls support, rather than replace, the proofs.

## External-status boundary

The checkpoint accurately distinguishes its elementary arguments from cited global results. OpenAI's [family-158 scope document](https://github.com/openai/math/blob/main/lean/docs/158.md) reports formalized bounds 6≤χ(R²)≤7 for arbitrary colorings. This independent review did not rebuild or independently audit that formalization. The [external local-coloring note, §7](https://github.com/decalion89/chromatic-number-of-the-plane/blob/main/notes/local_colourings.md) does report inconsistent section systems at q=5,11,17. No exhaustive novelty claim is supported or needed.


## Addendum: stress-free ceiling and sharpness control

**PASS** for the final Theorem 4 and revised next problem. For any subframework on at least two distinct vertices, two translations and the infinitesimal rotation are three linearly independent kernel vectors over the fixed residue field. Hence rank is at most `2v−3`. Row independence is hereditary, so every subgraph has fewer than `2v` edges and therefore a vertex of degree at most three. This proves 3-degeneracy and the four-color ceiling. For a k-critical n-vertex subgraph, minimum degree at least `k−1` gives stress dimension at least `ceil((k−5)n/2)+3`. In particular, a five- or six-chromatic framework cannot use full-row-rank Hensel lifting as its sole lifting mechanism. The revised next problem correctly requires compatible nonzero stress spaces and subsequent carries.

The writer's new sharpness verifier reproduced byte-for-byte in this separate directory. A separately implemented check was also added to `independent_checks.py`: it recomputes the eleven edges of the seven-point F_17² Moser spindle, verifies failure of all 2,187 three-color assignments and the supplied four-coloring, and checks all induced subgraphs. Its independently selected rigidity minor using columns 0 through 10 has determinant **15 modulo 17**, certifying rank eleven. Direct norm checks validate the supplied lifts modulo **289** and **4913**. Results are included in the refreshed `independent-results.json`; the exact witness data are retained in `rigidity-witness.json`.

The final coding formulation is also exact: every lift has error vector `c+Rh`, every correction vector is available, and failed equations are precisely its nonzero coordinates. Minimum section defect is therefore minimum Hamming weight in the affine code coset `c+im(R)`, with stress pairings giving its syndrome. The `|E|/p` cycle bound applies to the full residue graph, not an arbitrary selected framework; this qualification was requested explicitly. Sharpness at q=29 remains an open finite problem.
