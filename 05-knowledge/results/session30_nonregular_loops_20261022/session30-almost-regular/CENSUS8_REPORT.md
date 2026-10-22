# Complete almost-regular order-eight audit

## Scope and result

FINITE-EXACT universe: all85 isomorphism classes of almost-regular tournaments on8vertices, without identifying converses unless they are genuinely isomorphic. Their labelled mass is3,230,080. Every Hamiltonian count was computed by integer endpointDP and independently checked against a signed-path/linear-forest construction. Predictor arithmetic is double precision; it is not an interval or uniform asymptotic certificate.

Second order improves absolute y error on73 of85 classes and worsens it on12. Across unweighted isomorphism classes, first/second mean absolute errors are0.0830701/0.0383599 and RMSEs0.0964332/0.0478943. Labelled-weighted RMSEs are0.0947256/0.0448158. Thus the whole-universe improvement is real, but neither universal pointwise improvement nor universal ranking is true.

## Completeness certificate, independent of a source bank

For an ordered score sequence d, let F(d) be the number of labelled tournaments with that sequence. Removing a distinguished last vertex of outdegree d_last gives the exact recurrence: choose its outneighbour set O of size d_last, decrement every remaining vertex's required outdegree by1 iff it is not in O, and sum F on the resulting sequence. Sorting each remaining sequence is safe because relabelling gives a bijection. Reject negative scores, scores at least the remaining order, and incorrect degree sums; F(empty)=1.

This tiny memoized recurrence computes F(4,4,4,4,4,4,4,4,4)=3,230,080. Candidate regular9 representatives were discovered by deterministic circulants, C3[C3], and49 fixed-seed chain samples. Discovery randomness proves nothing, but the exact mass does: the15 pairwise nonisomorphic representatives have automorphism mass Σ9!/|Aut|=3,230,080. Therefore no missing isomorphism class exists. Isomorphism and automorphisms are computed by full edge-consistent backtracking with invariant pruning, not a spectral equality shortcut.

For each regular9 class, delete every vertex. Deduplication gives85 pairwise nonisomorphic almost-regular8 classes, whose independent automorphism mass Σ8!/|Aut| again equals3,230,080.

Why the masses are equal: every almost-regular8 tournament uniquely extends to a regular9 tournament by adding a distinguished vertex v which beats each score4 vertex and loses to each score3 vertex. There are four of each score. Conversely, deleting v from a regular9 tournament produces such an almost-regular8 tournament. This is a labelled bijection, so the second exact mass also certifies completeness directly.

All matrices, automorphism orders and mass totals are in census8_bank.json. No external graph generator, supplied class bank, or claimed universal mixing of the chain is assumed.

## Hamiltonian spectrum, by isomorphism class

621:2,627:2,629:1,631:4,633:2,637:1,639:4,641:8,643:6,645:4,649:4,651:10,653:12,657:15,659:4,661:6.

In particular, the interval between the minimum and maximum is not filled by all odd values. This is a spectrum for a fixed almost-regular order, not a claim about the unrestricted tournament Hamiltonian-count spectrum.

## Ranking audit

There are3263 unordered pairs of classes with distinct exact H. Ranking byρ produces1998 reversed pairs,181 exact ties and1084 correct pairs. Ranking by full first-order predictedκ produces636 reversed pairs; full second order produces465. No first/second numerical tie is within tolerance1e−10 for these distinct-H pairs.

Forρ, ties are decided exactly using the integer characteristic polynomial: ρ=P_S(8)/|P_S(8i)|, a rational number. The raw floatingρ ranking in the initial census8_summary.json is superseded by census8_ranking.json, which is the authoritative ranking audit. Floating rounding must not turn exact ties into wins or losses.

## Worst cases and direct formula control

Worst first-order absolute y error: class5,H659, exacty1.1830213647, first1.3756999272, second1.2818204610. Worst second-order absolute error: class60,H621, exacty1.5178474184, first1.3855395567, second1.3866388351. Its signed forest coefficients are(1,0,48,0,288,0,−1152,0). This hostile case retains error0.13120858 after second order; a finite predictor is not exact reconstruction.

The engine's complete first-order prediction was compared on all85 cases against this packet's independent closed formula y_regular+tree−n·lollipop. Maximum discrepancy4.22e−15. This checks conventions and first-order numerical implementation, not the separate second-order engine derivation.

## Reproduction

Run census8.py, predict_census8.py, census8_ranking.py. The first script is independent of producer machinery. The second uses the sibling decorated engine only for predictions; exact H and all forest coefficients come from independently implemented algorithms. The third uses exact integer characteristic polynomials forρ comparisons. Files contain every bit string, exact H, automorphism order, forest coefficient and predictor component.
