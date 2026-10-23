# Independent audit: exact affine hull of forward-step polynomials

Status: PASS. This is a separate audit; neither frozen producer packet was modified.

## All-size proof checks

The cycle expansion groups each nonempty collection R with |V(R)|−|R|=2s into 2^(|R|−1) C_s. This proves integral nonnegative coordinates. The exact vanishing orders2s at1 prove independence. The actual reversed-extreme-arc family has binomial coefficient matrix with diagonal1, so triangular integer subtraction proves the reverse lattice inclusion. Negative coefficients are correctly allowed only in affine combinations, not interpreted as tournaments.

Consequently every universal additive/affine-linear identity reduces to the actual B_s family. This statement concerns the full affine lattice, not just a rational span: unimodularity is essential and present. The linear-span conclusion follows by including A_n. Taking valuations after evaluation gives the minimum of the generator list, and its equality with the minimum over actual B_s follows in both directions from the integral change of basis. This does not require a nonactual integer combination itself to realize a tournament.

The exact odd-prime valuation compiler is proved. The further digit-sum closed formula remains conjectural, despite the finite exact table. The producer distinguishes these statements correctly.

## Structural hostile, independently checked

For n>=5, alpha_1=c_3 and alpha_2=c_5+2t_(3,3). Thus alpha_2=0 forbids both directed5-cycles and pairs of disjoint triangles. A strong component with at least5 vertices contains a5-cycle, so every strong component has size at most4. Two cyclic components would give disjoint triangles, hence at most one component is cyclic. Its triangle count is at most2. Therefore alpha_2=0 implies alpha_1<=2, and A_n+3C_1 cannot be realized at any n>=5.

The source is verified directly in Moon's1966 primary paper, Theorem1, which proves the stronger statement that each vertex of a strong tournament lies on cycles of every length3..n. Theorem3's extremizer is indeed the classical reversed-extreme-arc family. The publisher PDF is:

https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf

DOI: https://doi.org/10.4153/CMB-1966-038-7

The source credits L. Moser for non-vertex-specific pancyclicity and strengthens it in Theorem1. The producer's use is valid; invention of the reversed-arc family is not claimed.

### Entirely elementary order5 proof

The n5 hostile can be rejected even without the general pancyclicity theorem. A nonstrong5vertex tournament has at most one nontrivial strong component, of size3 or4. Four vertices have at most2 cyclic triangles: total outdegree6 minimizes the sum of binomial(outdegree,2) at2, so c_3=4−that sum<=2.

If the tournament is strong, take a longest directed cycle. A4cycle can absorb its remaining vertex: if that vertex has both incoming and outgoing arcs to the cycle, a boundary edge permits insertion; otherwise it cannot be in the same strong component. If the longest cycle has length3, neither outsider can have mixed adjacency, or insertion would create a4cycle. Strong connectivity forces one outsider below the entire triangle, one above it, and an arc from the lower to the upper. Replacing any triangle edge by the three-edge path through these outsiders yields a5cycle. Therefore three cyclic triangles force a5cycle. The candidate has c_3=3,c_5=0, a contradiction.

### Recovered known H=7 exclusion

This is a short recovery of a known forbidden count, not a new spectrum theorem. If a5cycle exists, its strong5vertex induced support contains at least3triangles (Moon, Corollary2.1), giving alpha_1>=3 and alpha_2>=1. If two disjoint triangles exist, alpha_1>=2 and alpha_2>=2. Either case yields H=1+2 sum alpha_s>=9. Otherwise there is no5cycle and no disjoint pair; the preceding component structure has at most one cyclic component of size at most4, so all higher coordinates vanish and H<=5. Hence H cannot equal7. This proof explains why the affine lattice admits H7 while actual incidence excludes it. It does not prove the21 exclusion or completeness of the spectrum.

## Independent computations

The checker imports only function definitions from the independent basic auditor, not producer code. It uses direct vertex permutations and cyclotomic field-norm determinants rather than the producer's subset DP and repeated uniformizer division.

-1099 complete labelled tournament coordinate reconstructions through order5.
-20 actual B_s family comparisons through order8.
-155 order/prime valuation cells through order31, at2,3,5,7,11; generator and actual-family minima agree.
-The positive palindromic order5 hostile is absent, and its affine certificate and coordinates are independently confirmed.
-Normal and optimized output agree byte-for-byte. No assertions or float literals occur.

The producer's frozen manifest and local dependency pins are checked separately in `source_audit.json`. No mathematical correction is required. The finite computations supplement the all-size proof, rather than establish it.

## Publication portability replay

The checkers also run from an unrelated working directory. The affine audit resolves its independent basic-audit sibling relative to its own file, accepting both workspace and published directory names. The basic audit has no runtime packet dependency. The producer affine checker was repinned after the same layout-only repair; no mathematics or numerical results changed. `portability.json` records normal and optimized executions of both producer packets and both independent audit packets in a mock published directory, all launched from `/tmp`, plus workspace-layout independent-audit replays.
