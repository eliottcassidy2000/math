# Independent audit: prime-plane unique-sum frames

## Verdict and exact domain

The argument is correct for an odd prime p, q=p^n with n>=2, a fixed two-dimensional F_p-linear subspace U of F_q, and A=C+U. A frame consists of two unordered p-element subsets B,D of F_q with the addition map B×D→A bijective. The total number is C_p q; the disjoint-factor count is C_p(q-p²), where C_p=(p+1)p^(p−1)−p(p+1)/2. Do not call the latter the count of all frames without specifying disjointness. The factorization theorem also holds at p=2, but the odd-characteristic disjointness formula does not.

## Proof checks

1. Choose b0 in B and d0 in D. Then B−b0 and D−d0 are subsets of U and sum uniquely to U, since b0+d0 belongs to A. This normalization requires no assumption that the original factors lie inside A.
2. For each nonzero linear form ell:U→F_p and primitive pth root zeta, the convolution identity gives (sum_b zeta^ell(b))(sum_d zeta^ell(d))=0. For either p-element factor E, write its character sum as P(zeta), with integer coefficients n_j=#(ell^-1(j)∩E). Since deg P≤p−1 and Phi_p=1+X+...+X^(p−1) is the minimal polynomial, vanishing forces P=c Phi_p; evaluating at 1 gives c=1. Thus ell is bijective on E, exactly saying that E misses direction ker ell.
3. Every direction is therefore absent from at least one factor. The two determined-direction sets are disjoint. If both factors were nonlinear, their combined size would be ≥p+3, exceeding the p+1 available directions. Thus at least one factor is an affine F_p-line.
4. For each of p+1 directions L, there are q/p choices for a line B of that direction in the ambient group. D must be contained in the uniquely prescribed coset of U and be a transversal of its p cosets of L. There are p^p such D. The line-first ordered incidence count is (p+1)q p^(p−1).
5. Both-line ordered incidences number (p+1)(q/p)·p·p=p(p+1)q. All mixed line/nonline frames occur once in the incidence count; both-line unordered frames occur twice. Subtract half the latter ordered count, giving C_p q. Equal factors cannot be frames since swapping two distinct elements gives repeated sums.
6. Orient mixed frames line first and both-line frames by a fixed direction ordering. This is preserved by the gauge action (B,D)→(B+t,D−t), so it defines an action on canonically oriented unordered frames. A stabilizer t lies in the line direction L, while the transversal D has no nonzero period in L. Thus each orbit has q elements, giving C_p orbits.
7. The difference map B×D→B−D is injective: b−d=b'−d' implies b+d'=b'+d, and unique sums force b=b', d=d'. Therefore |B−D|=p². A gauge is bad exactly when 2t=d−b for some b,d. Since 2 is invertible, exactly p² distinct gauges are bad, completing the disjoint count.

## Primary-source verification

Gábor Somlai, “A new proof of Rédei’s theorem on the number of directions,” Archiv der Mathematik 122 (2024), 575–580, Theorem 1.1, DOI 10.1007/s00013-024-01979-x. Publisher page: https://link.springer.com/article/10.1007/s00013-024-01979-x . Independently checked on 2026-10-25: its theorem explicitly applies to a p-point subset of F_p² not lying on a line and gives the lower bound (p+3)/2. The author gives a new proof; this is a primary research source, not merely a third-party summary. Preprint: https://arxiv.org/abs/2212.12823 .

## Important caveats

- Nonzero difference sets being disjoint does not, by itself, allow replacing vectors by projective directions. For example in F_5², {0,(1,0)} and {0,(2,0)} have disjoint nonzero difference vectors but the same direction. The p-element full-plane hypothesis plus the Fourier argument is what justifies directional disjointness here.
- A p-element affine F_p-line is not necessarily a three-term arithmetic progression, nor an affine F_q-line. For p>3 use “line” rather than “AP” unless AP is explicitly defined as the full cyclic progression.
- This is a classification/count consequence of the classical prime-plane factorization and direction theorem, not a claim that the at-least-one-line theorem is new.
- At p=2 the equation 2t=d−b cannot be inverted; the odd-characteristic count must not be extended there.

## Reproducible computational audit

`audit.py` uses explicit RuntimeError checks, not assert statements. It independently enumerates all p-element plane subsets for p=3,5, testing the direction bound and the exact cyclotomic coefficient criterion for each direction. It also constructs every line-first plane frame for p=3,5, checks sums and differences, deduplicates unordered frames, constructs canonical gauge orbits, and counts bad gauges after embedding each representative into F_p³. Run under both normal Python and python -O; outputs must match exactly. Numerical results are recorded in normal.json and optimized.json.

## Review of the complete construction packet and Hamilton paths

Reviewed `/tmp/session57-prime-plane-frames/prime_plane_frames.md`, including its gauge-equivalence formulation and section 8. All statements pass. Its intrinsic gauge-equivalence description is correct without choosing a canonical orientation. The generated group is a semidirect product with swap conjugating T_t to T_(−t); no swapped translation stabilizes an ordered frame because the two factors cannot be translates.

The complete bipartite skeleton has (p!)² unoriented Hamilton paths: the two initial-side choices are exactly canceled by reversal. The alternating bipartition of any such path uniquely recovers the unordered factor pair, so paths do not duplicate across frames. The exact line/line fraction simplifies to 1/(2p^(p−2)−1), and its conditional disjoint version is undefined at q=p², correctly guarded in the packet. The two-path obstruction follows from 4p−2<p² for p≥5; it concerns unique labels, as claimed.

`audit_paths.py` separately checks the Hamilton-path count and recovered bipartition for p=3,5, the fraction simplification, and p=5 size obstruction. Normal and optimized outputs match. The main audit normal/-O outputs also match byte-for-byte. At p=5 the independent line-first construction yields 93,375 distinct unordered plane frames, 375 both-line double-counts, 3,735 gauge classes, and exactly 25 bad gauges in every class embedded in F_5³. Completeness of the construction for p=5 uses the proven classification; the script does not exhaust every possible pair of five-point subsets.

The frozen packet's final section 9 also passes: affine bijections of a fixed F_p-line correspond exactly to AGL(1,p), order p(p−1), with index p!/[p(p−1)]=(p−2)!. Its two-side ordering counts are 2(p!)² total and 2[p(p−1)]² affine; for p=3 both equal 72, and for p=5 the one-side fraction is 1/6. Equality occurs only at p=2,3. No real-order interpretation is imported.
