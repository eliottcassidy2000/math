# Independent continuation audit: contextual collision and squarefree ports

2026-10-23. Verdict: **PASS** for the six-vertex collision, minimal order, all-word congruence, adjacency-spectrum identification, squarefree path/port identities, and exact odd-cycle explanation of the difference. This is separate from the frozen initial audit. No Git operations.

Audited sources: `session35-root/scan_marked_spectra.py`, `session35-root/marked_spectrum_scan.json`, `session35-recursions/CONTEXT-AND-PORTS.md`, and the producer's collision script/output. The independently reproduced range stops at six vertices. The n=7 scan's algorithm and integer bounds were reviewed, but its reported168/47 class totals were NOT independently rerun.

## 1. Targeted collision, genuinely independent counting

The two outneighbor lists are valid tournaments:

- T=[58,4,9,18,38,14].
- U=[22,12,8,17,38,15].

Independent exhaustive permutations of all six vertices give H(T)=29 and H(U)=25. This uses no subset/end-vertex dynamic programming. Their characteristic polynomial was independently reconstructed by interpolation of seven exact numerical determinants. Their marked numerator was reconstructed from six exact numerical linear solves and the characteristic polynomial; it uses no adjugate or matrix inverse. Results:

p=z^6+15z^4+47z^2+1,
q=2z(3z^4+34z^2+59).

The degree bounds make those interpolation checks exact polynomial identities rather than agreement at arbitrary sample points: deg(p)<=6 and deg(q)<=5. All solve points are positive real, so real skew symmetry guarantees invertibility.

Independent direct permutation counts for every induced vertex subset reproduce coefficient-difference totals0,0,8,9,3,1 at degrees1,...,6. On labels{0,1,2}, the counts are3 and1.

## 2. Census range, minimum order, and signature completeness

The independent census enumerated ALL labelled tournaments of orders1 through5 (counts1,2,8,64,1024), not merely those with a designated path. It then enumerated all1024 designated-path presentations at order6. Direct-permutation Hamiltonian counting and Python arbitrary-precision matrix products give:

- Orders1,...,6: marked-spectral class counts1,1,2,3,9,28.
- Orders1,...,6: H-ambiguous class counts0,0,0,0,0,3.

Thus six is indeed the smallest positive order with a same-(p,F), different-H pair. This minimal-order certification does not depend on the unrerun order7 scan.

For the producer's designated-path scan, every adjacent edge i→i+1 is fixed, while the other (n-1)(n-2)/2 edges are independently enumerated. Consequently the stated presentation counts2,8,64,1024,32768 for n=3,...,7 are correct. These are NOT isomorphism-class counts. Every tournament can be relabelled to contain the designated directed Hamiltonian path: an elementary inductive proof inserts a new vertex immediately before the first old path vertex it beats, or at the end if it beats none. The preceding vertex, if present, beats the new vertex. This justifies coverage without assuming random sampling.

The key consists of Tr(M^k) and m_k=1^T M^k1 for k=1,...,n. Newton identities recover the monic characteristic polynomial from the traces. The Laurent series F(z)=sum_(k>=0) m_k z^(-k-1), with known m_0=n, determines q=pF using only m_0,...,m_(n-1). Cayley–Hamilton determines every later moment. The retained m_n is redundant but harmless. Hence equal keys imply exactly equal(p,F), not merely a low-moment approximation; conversely equal(p,F) and p imply the same keys.

The int64 scan is safe at these sizes. For n<=7,k<=7, a crude entry bound is n^(k-1) and the all-entry absolute sum is bounded by n^(k+1)<=7^8=5,764,801. Intermediate matrix-product summands and their accumulation remain under the same elementary bounds, far below2^63. The independent census uses Python integers anyway.

## 3. All-word theorem and what it does NOT say

Each E or D maps the triple(n,p,F) deterministically to the next triple. Induction proves equal initial triples remain equal after EVERY finite word, not just bounded-length experiments. The p polynomial already fixes n, but retaining n explicitly is harmless.

E preserves H because its universal source must be first and its universal sink last in every spanning directed path. Hence the witnesses retain counts29 and25 after E^k. There is no supported assertion that the two H counts differ after every word involving D, and the addendum expressly avoids that assertion. The empty context alone already supplies nonrecoverability of H; unlimited contextual spectral measurements do not remove the ambiguity because the entire contextual family is equal.

## 4. Identification with ordinary adjacency cospectrality

The parity recovery is correct. Since M+J=2A+I,

g(z)=det(zI-M-J)=2^n p_A((z-1)/2)=p(z)-q(z).

Skew symmetry gives parity(-1)^n for p and the opposite parity for q. Taking the two parity projections recovers both p and q exactly, including cancellations in F=q/p. Thus the marked pair is equivalent to ordinary adjacency spectrum for same-size tournaments. Calling this a wholly new carrier stronger than adjacency spectrum would be incorrect; the addendum explicitly recognizes the identification and known predecessor.

## 5. Squarefree path lift and Schur ports

In R=Q[x_i]/(x_i^2), all positive-degree elements are nilpotent, so AX has a finite geometric inverse. The entry of X(AX)^r includes variables for the starting vertex and every subsequent visited vertex. Repeated visits vanish; every simple directed path contributes with positive coefficient. The empty-path convention belongs only to Phi=1+1^T P1, not to the port matrix itself.

The E factorization is exact for every subset, including the no-old-vertex cases. The D adjacency block is correctly [[A,A+I],[A,A^T]]; the off-diagonal identity corresponds to the orientation of the two copies of each vertex. The displayed Schur complement has the correct multiplication order and signs. Standard block inverse formulas are valid over this commutative ring because the relevant matrices have identity constant term.

Independent ring code checks all FOUR inverse blocks and EVERY port coefficient against direct simple-path enumeration for doubled two-vertex, transitive three-vertex, and cyclic three-vertex inputs. It also compares every coefficient of the scalar path polynomial with a separately constructed product over odd cycles. Twin variables are kept independent throughout.

Identifying twin variables would annihilate the full spanning coefficient since it visits both copies. The full port matrix recovers labelled adjacency at degree2; the scalar polynomial does not, since reversal preserves its path counts. The exponential-cost/nonminimality disclaimers are correct.

## 6. Odd-cycle expansion and exact difference4

The determinant lemma gives Phi=det(I+(J-A)X)/det(I-AX). For a tournament J-A=I+A^T. The factorization (I+A^T X)(I+X)=I+(I+A^T)X is valid because X^2=0. Transposition and det(I+UV)=det(I+VU) yield det(I+A^T X)=det(I+AX).

Finite logarithms over Q then give singleton terms plus twice the odd trace terms. No even trace survives; loops are absent; each surviving cycle of length k occurs in k rooted forms. Exponentiating selects vertex-disjoint odd cycles with weight2 each and optional singleton vertices. The claimed OCF follows directly, with no analytic convergence requirement.

Independent cycle enumeration gives c3=c5=6 for both witnesses. T has one disjoint triangle pair; U has none. At six vertices the only nonsingle odd-cycle collections are two triangles. Therefore H=1+2(c3+c5)+4D33 yields29 and25, and the whole discrepancy is exactly4. It is compatibility of cycle placements, not the individual odd-cycle totals, that is missing.

## 7. Reproduction and certification limits

Independent deliverables:
- `addendum_check.py` and `addendum_check.out`: full labelled census through5, designated-path census6, direct permutation witness/subset/cycle counts, determinant/solve interpolation.
- `ports_check.py` and `ports_check.out`: exact square-zero-ring Schur blocks, all port coefficients, and full odd-cycle polynomial controls.
- Corresponding `.optimized.out` files: python -O outputs compared byte for byte. All failure conditions raise explicit exceptions.

These finite controls supplement the algebraic proofs; they are not an all-order experimental certification. No order7 replay, independent literature-priority determination, or certification of inherited repository theorem labels was performed. The note's direct proofs are sufficient for the identities audited here.
