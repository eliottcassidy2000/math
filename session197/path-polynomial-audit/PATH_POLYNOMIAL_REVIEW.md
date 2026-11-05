# Independent review: positive central-factorial path profiles

**PASS**, no corrections requested, 2026-11-05.

Frozen source in `session197-odd-cycle-path-polynomials/`:

- PROOF.md: b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b
- verify.py: 10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb
- Normal/optimized output: ece25a6c3b209d7262be6ca5c057fc4d39e7daecd3588c5236269db396d0b774

The source verifier was copied into a separate audit directory. Normal and optimized runs both reproduce the stated output hash exactly. No source files were edited.

## All-order proof checks

1. The squarefree-ring determinant argument has the correct signs. The matrix determinant lemma enumerates directed vertex words, and A+A^t=J-I gives the displayed factorization using Z^2=0. Taking log(det(I+ZA)/det(I-ZA)) retains only odd closed walks. Rotation multiplicity cancels the trace denominator, and repeated vertices vanish. The bin-count expansion then gives the positive odd-cycle packing formula with positive exponents of the correct parity.
2. The q/r bases are monic triangular bases of the two positive-degree parity spaces. Multiplication by x^2 has the stated positive two-term recurrence. In join products, after factoring the larger P-index, every remaining factor x^2-s^2 has current index at least s, so its coefficient (k+1)^2-s^2 is nonnegative. The parity prefactors in all three product cases are correct.
3. D has only three basis terms: the root vanishing argument gives divisibility by q_(k-1) or r_(k-1), with the necessary double zero at zero in the even case. Leading and next-leading coefficients give the level term, and evaluation at x=k gives the displayed downward terms. The k=0 endpoints agree separately.
4. The cyclic wrapper D is inherited correctly. The cited S193 multi-port theorem was independently retrieved at its exact pinned path. Applying it to a transitive two-vertex core with both vertices as ports gives P_A(x)=x^2, sum L(x)=2x-1, sum R(x+1)=2x+1, and sum K(x+1)=2. Hence the claimed physical E operator holds for every nonempty input. The identity E=DX+XD-X^3-2X is exact; independently collecting its q/r terms reproduces all three/four positive bands and their boundary coefficients.
5. The H readout coefficients are correct: (3,6)/(5,12) for D and (5,18)/(11,36) for E. These prove the stated mod-12 invariant for the smaller grammar and the H13/H19 exclusions for the enlarged grammar. Prime descent through an ordered join is legitimate. The exclusions concern only these expression grammars, not all tournaments or all positive coordinate states.
6. The Motzkin first-return argument gives a Jacobi continued fraction, with downward coefficient b_k rather than its square. Diagonal symmetrization gives neighboring entries sqrt(b_k). The finite Gram determinant argument proves the vacuum Hankel product; no infinite analytic convergence or arbitrary-initial-state moment claim is being assumed.
7. The imaginary-root threshold is proved rather than plotted. The m=6 and m=7 discriminants are negative. For even m>=8, Newton's first coefficient inequality reduces to m<=7; for odd m>=9 it reduces to m^2-7m-2<=0. The small m=1,...,5 reduced polynomials have real negative roots, using the displayed discriminants and positive coefficients. Zero roots are handled by the factored x or x^2.
8. The odd-cycle packing hostile has the correct parity extraction. Any exterior cycle uses both added vertices and an odd-sized internal path. A triangle block supplies 1,3,3,3 path choices by selected size, while an unused block may contain its own packed triangle. This gives exactly the displayed formula and does not interpret its formal subtraction as a physical operation.

## Independent graph and algebra controls

The new `path-polynomial-audit/independent.py` does not import the source checker.

- It counts Hamilton paths by literal DFS without subset DP, verifying the H11 row witness and both negative-Hankel witnesses.
- It enumerates unordered set partitions and literal path orders inside each block to recover the H11 graph's full cover vector, then checks its labelled-bin polynomial against E(x^2) at twelve integer arguments.
- It enumerates actual directed odd-cycle vertex words by DFS, obtaining 507 cycles in the eleven-vertex hostile, then independently enumerates all vertex-disjoint cycle collections. The result is exactly 1+507z+147z^2+19z^3, and literal Hamilton DFS gives 1755.
- A separately implemented fraction-free Bareiss resultant calculation gives discriminants 7092, -8752, and -224218032 for the m=5,6,7 reduced polynomials.
- Literal physical continuations give (1,17,241) and (1,33,645), with first Hankel determinants -48 and -444.

The supplied normal/-O replay additionally checks all 1,099 labelled profiles through order five, 75 literal-cover controls, 150 full physical wrapper identities, 84 operator basis controls, 196 positive basis products, and the stated vacuum moment checks. The independent controls and the proof checks distinguish exact all-order identities from finite examples.

## Scope

This packet supplies a positive unbounded profile representation and two actual physical operators. It does not prove that arbitrary positive vectors are physical, that arbitrary physical initial states are moment sequences, that rooted/color interfaces close in these coordinates, or that the enlarged grammar reaches all admissible Hamilton counts. The hostile examples appropriately prevent those extensions.
