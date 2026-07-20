        # Message: mac-mini-2026-07-20-S130: THM-1440 the two 'odd n only' facts are ONE parity fact and it is the parity of sin/cos + the skew-Seidel spectrum is a COMPLETE switching invariant iff n<=6, first failing at n=7 with one cospectral pair p(x)=x(x^2+7)(x^4+14x^2+17)

        **From:** mac-mini-2026-07-20-S?
        **To:** all
        **Sent:** 2026-07-20 09:17

        ---

        OWNER: 'think seidel spectra and even graph bijections that are odd n only, and how these two interlaced can be seen in sin and cosine.'

THE ANSWER: they are the SAME parity fact, and the coat it wears is literally the parity of sin and cos.

(A) THE CHARACTERISTIC POLYNOMIAL HAS THE PARITY OF sin/cos. S = A - A^T is skew, so p(-x) = det(-xI-S) = (-1)^n det(xI+S) = (-1)^n det(xI-S^T) = (-1)^n p(x). So p is an EVEN polynomial at even n (like cos) and an ODD polynomial at odd n (like sin) -- and that is exactly why odd n FORCES a zero eigenvalue: sin(0)=0, cos(0)!=0. Concretely p = prod(x^2+lambda_k^2) at even n and x*prod(x^2+lambda_k^2) at odd n, the same extra linear factor separating sin z = z*prod(1-z^2/k^2pi^2) from cos z = prod(1-z^2/(k-1/2)^2pi^2). Verified over all iso classes n=3..7: eigenvalues purely imaginary, parity exact, zero present at every odd n, rank always even (skew rank is even, so odd n forces corank >= 1).

(B) THE NEW CONTENT. Switching sends S -> DSD with D = diag(+-1) and D = D^{-1}, so it is a SIMILARITY and the spectrum cannot move. Sharp contrast with S129, where H, min-FAS and the cyclic-triangle count are NONE of them switching-invariant. Then:
  n            3   4   5   6    7
  distinct spectra   1   2   2   6   11
  A049313            1   2   2   6   12
THE SKEW-SEIDEL SPECTRUM IS A COMPLETE INVARIANT OF THE SWITCHING CLASS (= oriented two-graph) FOR n <= 6, AND FIRST FAILS AT n = 7 WITH EXACTLY ONE COSPECTRAL PAIR.

THE PAIR, EXHIBITED. Independent recomputation gives 456 iso classes, **12 switching classes -- an independent verification of A049313(7)** -- and 11 spectra. The two cospectral classes contain EXACTLY 64 iso classes each and share
  spec(iS) = {0, +-sqrt(7), +-sqrt(7-4sqrt2), +-sqrt(7+4sqrt2)},  p(x) = x(x^2+7)(x^4+14x^2+17)
verified to the integer (coefficients 1,0,21,0,115,0,119,0 -- supported on degrees 7,5,3,1, ALL ODD, which is (A) at work). Splitting field Q(sqrt 2); sum lambda^2 = 42 = n(n-1) as forced by tr S^2 = -n(n-1). sqrt(7) appears but the pair is NOT Paley -- Paley on 7 vertices is the degenerate {0, +-i sqrt7} with multiplicity 3, whereas all three magnitudes here are distinct. Representatives have scores [0,1,3,4,4,4,5] and [0,2,2,4,4,4,5] with 4 and 5 cyclic triangles, but NEITHER separates the classes, nor does H (value sets overlap heavily) -- none of them is switching-invariant.

(C) WHY THE EVEN-GRAPH PROJECTION IS ODD-n ONLY -- PROVED TWICE.
  C1. For T_cycle = (I+L(K_n))T mod 2 the image at e={u,v} is T_e + d_u + d_v, so the image degree at w is (n-1)*d_w + sum_x d_x = (n-1)*d_w (mod 2), since sum_x d_x = sum_f 2T_f = 0. THE IMAGE IS AN EVEN GRAPH IFF n IS ODD. Verified n=3..8.
  C2. A cut delta(S) has all degrees even iff |S| and |S^c| are both even -- impossible at odd n. So the BICYCLE space Cut & Cycle of K_n is 0 iff n is odd, and has dimension n-2 when n is even (verified n=3..9: 0,2,0,4,0,6,0). That is exactly when Cut (+) Cycle is a genuine direct sum.
THIS COMPLETES THM-1405, whose gauge-bits/holonomy-bits splitting I flagged as a non-canonical choice of complement precisely because of the bicycle space: the splitting is CANONICAL EXACTLY AT ODD n (HYP-8295).

(D) CIRCULANTS MAKE THE SINE LITERAL. For connection set C with C u (-C) = Z_n\{0} disjointly,
  mu_j = sum_{k in C} (w^{jk} - w^{-jk}) = 2i * sum_{k in C} sin(2 pi j k / n)   -- a pure SINE sum,
and z_j = sum_{k in C} w^{jk} obeys z_j + conj(z_j) = -1, so Re(z_j) = -1/2 for EVERY j != 0: the cosine part is PINNED to a critical line and all the content of a circulant tournament sits in its sine part. Verified < 2e-14 over all 4/8/16/32 circulant tournaments at n=5,7,9,11. Paley: spectrum exactly {0} u {+-i sqrt q}, verified q=3,7,11,19,23.

(E) INTERLACING IS THE sin/cos ZERO INTERLACING. iS is Hermitian, so Cauchy interlacing holds under vertex deletion -- verified on all 16/60/336/3192 deletions at n=4..7. Odd n: spectrum symmetric about 0 and CONTAINS 0, like the zeros of sin (0, +-pi, ...). Even n: symmetric and OMITS 0, like the zeros of cos (+-pi/2, ...). Deleting a vertex flips the parity and the spectra must interlace -- precisely the classical fact that the zeros of sin and cos interlace.

(F) THE THREE-TIER PICTURE, with THM-1420 from S129:
  any F_2-linear functional : iso-invariant NO  (THM-1420 -- none exist at all)
  H, min-FAS, triangles     : iso-invariant YES, switching-invariant NO
  skew-Seidel spectrum      : iso-invariant YES, switching-invariant YES
Strictly nested, spectrum coarsest -- which is exactly why it can be a COMPLETE invariant of switching classes at small n while H is not even defined on them.

PRIOR ART CHECKED AND CREDITED, and I want this on the record given last session's d_sat rediscovery: THM-474 (tilings ARE switching classes) is the framework; the determinant-lens reflection ALREADY carries the odd/even-function axis (graph = even/symmetric, tournament = odd/skew, S = A - A^T) that the whole sin/cos reading is built on -- THE FRAMING IS NOT CLAIMED AS NEW. (A) and (D) are classical facts about skew matrices, circulants and Paley tournaments, re-derived and verified, not claimed. THM-1430 is the undirected companion (graph switching = E_n = A002854). ONLY (B) and (C)'s bicycle identification are claimed new.

HANDOFF -- three:
(i) HYP-8290. The n=7 completeness boundary coincides with THM-499/500's ADJACENCY-spectrum boundary (H non-spectral from n=6, odd-cycle count from n=7). But those use a DIFFERENT matrix: S = 2A + I - J, and only the S-spectrum is switching-invariant. ONE MECHANISM OR TWO? Nothing here connects them. Next test: n=8 against A049313(8) = 79 -- if the cospectral-pair count grows like the adjacency defect, that is evidence for a shared cause.
(ii) HYP-8295. At even n there is an (n-2)-dimensional bicycle space with no canonical complement. Nobody has asked what it IS combinatorially, or whether a natural even-n substitute for the projection exists.
(iii) The cospectral pair's splitting field is Q(sqrt 2), and sqrt 2 is the repo's hypotenuse/leg ratio. Suggestive, unexamined, and I am deliberately NOT claiming it means anything.

Artifacts: THM-1440; 04-computation/seidel_spectra_sine_cosine_macmini_S130.py, seidel_cospectral_pair_n7_macmini_S130.py (+outs); HYP-8290, HYP-8295.

        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
