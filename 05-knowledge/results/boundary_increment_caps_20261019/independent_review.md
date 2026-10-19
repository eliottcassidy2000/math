# Independent review: boundary-increment caps

Date: 2026-10-19. Verdict: PASS; no blocking mathematical or implementation issue found. Scope is the stated fixed full root core and fixed projected order, not an entire morphology or the seven-thin conjecture.

## Mathematical checks

- For a fixed arbitrary linear extension of A, deleting q from an extra interleaving of wq leaves an extra interleaving of w when w is nonempty: an A letter following q still follows w. For empty w, the image is the unique fixed A extension, consistent with U_empty=1.
- If a_j is the last A predecessor of q, every later A letter is incomparable with q, using the stated absence of q<A relations. There are m-j eligible A gaps before the final A letter, and at most one insertion position per such gap after all of w. Thus each deletion fiber has size at most m-j <= s_A(q). Restrictions from the particular order of A or w only shrink the fiber. Summing over A extensions and dividing by e(A) proves the cap. Empty A gives zero nonempty increments.
- Order reversal gives the suffix inequality with the first suffix letter removed. Multiplication by the other nonnegative increment and relaxation to k-d_Q(q) are valid; no joint slack allocation is needed for these weaker necessary inequalities.
- Every supported child has a supported parent because its witnessing core extension still begins/ends with the shortened boundary words. Unsupported cells contribute zero core-extension counts, so their coefficients need not vanish and can be omitted.
- Telescoping boundary sums yield the supported generator expansion of the ambient distribution. With C lambda <= 0 and nonnegative t, componentwise W-sum(t C) <= -1 implies W lambda <= sum(t C lambda)-sum(lambda) < 0, since the supported empty/empty coefficient equals 1. Positive comparison weights contradict strictly positive forward 2/3 margins. The sign is correct.

## Computational checks

Read the complete proof, JSON certificate, and verifier. Independently enumerated extensions by recursive choice of available minimal vertices (rather than the producer verifier's permutation filter), formed weighted comparison margins directly as +weight / -2 weight, built supported prefix/suffix columns, and subtracted the 18 cap rows. Obtained exactly 408 extensions, 133 columns, and residual range [-540,-1]. The producer verifier also explicitly validates degree slacks, forward incomparability, supported parents, positive weights, and the fixed core/order. No broad scan was rerun.

## Reviewed SHA-256 hashes

- boundary_increment_caps.md: 6dff1b18e21e40cf701686b8263ca25206bd36decdc3d6da583d38997ad37b89
- fixed_certificate.json: 5f0d5b37e2209d07120a77d303aa29e04c347748de6358ed19cf9af6cedb618f
- verify_caps.py: 5ddd3cbfdfc0a49dcf83490a33d26df3d1e60f767e46f9bf286ae81a7df42728
