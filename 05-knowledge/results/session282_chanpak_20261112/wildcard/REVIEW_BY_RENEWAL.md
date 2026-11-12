# Independent algebraic review of the rational-carry lane

Reviewer: Session282 renewal lane. **PASS.** The standalone verifier completed successfully after recovery; its complete output is retained in REVIEW_REPLAY.json.

I checked the two telescoping carry identities directly from the weighted U,V sums. Odd h is not used for rational carries: the relevant hypothesis is gcd(p,q)=1 after dividing h,Q by g. The bound |x|<g follows from the physical digit range, and is strict at every finite depth. The forward/reverse orientation in equations(3)–(5) agrees with the word-matrix convention.

Expanding the triangular matrix gives exactly the direct, arbitrary-duration Q-channel and adjacent signed-channel terms in(6). Normalizing the Q-channel by h^(k−2) leaves powers(Q/h)^(j−i−1). Substituting delta c_j=q x_j−p x_(j−1) and summing makes all interior terms cancel, leaving−p x_i because terminal x_k=0. The reversed a calculation likewise leaves−p y_(j−1). Thus the local action(7) has the correct factor p and signs.

The h^(k−2) divisibility is valid for k≥2; the separate k1 interpretation avoids a false integer-power assertion. The general varying-letter valuation budget follows because a matrix path omits at most two diagonal factors, each other diagonal being h_i or Q_i (or zero, which only removes a term). Subtracting the two largest valuations is therefore a safe bound, not an assumed equality.

The finite weighted-automaton claim retains both carries, the preceding letter pair and endpoint flags. Its Laurent variable tracks normalized difference, not absolute H. This scope is necessary and is retained. The central E14 identity scales under arbitrary contexts by h(L)h(R), including singleton contexts; it does not give a free operation on unrelated roots.

No substantive correction requested. The h169 numerical pair and reproduction controls are checked by the lane's standalone verifier; this review does not replace independent auditing of the inherited fusion theorem.
