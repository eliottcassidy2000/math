# Independent audit of the Mahler reset-cost and native-clock theorem

**Verdict: PASS, no mathematical repair requested.** The cost spectrum, native-height conversion, sharp rejecting family, finite-cycle lag law, and stated fair-coin limit are valid with the declared hypotheses. The real companion and the Rule30 spatial inverse observer are correctly prevented from substituting for the ordinary infinite carry carrier.

Audited producer: [proof](continuing8_20260906_mahler_reset_cost.md), [source](../../04-computation/continuing8_20260906_mahler_reset_cost.py), [certificate](continuing8_20260906_mahler_reset_cost_certificate.json). The canonical affine-cylinder and strict follower statements were read in THM-4072; the Rule30 generator/reset convention is that of THM-4204. No producer implementation is imported by the referee.

## Analytic audit

For a nonzero strict-safe prefix ending at follower zero, the last carry one at index \(j-1\) leaves the oriented finite phase \(D/3^j\) reduced: its numerator is \(2^{j-1}\) modulo three. Strict suffix safety gives \(0<D/3^j<1/2\). Appending zero carries changes neither reduced numerator nor denominator.

The literal ceiling iterate satisfies \(2^nT^n(A)=3^nA+C_n\). Thus the compatible launch is one residue class modulo \(2^n\), and its representative is
\[
A_n=(2^n b_n-D)/N,\qquad b_n=[2^{-n}D]_N,\qquad N=3^j.
\]
It is positive since a zero launch has only zero carries. Every positive ordinary completion is \(A_n+2^n k\), \(k\ge0\); its exit changes by \(3^n k\). The bound \(2^n\le NX+D\) for a compatible launch at most \(X\) follows without dropping the numerator shift. At \(b_n=1\) it is an equality for \(X=A_n\).

The unit-cycle proof is complete: \(v_3(4^{3^a}-1)=a+1\) follows by the binomial expansion and induction, and multiplication of the exponent by an integer prime to three does not change that valuation. Every order of two modulo \(3^j\) must be even, so it is exactly \(2\cdot3^{j-1}\), the number of units. Hence all unit residues occur, not just a cyclic subgroup selected by an unverified primitive-root assumption. The periodic object is \((A_n+D/N)/2^n\); the unshifted cost only has the corresponding subsequential limits.

For the native clock, write \(b_n=2^r c\), \(c\) odd, and \(H=\operatorname{bitlength}(A_n)\). The strict inequalities
\[
N2^{H-1}<NA_n+D=2^{n+r}c<N2^H
\]
give
\[
n+r-H=\left\lfloor\log_2(N/c)\right\rfloor.
\]
The strict upper inequality uses \(D<N\), not an asymptotic approximation. The ratio \(N/c\) cannot be a power of two. This proves the exact integer identity and its ceiling \(R=\lfloor\log_2N\rfloor\). Only when \(H\ge m\) does the interval from the actual native termination to the next odd state lie entirely in the appended-zero part at follower zero. That condition is retained in the report. At a terminal depth \(H\), the ordinary launch is its own least representative; applying the identity with \(n=H\) proves the claimed all-orbit bound within this fixed-prefix class.

For the phase distribution, each odd unit \(c\) with \(N/2^{k+1}<c<N/2^k\) admits exactly the \(k+1\) unit residues \(c,2c,\ldots,2^kc<N\), all with lag \(k\). This proves the frequency formula. In particular there are exactly \(R+1\) maximal-lag phases. The arithmetic gaps are real: at \(N=9\), only lags zero and three occur, with frequencies two and four.

The count of integers in the two residue classes \(1,5\bmod6\) differs from one third of an interval's length by at most two. This yields the stated pointwise discrepancy \(3(k+1)/N\). Summing through \(R\) and adding the exact limiting tail \((R+3)/2^{R+2}\) gives the reported total-variation bound \((R+2)^2/N\), which tends to zero. The probability space is uniform phase on a finite unit cycle of **different ordinary launches**, at sufficiently large depths. The limit \((k+1)/2^{k+2}\) is indeed the number-of-failures law before the second head of independent fair coins. It proves neither independence of the native bits nor a distributional assertion along one selected infinite ordinary orbit.

For the sharp family, \(E=2\cdot3^{j-1}t\), \(t\ge2\), makes the launch integral. Since \(E\ge R+3\), the declared height is exactly \(H=E+j-1-R\), and \(H\ge j+2\). Direct iteration gives the entire word
\[
0^{j-1}1\,0^{E-1}11.
\]
The first one is followed by the safe reset loop \(100\), every further zero is a reset, and the final second one is the first rejected symbol. At depth \(H\) there are precisely \(R\) zero resets left; rejection occurs at index \(H+R+1\). This proves universal sharpness and universal eventual rejection for this family. The hypothesis \(t\ge2\) matters to the reset interpretation: for \(j=t=1\), the terminal follower is at state one, so the next zero is a match rather than a reset.

The real companion \(A_n+D/N\) has strict fractional phases equal to the oriented remaining finite tails through time \(n\), followed by \(r\) integral phases and then an odd half-integer at time \(n+r+1\). This is its first strict real failure. At that same endpoint the ordinary carry word has only just read its first post-runway one from follower zero, which is accepted. The producer explicitly retains this mismatch and uses the separately proved final \(11\) only for the sharp family.

The Rule30 quotient obstruction is also exact. Independent integer matrix multiplication gives the same constant observer for \(00100\) and \(0010011\), while exact suffix phases give opposite finite safety predicates. The stronger launches \(180\) and \(148\) have the declared common depth, bit length, follower state and saturated observer, but exit states \(4617\) and \(3798\), hence opposite next reset decisions. The arithmetic denominators \(27\) and \(729\) are correctly retained as lost data. This refutes this specific quotient map; it is not a general impossibility of connecting Mahler and Rule30.

## Independent arithmetic path

The [referee source](../../04-computation/continuing8_20260906_mahler_reset_cost_audit.py) constructs addresses by choosing one native bit at a time. It maintains \(A<2^i\) and \(U=T^i(A)\), then adds either zero or \(2^i\) to the launch so the parity of \(U\), shifted by the odd number \(3^i\), matches the desired carry. It never calls the producer or uses its inverse-power implementation for the address.

The complete native-first universe contains **4,094 launch/depth pairs** through depth eleven. Literal ceiling iterations establish the full level bijection, oriented affine identity, and agreement of rational suffix safety with the follower; 388 nonzero safe rows additionally check the exact denominator, cost and native-clock formulas.

The referee independently reconstructs every saved unit-cycle entry and lag histogram for all seven distinct declared prefixes. It checks all-completion affine controls, direct postterminal runs, direct odd-part interval counts, and fresh complete unit distributions through denominator \(3^8\). Twenty-one fresh sharp-family rows use \(1\le j\le7\) and \(t=2,3,5\), including values not in the producer's bank. The real endpoint is checked by repeated exact rational multiplication, and the inverse observer by four-by-four integer matrices. The \(t=1\) reset-state hostile is explicit.

All **19,102 always-active exact gates** pass. [Normal output](continuing8_20260906_mahler_reset_cost_audit.out) and [optimized output](continuing8_20260906_mahler_reset_cost_audit_optimized.out) were captured as raw subprocess bytes and are identical, with no carriage returns. The certificate is read only and its full source/data hashes are pinned.

~~~text
python continuing8_20260906_mahler_reset_cost_audit.py
python -O continuing8_20260906_mahler_reset_cost_audit.py
~~~

Referee source SHA-256: ecf248959988904fbe0851887bbffec4f7b64fe0df5be16070de5d4fbabfe9ae.

Referee output SHA-256: f8ae3455f58452d90d95ca6661e8b64ff0dd748c75e23d722fe00c0b82737d5c.

Pinned producer source: 77757ed311aac5246cf3a420dea282e84cb7aaaf428abcc928c633e8865f5646; certificate: ac3f7f74a5ed80351218669c628a4bf068384f7e8132fa41af0b0828b6646d6c. The full producer report with source/output pins was read; no change is requested.
