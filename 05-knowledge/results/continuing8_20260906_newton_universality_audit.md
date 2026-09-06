# Independent referee of exact Newton circuit-ratio realization

**Full analytic proof, source, and producer prose: PASS. No mathematical repair requested.** This referee establishes the precise theorem below independently from the declared construction and frozen producer source. The complete producer report has now been read, including the strict integer-root alternative and both moment-sidecar exclusions.

The producer is continuing8_20260906_newton_universality.py, in the adjacent continuing8_20260906_root folder outside the repository. The referee imports no producer implementation. The coefficient construction is independently solved as a triangular product, and exact root counts use Sturm sequences in addition to a direct proof of strict term dominance.

## Exact scope and independent proof

For every integer \(d\ge2\) and every positive real vector \((c_2,\ldots,c_{d-1})\), there is a monic degree-\(d\) polynomial
\[
N(n)=\prod_{i=1}^d(n+r_i)
\]
with distinct \(r_i>0\), \(\sum_i r_i=d\), whose normalized Newton ratios \(R_k\) satisfy \(R_k/R_{k-1}=c_k\) exactly. Empty vectors at \(d=2\) are allowed. Every positive root sum can replace \(d\) by common scaling. Rational target ratios admit rational coefficients and rational sign brackets; the roots themselves are not asserted to be rational.

Let \(p_1=1\), \(p_k=\prod_{j=2}^k c_j\), and
\[
\kappa_k=\frac{\binom d{k-1}\binom d{k+1}}{\binom dk^2}.
\]
Choose
\[
\lambda\ge9\max_{1\le k\le d-1}\frac{\kappa_k}{p_k},
\qquad R_k=\lambda p_k.
\]
The recurrence \(h_0=h_1=1\), \(h_{k+1}=h_k^2/(R_kh_{k-1})\) has the closed solution
\[
h_k=\lambda^{-\binom k2}
\prod_{j=2}^{k-1}c_j^{-\binom{k-j+1}2}.
\tag{1}
\]
This follows either by taking second differences of logarithms or by substituting the exponents directly. Consequently \(a_k=\binom dk h_k\) is positive, \(a_0=1,a_1=d\), and
\[
\frac{a_k^2}{a_{k-1}a_{k+1}}
=\frac{R_k}{\kappa_k}\ge9.
\tag{2}
\]
The division by \(\kappa_k\), rather than multiplication, is essential to the normalization.

Set \(q_k=a_{k-1}/a_k\); (2) gives \(q_{k+1}\ge9q_k\). For \(P(z)=\sum_{j=0}^d a_jz^j\), inspect \(P(-3q_k)\). Relative to the magnitude of the \(k\)th term, the nearest term on either side has ratio at most \(1/3\). Each subsequent step away from that term has ratio at most \(1/3\) as well. The finite left tail is therefore strictly below one half of the \(k\)th term, and the finite right tail is strictly below one half; a missing tail contributes zero. Their sum is strictly smaller than the central term even if equality holds in (2).

Thus \(P(-3q_k)\) has sign \((-1)^k\). At zero its sign is positive, and the sample magnitudes are strictly increasing. There is a real root in each of the \(d\) disjoint intervals between consecutive samples. Since the degree is \(d\), these roots exhaust it and are all simple and negative. This argument proves the required real-root statement directly; it invokes no external coefficient criterion and makes no sharpness claim for the conservative constant nine.

Reversal \(N(n)=n^dP(1/n)\) has the same distinct negative-root property and is monic. Its elementary symmetric root parameters are exactly the \(a_k\), so its normalized ratios are the prescribed \(R_k\), and all circuit ratios are exactly \(c_k\). The coefficient \(a_1=d\) fixes the root sum. For rational targets the choice of \(\lambda\) above is rational, as are all coefficients and samples. The reciprocal of the zero sample is an infinite endpoint for the reversed polynomial; if finite rational endpoints for \(N\) are desired, use \(-d\) for that end, since every \(r_i\) is strictly smaller than its positive sum \(d\).

For a prescribed word with letters \(-,0,+\), take respectively \(c_k=1/2,1,2\). The choice \(\lambda=9\cdot2^{d-2}\) works for every word because \(p_k\ge2^{-(d-2)}\) and \(0<\kappa_k<1\). Every ternary sign word, including exact ties, is therefore realizable. Positive circuit ratios never vanish: a zero sign here means \(c_k=1\).

## Predicate and lost information

The target is the complete vector of normalized circuit ratios, not just its signs. Distinct negative roots and a fixed first elementary moment are preserved in the construction. The free common magnitude \(\lambda\) of the Newton ratios is the missing coordinate that makes realization possible.

The second moment is not preserved. For root parameters \(1/2,1,3/2\), the first normalized ratio is \(12/11\). If the same first two elementary moments are imposed, a target \(c_2=1/2\) would force \(R_2=6/11<1\), contradicting the Newton inequality for positive real root parameters. Likewise the actual degree-five anchored model with sum thirteen and sum of squares fifty-nine has \(e_2=55\) and \(R_1=338/275\); the same target would force \(R_2<1\). These are exact nonempty moment-fibre hostiles. The theorem therefore does not close the anchored Laurent model, preserve an original factorial row, or prove anything about a fixed wall-stripped polynomial.

This is compatible with the separated-cluster theorem: that theorem imposes quantitative width and gap restrictions on its actual roots, while this realization changes the whole coefficient row. It supplies a stronger answer for unrestricted positive-root circuit words, with a correspondingly explicit loss of the moment sidecar.

## Independent exact verification

The [referee](../../04-computation/continuing8_20260906_newton_universality_audit.py) uses (1), not the producer's recursive \(h\) construction. It clears denominators and contents to primitive integer coefficient rows before independently reconstructing normalized Newton ratios. A separate product formula gives the sample addresses; direct sums verify every strict central-term dominance.

For every ternary word at each degree two through eight, it reconstructs the entire producer record and its complete-bank SHA-256 digest: **1,093 rows**, no orbit quotient or sign-only comparison. It also reconstructs all saved arbitrary positive targets, uses three fresh nonternary vectors, and independently scales the root sum to \(13/2\). There are **45 exact Sturm root counts and squarefreeness checks**, including the small complete universes and fresh targets. Both second-moment hostiles are retained.

All **25,328 always-active gates** pass in normal and optimized Python. Outputs were captured as raw subprocess bytes; they are identical and contain no carriage returns.

~~~text
python continuing8_20260906_newton_universality_audit.py
python -O continuing8_20260906_newton_universality_audit.py
~~~

The outside default locates the sibling producer folder. A --producer-dir argument is available; after filing in 04-computation, the default source and results paths are repository-relative.

Referee source SHA-256: 5ec2f72128e336b3af97010d52a19efbaac7272d7a15b98ce918e0bf52e9a379.

Referee output SHA-256: ce48e09f8adc41ea333a9717b92105e52534e2e07871bfd1bcbe48aed3281fe7.

Pinned producer source SHA-256: 1fc809f98a154e75dff50bf20c43e8a68418a438a0af9ec08d19377f10e13700.

Pinned producer certificate SHA-256: 8d11f09b961beb042f4146f55459141a2372275995b7b476cc85377633f6c4b0.
