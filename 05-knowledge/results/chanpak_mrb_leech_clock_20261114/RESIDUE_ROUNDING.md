# A late corridor has one inverse candidate, not a whole search

The earlier [exact frame decoder](https://github.com/eliottcassidy2000/math/blob/edf094e6bc983f5e8c8599aecdd2e0c00783773c/05-knowledge/results/chanpak_frame_decoder_pressure_20261114/RESETS_AND_DECODER.md)
fixes the corridor lengths and solves for outer widths. Here fix the
left/right reservoir profiles \(u,v\), nonnegative nonzero rows,
and instead solve for a missing corridor length \(q\):
\[
                        H_q=u C_q v^t .
\]
The two reservoir profiles can be arbitrarily large. No loss gate or
uncontrolled approximation is used in the final verification.

## Rank-one residues at all three poles

Let \(p(r)=r^3-r^2-r-1\), with largest real root \(\rho\) and complex
roots \(\alpha,\bar\alpha\). The complex pair has absolute value
\(\rho^{-1/2}\). The corridor's matrix generating function is
\[
 \sum_{q\ge0} C_q z^q =
 \frac{1}{1-z-z^2-z^3}
 \begin{pmatrix}1+z^2&1+z\\1+z&1+z+z^2\end{pmatrix}
 = \frac{N(z)}{1-z-z^2-z^3}.
\]
The numerator obeys
\(\det N(z)=-z(1-z-z^2-z^3)\). Thus **each** partial fraction
contributes a rank-one matrix, even at the nonreal poles. Explicitly
\[
 C_q=\rho^q P+\alpha^q R+\bar\alpha^q\bar R,\qquad
 P=\frac{\rho^2 N(1/\rho)}{p'(\rho)},\quad
 R=\frac{\alpha^2 N(1/\alpha)}{p'(\alpha)}.
\]
All entries of \(P\) are positive. The dominant contribution
\(uPv^t\) therefore factors into separate positive functions
of the two reservoirs. This is an exact rank-one **leading term**,
not an exact product formula for a finite ordinary count.
Only lengths 0 and 13 have that latter rank property.

For explicit uniform bounds, \(1.83<\rho<1.84\), \(|\alpha|<0.74\),
and \(\operatorname{Re}\alpha=(1-\rho)/2\).
From \(|\alpha|^2=1/\rho\) follows
\(\operatorname{Im}\alpha>0.6\) for the root in the upper half-plane.
Consequently
\[
 |p'(\alpha)|=|\alpha-\rho||\alpha-\bar\alpha|>2.4,\qquad
 |\alpha|^2+|\alpha|+1<2.4.
\]
Thus every entry of \(R\) has absolute value below 1. The
smallest numerator entry of \(P\) is \(\rho^2+1\), and
\(2(\rho^2+1)>p'(\rho)>0\), so every entry exceeds \(1/2\).
The two nonnegative profiles cannot amplify the *relative* error:
\[
 \left|\frac{H_q}{\rho^q uPv^t}-1\right|<4\rho^{-3q/2}.
\]

## Uniform inverse and one exact check

For all \(q\ge4\),
\[
 4\rho^{-3q/2}< 1-\rho^{-1/2} < \rho^{1/2}-1 .
\]
The first inequality follows, for example, from
\(4/(1.83)^6<1/8\) and \(1-\rho^{-1/2}>1/4\);
the last is \(x+x^{-1}>2\) at \(x=\sqrt\rho>1\).
This gives
\[
 q-\tfrac12 < \log_\rho\!\left(\frac{H_q}{uPv^t}\right)
                         < q+\tfrac12 .
\]
To decode a prescribed positive integer \(N\), test \(q=0,1,2,3\)
and the single remaining candidate
\(\operatorname{round}\log_\rho(N/(uPv^t))\) if it is at least 4.
There is no floating-point correctness assumption: compare
\(N^2\) to \((uPv^t)^2\rho^{2k\pm1}\) in the cubic field, and
then check the candidate with the exact integer matrix recurrence.
A half-integer tie can never correspond to a genuine late witness.

The conceptual link with the historical MRB inverse-residue work is
the distinction between *predicting the scale* and *deciding the
integer*. A leading analytic term here genuinely determines the
only possible late index, because the error is proved uniformly
smaller than half a step. It does **not** determine whether the
candidate is an exact realization: that last test is an integer
constraint. Also the method fixes both reservoirs; choosing them
for every target remains open. In a frame with both widths and a
corridor unknown, it combines with the quadratic wing inverse
after conditioning on complementary coordinates, but does not
remove their joint search cost.
