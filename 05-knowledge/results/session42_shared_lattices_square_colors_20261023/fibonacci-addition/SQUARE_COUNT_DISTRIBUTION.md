# The full forced-A square-count distribution

This is a short consequence of RESULTS.md, not an additional equidistribution assumption. Let alpha=phi^(-2) and X_k=U_A(k^2)/k^2, where U_A counts unordered distinct canonical-A pairs.

The uniform logarithmic error and diagonal correction imply

X_k−ell({k^2 alpha})/2 ->0.

For every bounded continuous test function f,

lim_(K->infinity) (1/K) sum_(k<=K) f(X_k)
 =(1−2alpha)f(0)+4 integral_0^(alpha/2) f(x) dx.

Indeed ell is zero on an interval of length1−2alpha. On its support it rises with slope1 and falls with slope−1; the map ell/2 therefore has two inverse branches of derivative magnitude2. Square-phase equidistribution gives the asserted pushforward measure. It is a probability measure because (1−2alpha)+4(alpha/2)=1.

In particular, for every fixed0<t<alpha/2,

lim_(K->infinity) #{k<=K: X_k>t}/K=2alpha−4t.

The mean is alpha^2/2, the second moment alpha^3/6, and the variance alpha^3/6−alpha^4/4. The limiting distribution conditional on its positive part is uniform on(0,alpha/2). This distributional conditioning is not a uniform bound for every individual representable target; targets near the phase boundary can have much smaller counts.

The complete nine-entry vector R_cd(k^2)/k^2 likewise converges in distribution to the vector g_cd(Theta), with ONE uniform phase Theta shared by every entry. Its mean matrix is (lambda_c lambda_d), but the entries are not independent. This is a precise distinction between an averaged product formula and compatible per-target additive data.

For noncanonical branches, canonical A remains a subset, so the canonical count gives a lower envelope. The distribution above is asserted for the canonical count, not for every branch or for prime-restricted summands.
