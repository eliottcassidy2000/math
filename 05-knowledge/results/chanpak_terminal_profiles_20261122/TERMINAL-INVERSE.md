# Exact terminal fibers and a rank-one zero-run reset

This is an absolute inverse inside the transitive-gap [backchain](../chanpak_absolute_descent_20261122/BACKCHAIN.md). At a positive boundary write the physical state as \((D,0,G)\), with height \(F=D+G\). Then \(G\) is the immediately preceding prefix height, \(D\) is positive even, \(G\) is positive odd, and \(D\ge2G\). Put \(Z(D,U,V)=(D+2U,V,D+U)\).

## The last positive width needs no guessing

Suppose the target pair is \((N,G')\) at its last positive gap. The state immediately before it has height \(G'\) and middle coordinate \(U_*\). Therefore
\[
 L=N-G'=2^m(G'+U_*),\qquad 0\le U_*<G'.
\]
The strict upper bound follows from \(D_*+V_*>0\) on every state following an initial positive. So there is **at most one** possible ordinary integer width:
\[
 m=\left\lfloor\log_2(L/G')\right\rfloor.
\]
Reject if \(m<1\) or \(2^m\nmid L\). Otherwise \(U_*=L/2^m-G'\) is forced. This extends the earlier odd-part rule to arbitrary zero runs. The middle coordinate vanishes exactly for an immediate reset (or the first positive gap), since a nonempty zero run starting at \(D\ge2G>0\) has positive middle coordinate. In the reset case one reverts to the strictly smaller scalar \(G'\); its earlier second count is immaterial.

For a nonempty zero run of length \(j\), write
\[
 (a_j,c_j)=(\boldsymbol 1^\top Z^j e_1,(Z^j e_1)_2),
 \quad
 (b_j,d_j)=(\boldsymbol 1^\top Z^j e_3,(Z^j e_3)_2).
\]
Then the prior pair has \(D,G\) solving
\[
 a_jD+b_jG=G',\qquad c_jD+d_jG=U_*.
\]
If the determinant is nonzero this gives a unique rational solution; test it for integral parity and the positive-boundary cone \(D\ge2G\), then recurse. If it is zero, test compatibility and list the bounded integer solutions on the remaining line. Since \(Z^3\ge2I\), only \(O(\log G')\) values of \(j\) are possible, and the recovered older height never exceeds \(G'<N\). The direct reset can be handled by a scalar recursive call. A height-only target can also end in \(j\) zeros: for this case start with the single exact equation \(a_jD+b_jG=N\) in the same cone, and try the legal pairs it gives. Again \(j=O(\log N)\); the leading-zero normalization in [RESET-PRESSURE.md](RESET-PRESSURE.md) makes these cases exhaustive. This enumerates literal witnesses or certifies no witness in the grammar; it does not guarantee an accepting pair.

## Thirteen zeros really do lose the second count

The determinants start
\[
 0,2,-3,1,4,-8,5,7,-20,18,9,-47,56,0,\ldots
\]
and follow \(\Delta_{j+3}=-\Delta_{j+2}-\Delta_{j+1}+\Delta_j\).
Indeed the exterior-square transfer for \(Z\) has that characteristic polynomial; the determinant is a fixed linear read on its iterates. In particular at \(j=13\) the equations above have the proportional rows
\[
 (a_{13},b_{13})=47(56,47),\qquad
 (c_{13},d_{13})=9(56,47).
\]
Set \(K=56D+47G\). A thirteen-zero terminal run sends the height to \(47K\); with a following positive gap of width \(m\) it sends the pair to
\[
 (N,G')=((47+56\cdot2^m)K,\ 47K).
\]
This is a second genuine rank-one reset, distinct from immediate positive-to-positive multiplication. Any later linear suffix preserves the factor \(K>1\), so no prime height can pass through this run followed by a positive; nor can a prime stop at this run. A longer run may recover linear information from the unused coordinate, so the statement does not exclude longer uninterrupted runs on this ground. The finite comparison and this reset certify why a terminal inverse cannot generally recover the preceding second count just from a closed-form rational formula for each run.

## With isolated zeros, the gcd is the earlier product

In the normalized isolated-zero grammar the coupled step with \(w=2^m\) has
\[
\gcd(2(1+w)F-G,\,2F-G)=\gcd(2wF,2F-G)=\gcd(F,G),
\]
since the second argument is odd. A direct positive step instead sets the new gcd to the earlier height \(F\). Therefore the terminal gcd is precisely the product of all completed earlier atoms in this restricted grammar. Divide the pair by this gcd. The remaining primitive suffix has a unique exact backward history if it has any: when its second read is 1 it must be the initial positive step and the first read minus one must be a power of two of size at least 2. Otherwise the last step is coupled and its preceding state is forced by
\[
 m=v_2(F-G)-1,\qquad F_{\rm old}=(F-G)/2^{m+1},\qquad
 G_{\rm old}=2F_{\rm old}-G.
\]
Reject when parity, positivity, or the cone \(3G_{\rm old}\le F_{\rm old}\) fails. The terminal pair is realized exactly when this normalized descent succeeds and its gcd is a product of scalar isolated-zero atoms. Earlier scalar factorization need not be unique. This identifies what information the compatible second count carries: it separates earlier multiplicative pieces from the last coupled suffix.

## Independent replay

The backward and forward solvers agree on every positive boundary fiber through 1001. The independent census through 100001 gives 13,263 heights and 52,351 positive-boundary pairs; 1,172 of the 9,592 primes in range occur. These finite results are checks of the inverse; the sparsity theorem uses the infinite pressure estimate instead.
