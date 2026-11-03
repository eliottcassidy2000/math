# Independent audit: finite data and right-half-plane decay do not determine the MRB sum

Status: PROVED analytic statement; finite symbolic controls separately VERIFIED.

For an integer K>=0 define P_K(z)=product_{j=1}^K(z-j), E(z)=(1-exp(-z))/z with E(0)=1, and h_K(z)=exp(-z)P_K(z)E(z)^(K+1). E is entire by its power series or E(z)=integral_0^1 exp(-tz)dt. Thus h_K is entire, has zeros at 1,...,K, and has no branch jump. On the positive real axis it is exponentially small times O(1/x), so it does not alter any algebraic/logarithmic asymptotic expansion.

For Re z>=0 and |z|>=1, |1-exp(-z)|<=2 and |P_K(z)|<=(K+1)^K |z|^K (empty product 1). Consequently |h_K(z)|<=2^(K+1)(K+1)^K/|z|, uniformly in this half-plane, with the K=0 convention. This does NOT assert decay on the whole complex plane. It does not preserve the exact functional equation specifying z^(1/z), which is information beyond these observations.

Using E(n)^(K+1)=integral_[0,1]^(K+1) exp(-n sum(t_i))dt, absolute convergence permits summation inside the integral: polynomial-in-n times exp(-n) is summable uniformly on this cube. For a=1+sum(t_i),

sum_{n>=1}(-1)^n exp(-an) product_{j=1}^K(n-j)
= (-1)^(K+1) K!/(exp(a)+1)^(K+1).

Indeed the first K terms vanish; product=(K!)binom(n-1,K) for n>=K+1, and the geometric generating function gives the identity. Integrating proves the claimed nonzero alternating sum, of sign (-1)^(K+1). Scaling by any real number therefore changes the alternating functional by any prescribed real amount while retaining branch jumps, a chosen finite prefix, complete positive-axis power/log asymptotics, and the uniform right-half-plane O(log|z|/|z|) decay class of the MRB kernel.

The earlier Session31 source explicitly did NOT claim that its perturbation preserves this growth class. Its suggestion of a missing 'global growth/profile sidecar' was not a theorem that growth alone determines the function. The present result strengthens its obstruction; it does not refute an earlier uniqueness theorem. A full profile or a genuinely determining functional relation could still restore uniqueness.

Independent exact checks: check_mrb_boundary.py derives the rational generating functions for K=0,...,8 by applying (q d/dq-j) successively to -q/(1+q), rather than summing the producer formula. It checks prefix roots and sign and rejects a reversed-sign hostile. Normal and Python -O outputs agree. The all-K proof above is analytic, not inferred from these finite cases. No MRB irrationality/transcendence, modular identity, or causal relation to ETH or Jacobian counterexamples is asserted.

## Stronger smooth-density version, independently checked analytically

Take any nonnegative, nonzero smooth real φ compactly supported inside (1,2), and set h(z)=P_K(z)∫φ(t)exp(-tz)dt. Compact integration makes h entire. For each integer N>=0, integrate by parts N times; all boundary terms vanish. Uniformly for Re z>=0, the transform is bounded by |z|^-N∫|φ^(N)(t)|dt. Multiplication by P_K costs at most K powers; hence h is smaller than every inverse power uniformly in the closed right half-plane. There is no contradiction with Liouville because no left-half-plane bound is asserted. The same absolutely convergent geometric identity gives sum(-1)^n h(n)=(-1)^(K+1)K!∫φ(t)/(exp(t)+1)^(K+1)dt, nonzero with its stated sign. Thus even all-order uniform half-plane flatness leaves ambiguity. This is a mathematical existence proof using a standard smooth bump, not a finite-bit computation of its integral.
