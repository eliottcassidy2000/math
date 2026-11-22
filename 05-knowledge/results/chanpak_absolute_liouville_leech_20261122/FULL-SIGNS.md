# Squareclass separation, prime-power removal, and full signs averaged over quadratic twists

**Status: PROVED** for the squareclass, prime-power, random-sign, and quadratic-twist statements below. This does **not** give cancellation of the *untwisted* Liouville function (nor primes) on the literal pairs. Rather, it removes an absolute thin-support error in the signed reduction, and isolates a deterministic large-prime sign problem that survives even when squareclass degeneracies are negligible.

## 1. Work on the literal fixed-gap family

Use the [literal fixed-gap proof](../chanpak_prime_parity_inverse_20261122/FIXED-GAP-SIEVE.md). The independent binary blocks AB and BA have common scale \(a=490032\) and digit translations differing by \(2424\), and the fixed continuation gives the actual Hamilton counts \(N_w=F(v_w)\) and \(N_w+2\). At depth \(n\), the \(2^n\) words \(w\) give distinct \(N_w\), with \(N_w\asymp X=a^n\). The earlier sieve supplies a subfamily \(\mathcal W_n\) of size \(\asymp 2^n/n\) such that every prime factor of each endpoint exceeds \(z_n=\exp(c\sqrt n)\) for some \(c>0\). No step below changes the output map or replaces \(N_w+2\) by a decoded number. All estimates also apply to the unsifted word family when stated with its size.

If the first \(n-\ell\) binary choices are fixed, the \(2^\ell\) numbers in that cylinder lie in an interval of length \(O(a^\ell)\), by the common-scale digit expansion and the fixed final affine map.

## 2. An absolute thin-support bound for prime powers

**Claim.** Among words of depth \(n\), only \(O(2^{n/2})\) can give a proper integer power in either endpoint.

Indeed, for a fixed exponent \(r\ge2\), distinct positive \(r\)th powers of height \(\asymp X\) are spaced \(\gg X^{1-1/r}\). The constant is uniform in \(r\): increasing the positive base by one adds at least \(r\) times the smaller base to the power \(r-1\). In a cylinder with \(\ell=\lfloor(1-1/r)n\rfloor\), at most \(O(1)\) words can attain such an \(r\)th power. There are \(O(2^{n/r})\) cylinders. Exponents are \(O(n)\) unless the output is 1, which never occurs here. Summing gives \(O(2^{n/2}+n2^{n/3})=O(2^{n/2})\).

This is a deterministic gain for the untwisted signed equation. In any normalized rough-word average, replacing \(\Lambda(N_w+j)\) by \(\theta(N_w+j)=\log(N_w+j)1_{\mathrm{prime}}(N_w+j)\) costs \(O(n^2 2^{-n/2})\), even in the presence of a bounded Liouville weight at the other endpoint. An additional factor of \(\Lambda\) costs only another \(O(n)\). The old ambient prime-power bound can swamp this exceptionally thin word family; the cylinder-scale argument does not. In addition \(N_w(N_w+2)\) can never be a proper power: its odd factors \(N_w,N_w+2\) are coprime, and two distinct positive powers of the same exponent \(r\ge2\) cannot differ by two at this height.

## 3. Three parity probes and their squareclasses

For a positive integer \(B\), write \(\operatorname{sf}(B)=\prod_{v_p(B)\ {\rm odd}}p\). Use any of the three probes \(B_w=N_w,\ N_w+2,\ N_w(N_w+2)\). For either single-endpoint probe, any fixed class \(d\) gives values \(d u^2\); consecutive values at height \(\asymp X\) are spaced \(\gg\sqrt X\), uniformly in \(d\). Splitting into cylinders with \(\ell=\lfloor n/2\rfloor\) gives at most \(O(2^{n/2})\) words per class.

**The paired probe has no repeated squareclass at a fixed depth.** Set \(x=N_w+1\). For a fixed squareclass \(d\), the equation is
\[
 x^2-dy^2=1,\qquad x\asymp X,\quad x,y>0.
\]
There is no \(d=1\) solution at this height. If \(d>1\) and two solutions are ordered by \(x\), their units \(x+y\sqrt d\) have ratio a real unit \(r+s\sqrt d>1\), with integral coefficients and norm 1. Its trace satisfies \(r=(\eta+\eta^{-1})/2\ge2\), so the ratio is at least \(2+\sqrt3\), uniformly in \(d\). In the literal alphabet at one depth, the ratio of maximum to minimum increments is at most \(362956/360532<1.01\), since each base digit lies between those two translations. The fixed affine read and the addition of 1 are positive, so \(x_{\max}/x_{\min}<1.01\) as well. The units lie between \(x\) and \(2x\), making their ratio less than \(2.02<2+\sqrt3\). Therefore at most one word lies in each product squareclass, without any assumption on a fundamental Pell unit or on factorization of the outputs.

For any subfamily of size \(T_n\), the number of ordered pairs of words with the same class is therefore \(O(T_n2^{n/2})\) for a marginal and exactly \(T_n\) for the paired probe. In particular on the rough family the squareclass parity vectors of the primes above \(z_n\) span an \(\mathbb F_2\)-space of rank at least \(n/2-O(\log n)\) for either endpoint and \(n-O(\log n)\) for the product: a vector space of rank \(r\) has at most \(2^r\) distinct vectors. This is an information statement about the multiplicative group modulo squares, not cancellation of its particular all-primes sign character.

## 4. Full randomized prime signs

Let \(\xi_p\) be independent fair signs, one fixed choice per prime, extended completely multiplicatively to the positive integers. For arbitrary fixed complex weights \(u_w\) of \(|u_w|\le1\), orthogonality gives
\[
 \mathbb E_\xi\left|\frac1{T_n}\sum_w u_w\,\xi(B_w)\right|^2
 \le \frac{\#\{(w,w'):\operatorname{sf}(B_w)=\operatorname{sf}(B_{w'})\}}{T_n^2}.
\]
The two marginal variances are \(O(n2^{-n/2})\) and the joint variance is \(O(n2^{-n})\) on rough survivors. The same estimates hold after fixing the values \(\xi_p\) at primes up to \(z_n\), because the surviving integers contain none of them. One may, for example, take \(u_w=\lambda(B_w)\), so the actual Liouville value is included with all its prime factors before averaging over the twist. Taking fixed \(\xi_p\) at all primes simultaneously, Chebyshev and Borel–Cantelli show that for almost every such twist the three full correlations on these deterministic word families tend exponentially to zero. Fourier inversion on the two signs gives all four joint sign patterns with limiting fraction \(1/4\). No small-factor truncation appears in this averaged assertion.

## 5. Full Liouville cancellation for most genuine quadratic twists

Let \(b\) range over positive fundamental discriminants \(b\le D\) with \(b\) squarefree and \(b\equiv1\pmod4\). The associated quadratic character at positive odd \(t\) is the Kronecker symbol \(\chi_b(t)=(b/t)\). For each fixed \(\epsilon>0\), uniformly for \(D\ge X^{s/2+\epsilon}\), the mean square over those discriminants of \(\frac1{T_n}\sum_w u_w\lambda(B_w)\chi_b(B_w)\) is
\[
 \ll_\epsilon \frac{L_s}{T_n}+O_\epsilon(X^{-\delta}),
 \qquad
 (s,L_s)=(1,2^{n/2})\ \text{for an endpoint},\quad
 (s,L_s)=(2,1)\ \text{for the product},
\]
for some \(\delta>0\). The word family can be any deterministic subfamily; the rough family is where the negligible-zero statement below is uniform. Thus for \(D\ge X^{1+\epsilon}\) the three full Liouville correlation bounds hold simultaneously for a fraction \(1-o(1)\) of these genuine twists. This is not a theorem about a single fixed discriminant as \(n\) varies.

For the proof, expand the square. If \(B_wB_{w'}\) is a square use the bound 1 and Section 3. Otherwise \(q=\operatorname{sf}(B_wB_{w'})>1\) is odd and squarefree, with \(q\ll X^{2s}\), and \((b/q)\) is a primitive nonprincipal character of conductor \(q\) as a function of \(b\). Primes in the even-exponent part not dividing \(q\) impose only a coprimality condition; let their squarefree product be \(R\). Expand \(1_{(b,R)=1}\) by divisor inclusion–exclusion and impose squarefreeness of \(b\) with \(\mu^2(b)=\sum_{\ell^2\mid b}\mu(\ell)\). For an inclusion–exclusion divisor \(r\mid R\), the combined divisor is \({\rm lcm}(\ell^2,r)\), not necessarily their product. The remaining sum is over one residue class modulo 4, which is a translate of a consecutive character sum after extraction of that unit step, since \(q\) is odd.

Use [Burgess's cubefree-modulus bound](https://doi.org/10.1112/plms/s3-13.1.524),
\[
 \left|\sum_{M<h\le M+H}\chi(h)\right|
 \ll_{\nu,\eta} H^{1-1/\nu}q^{(\nu+1)/(4\nu^2)+\eta},
 \qquad \nu\ge3.
\]
Length \(H>q\) reduces to the remainder interval by vanishing complete sums. Nonunit extracted divisors give a zero summand. Since \({\rm lcm}(\ell^2,r)\ge\ell^2\), the sum of the length penalties from \(\ell\) converges (as \(2(1-1/\nu)>1\)), and the \(r\) sum costs at most \(\tau(R)\ll_\eta X^{O(\eta)}\). The normalizing number of fundamental discriminants is \(\asymp D\), by the same elementary squarefree expansion in \(1\bmod4\). This bounds the normalized nonsquare correlation by
\[
 O_{\nu,\eta}\left(D^{-1/\nu} X^{\,s(\nu+1)/(2\nu^2)+O(\eta)}\right).
\]
For \(D\ge X^{s/2+\epsilon}\), take \(\nu>s/\epsilon\) and then \(\eta\) small enough. This is a uniform power saving, so summing over all off-diagonal word pairs incurs no further loss after dividing by \(T_n^2\).

On rough words each \(\chi_b(N_w+j)\) is zero only when some prime factor of the endpoint divides \(b\). The normalized number of admissible \(b\) divisible by \(p\) is at most \(O(1/p)\) by bounding them by all integers divisible by \(p\). Since \(\Omega(N_w(N_w+2))=O(\sqrt n)\) on rough words, the average over twists and words of the zero rate is \(O(\sqrt n/z_n)\). Therefore for most twists the signs are nonzero on \(1-o(1)\) of the rough pairs, and the three correlation estimates imply that their four full joint sign frequencies tend to \(1/4\).

**Exact remaining obligation.** The estimates average over independent sign systems or over conductors growing with the output height. Neither implies that the particular untwisted character \(b=1\) is not exceptional. Even the absence of every squareclass repetition would not imply that: a large subset of distinct vectors can all lie in one affine hyperplane of the squareclass space. The unresolved step is the correlation of the deterministic all-primes parity functional with these actual tournament word vectors. The thin-support replacement of prime powers is unconditional and already applies to that fixed untwisted problem.
