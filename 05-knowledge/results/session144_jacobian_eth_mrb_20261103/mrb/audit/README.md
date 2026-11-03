# Independent MRB precision and Abel finite-part audit

This audit is independent of the inherited programs' numerical output. It does not edit or execute the source repository. The only inherited mathematical input to the precision certificate is the explicitly stated Euler/Cauchy remainder estimate below. No ETH claim is made.

## 1. Periodic cutoff and Abel finite parts: the claimed sign is correct

Put `f(n) = exp(log(n)/n) - 1`. Let `a_n` be a fixed real L-periodic sequence with mean `mu`, and define

\[
F(x)=\sum_{1\le n\le x}a_n f(n).
\]

Then there is a real constant `K_a` such that, for real `x -> infinity`,

\[
F(x)=\frac\mu2\log^2x+K_a+O_a(\log^2x/x).
\]

Proof. Write `a_n=mu+b_n`, so the partial sums of `b_n` are bounded. With `h(n)=log(n)/n`, Taylor's theorem gives `f(n)-h(n)=O(log^2(n)/n^2)`, and its tail is `O(log^2(x)/x)`. The integral-test sum/integral difference for `h` converges, with tail `O(log(x)/x)`, since `h` is eventually decreasing to zero and `integral h = log^2(x)/2`. Summation by parts gives convergence of `sum b_n h(n)`, with tail `O_a(log(x)/x)`. Combining these facts proves the stated formula at integer cutoffs. Replacing `floor(x)` by `x` costs `O_a(log(x)/x)`. This proves existence of exactly the cutoff constant asked for, without invoking an analytic continuation convention.

For `t>0`, the exponentially damped sum is absolutely convergent. Its Stieltjes integration by parts is

\[
A_a(t)=\sum_{n\ge1}a_nf(n)e^{-tn}
=t\int_1^\infty e^{-tx}F(x)\,dx.
\]

There is no lower-end boundary term: `f(1)=0`, and `F(x)=0` on `0 <= x < 1`. Put `ell=log(1/t)` and substitute `u=tx`. Replacing `F(x)` by its displayed asymptotic yields the main integral

\[
\int_t^\infty e^{-u}\{\tfrac\mu2(\ell+\log u)^2+K_a\}\,du.
\]

Extending this integral down to zero changes it by `O_a(t)`: substitute `u=tv` on the omitted interval, where `ell+log(u)=log(v)`, and use integrability of `log^2(v)` on `(0,1)`. The remainder contributes

\[
O_a\!\left(t\int_1^\infty e^{-tx}\frac{\log^2x}{x}\,dx\right)
=O_a(t(1+\ell^3))=o(1).
\]

Alternatively, only `F(x)-mu log^2(x)/2 -> K_a` is needed: split its bounded remainder at a fixed large cutoff and use the exponential kernel's unit mass.

Finally,

\[
\int_0^\infty e^{-u}\,du=1,\qquad
\int_0^\infty e^{-u}\log u\,du=-\gamma,\qquad
\int_0^\infty e^{-u}\log^2u\,du=\gamma^2+\pi^2/6.
\]

These follow by differentiating the defining Gamma integral at 1. Therefore

\[
\boxed{\lim_{t\downarrow0}\left[A_a(t)-\mu\{\ell^2/2-\gamma\ell\}\right]
=K_a+\mu\{\gamma^2/2+\pi^2/12\}.}
\]

In particular the constant correction is positive for `mu>0`; neither the sign of the `gamma*ell` term nor the sign of the finite correction should be reversed. When `mu=0`, this is ordinary Abel convergence to the convergent periodic weighted sum. This is a fixed-filter analytic theorem, not a uniform complexity statement for arbitrary unspecified real filters.

## 2. Exact Euler weights have only linear absolute-error amplification

For the fixed constant

\[
C=\sum_{n\ge1}(-1)^n f(n),
\]

the inherited Euler approximation is the prefix through `M-1` plus `(-1)^M` times

\[
S_{M,N}=\sum_{j=0}^{N-1}\frac{(-1)^j\Delta^j f(M)}{2^{j+1}}
=\sum_{k=0}^{N-1}(-1)^k w_kf(M+k).
\]

Here the weights are **exact dyadic rationals**:

\[
w_k=\frac{W_k}{2^N},\qquad
W_k=\sum_{j=k}^{N-1}\binom jk2^{N-j-1}
=\sum_{j=k+1}^{N}\binom Nj.
\]

The last identity follows from Pascal's identity (or by induction on `N`). Consequently

\[
0<W_k<2^N,\qquad \sum_{k=0}^{N-1}W_k=N2^{N-1},\qquad
\sum_{k=0}^{N-1}|w_k|=N/2.
\]

Thus replacing each evaluated `f(n)` by a value with absolute error at most `epsilon` perturbs the entire prefix-plus-tail approximation by at most

\[
(M-1+N/2)\epsilon.
\]

There is no exponential error amplification in this flattened representation. This is an absolute-error statement, not a claim of relative accuracy in individual small differences or a claim that arbitrary floating-point difference-table implementations attain it. Compute `W_k` by exact binomial recurrence and accumulate exact dyadic products to realize the bound.

## 3. A universal precision rule

The inherited analytic tail bound, valid for `M>=4, N>=1`, is

\[
B_{M,N}=2N!e^b\left[\frac{L}{M^{N+1}}+
\frac{L/N+1/N^2}{M^N}\right],\quad
L=\log(3M/2)+\pi/6,\quad b=2L/M.
\]

Let the requested precision be an integer `p>=1`, and choose

\[
N=p+4,\qquad M=4N,\qquad
q=p+\lceil\log_2(M+N)\rceil+3.
\]

Since `N!<=N^N`, `L/N<3`, and `e^b<6`,

\[
B_{M,N}
\le2e^b4^{-N}\left[\frac54\frac LN+\frac1{N^2}\right]
<57\,4^{-N}<64\,4^{-N}
=2^{-2p-2}\le2^{-p-3}.
\]

For completeness, `log(6N) <= log 6 + N-1` and `log 6+pi/6<3` give `L/N<3`; then `b=L/(2N)<3/2` and `e<3` imply `e^b<sqrt(27)<6`.

If each term is replaced by a dyadic approximation within `2^-q`, the total arithmetic input error is less than `2^(-p-3)`, because `M-1+N/2 < M+N`. Exact integer accumulation introduces no additional error. Round the result to the nearest multiple of `2^(-p-2)`; this contributes at most `2^(-p-3)`. The final output obeys

\[
|\widetilde C-C|<3\,2^{-p-3}<2^{-p}.
\]

This rule is intentionally simple, not parameter-optimal.

## 4. Explicit certified primitive implementation and bit model

`exact_precision_audit.py` implements the above evaluator using Python integers and exact rational arithmetic only. It does not treat a call to floating-point `log`, `exp`, or `mpmath.iv` as a free certified primitive. Floating point occurs only when formatting a diagnostic ratio and elapsed time, not in any certificate or acceptance decision.

For each integer `n`, range-reduce `n=2^k*y`, `1<=y<2`, and use

\[
\log y=2\sum_{j=0}^{J-1}\frac{z^{2j+1}}{2j+1}+R_J,
\quad z=(y-1)/(y+1)\in[0,1/3],\quad 0\le R_J<9^{-J}.
\]

Use the same series for `log 2` with `z=1/3`. With `s=ceil(log2(k_max+1))` and `J=q+s+4`, the resulting log enclosure has width less than `2^(-q-12)`. Divide by `n`. Both endpoints lie in `[0,1]`: the log enclosure excess is less than 1, whereas `n-log n>=1`.

For `0<=x<=1`, use `T=q+6` and

\[
\sum_{j=0}^{T-1}\frac{x^j}{j!}\le e^x
\le\sum_{j=0}^{T-1}\frac{x^j}{j!}+\frac2{T!}.
\]

Evaluating the lower polynomial at the lower endpoint and the upper polynomial at the upper endpoint gives an `exp(x)-1` enclosure of width less than

\[
3\cdot2^{-q-12}+2^{-q-4}<2^{-q-3}.
\]

Here `T!>=2^(T-1)` and `e^x<=e<3` suffice. Round the rational interval midpoint to a `q`-fractional-bit dyadic. The endpoint distances are strictly less than `2^-q`; the program checks that exact condition for every term.

**Bit model.** The fixed object to be evaluated is `C`; the input is unary precision `p` (equivalently, complexity is measured as a function of the requested output bit count, not the binary encoding length of `p`). `M,N,q,J,T=O(p)` and `n=O(p)`. The atanh sums have rational numerator/denominator bit lengths `O(p log p)`: a common denominator can be bounded by a power of the rational argument denominator times `(2J-1)!`. The exponential sums then have numerator/denominator bit lengths `O(p^2 log p)`, using a power of that denominator times `(T-1)!`. All intermediate reduced fractions obey the same order bound. There are `O(p^2)` exact rational operations across `O(p)` sample values, plus at most `O(p^2)` integer operations in the independent weight checks. Ordinary quadratic-time integer multiplication/division and rational gcd operations therefore give a conservative bit bound `O(p^6 log^2 p)` (replace `log p` by `log(p+2)` at small inputs). This is more than sufficient to prove polynomial-time computability of this fixed real constant. It is not an optimized complexity claim. A faster certified fixed-point library could improve it, but is not needed for the theorem.

This statement does **not** claim a uniform polynomial algorithm for arbitrary real periodic coefficients, unrestricted analytic functions, symbolic identity recognition, or integer-index recovery at any arbitrary residual scale. Those would need separate input representations and conditioning/error analysis.

## 5. Stronger nonuniqueness obstruction: audit passed

Let `K>=1`, `P_K(z)=product_{j=1}^K(z-j)`, and choose any nonzero nonnegative `phi` in `C_c^infinity((1,2))`. Define

\[
h_K(z)=P_K(z)\int_1^2\phi(t)e^{-tz}\,dt.
\]

This is entire and real on the real line. It vanishes at `1,...,K`. If `a=inf(supp phi)>1`, then on the positive ray it is `O_K((1+x)^K exp(-a*x))`. For every integer `A>=0`, integrating by parts `m>=K+A` times, with no endpoint terms, gives

\[
\left|\int\phi(t)e^{-tz}\,dt\right|
\le\|\phi^{(m)}\|_1|z|^{-m}\qquad(\Re z\ge0,\ z\ne0).
\]

Since `P_K(z)=O_K((1+|z|)^K)`, enlarging the constant to cover `|z|<=1` proves

\[
|h_K(z)|\le C_{K,A}(1+|z|)^{-A}\qquad(\Re z\ge0).
\]

There is no contradiction with entire-function uniqueness: no comparable bound is imposed on the left half-plane. Absolute convergence justifies Fubini in the alternating sum. The generating-function identity is

\[
\sum_{n\ge1}(-1)^nP_K(n)q^n
=(-1)^{K+1}K!\frac{q^{K+1}}{(1+q)^{K+1}},\qquad 0<q<1.
\]

Consequently

\[
\sum_{n\ge1}(-1)^nh_K(n)
=(-1)^{K+1}K!\int_1^2\frac{\phi(t)}{(e^t+1)^{K+1}}\,dt\ne0.
\]

Thus `f+lambda*h_K`, for any nonzero real `lambda`, has the same branch jump, the same first `K` sampled values, and the same complete algebraic-order uniform asymptotics in the right half-plane, but a different alternating sum. Entire addition preserves the relevant slit-plane holomorphy and reality. These statements preserve the analytic/prefix/asymptotic data under discussion, not an additional governing functional equation or extra positivity restrictions unless those are separately checked.

The elementary family also checks out: for any integer `r>=1`,

\[
H_{K,r}(z)=e^{-z}P_K(z)\left(\frac{1-e^{-z}}z\right)^{K+r}
\]

is entire (the quotient takes value 1 at zero), vanishes at `1,...,K`, is exponentially flat on the positive ray, and is uniformly `O_K,r(|z|^-r)` on the right half-plane as `|z|->infinity`. Its alternating sum is

\[
(-1)^{K+1}K!\int_{[0,1]^{K+r}}
\frac{dt_1\cdots dt_{K+r}}{(e^{1+t_1+\cdots+t_{K+r}}+1)^{K+1}},
\]

which is nonzero. This follows from `(1-e^-z)/z=integral_0^1 exp(-tz) dt`. The numerator exponent `K+r` is the cube dimension; the denominator exponent stays `K+1`.

## 6. Execution and scope

Run the exact evaluator normally and under optimization:

```
python exact_precision_audit.py --precisions 1 8 16 32 --output exact_precision_results.json
python -O exact_precision_audit.py --precisions 1 8 16 32 --output exact_precision_results_optimized.json
```

The script uses explicit exceptions, no assertion statements. Each recorded row is a fully executed exact-rational enclosure certificate, conditional only on the inherited analytic Euler remainder theorem. The all-precision complexity statement is a proof about the algorithm, not a claim to have run every precision. The regularization and entire-bump results above are analytic proofs, not numerical interval experiments.

Final executed scope: `p=1,8,16,32` in normal and optimized Python. Eight hostile inputs were rejected in both modes: nonpositive precision; negative and oversized atanh arguments; negative and oversized exp arguments; corrupted weight; corrupted rounding; and an interval too wide for the claimed term precision. The optional `p=64` runs were stopped before completion and are not claimed as executed certificates. `compare_receipts.py` removes only elapsed-time and interpreter-optimization fields, requires the two mathematical receipts to match exactly, and emits a stable SHA-256 digest of the normalized content.
