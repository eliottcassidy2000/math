# Exact spectral compilation of odd square-class lifts

Session120, 2026-10-30. Status: proved uniform identities plus exhaustive finite census; no Ramsey-number improvement is claimed.

## Scope and provenance

This continues Session119 `SQUARE_QUOTIENT_CLASSIFICATION.md` and its valuation-stratum classification, not merely the unit stratum. Let

\[
N=\prod_p p^{A_p}>1,\qquad p\equiv3\pmod4.
\]

For every divisor m>1 choose an arbitrary antipodal odd Boolean function
\(F_m:\{\pm1\}^{\omega(m)}\to\{\pm1\}\).
The symbol, including its actual arcs, is

\[
\sigma(0)=0,\qquad
\sigma((N/m)t)=F_m((\chi_p(t))_{p\mid m}),\quad t\in U(m).
\]

Arc a→b means σ(b−a)=+1. These are exactly the square-unit-invariant circulant tournament symbols in the stated prime class. Counts here are labelled symbols, not isomorphism classes.

The ingredients below are classical finite Fourier analysis, Boolean character expansion, CRT, Ramanujan sums, and quadratic Gauss sums. The contribution is their complete valuation-preserving compiler and its exact census application to the newly opened nonlinear family. For classical definitions and Fourier/Ramanujan identities see [NIST DLMF §27.10](https://dlmf.nist.gov/27.10); the proof below supplies all local evaluations needed, without relying on a choice of square-root sign.

## Theorem 1: uniform local-to-global spectral formula

Write
\[
F_m(s)=\sum_{E\subseteq\{p:p\mid m\}}\widehat F_m(E)\prod_{p\in E}s_p,
\quad
\widehat F_m(E)=2^{-\omega(m)}\sum_sF_m(s)\prod_{p\in E}s_p.
\]
Only odd-cardinality E have nonzero coefficient, by replacing s with −s.
For p^a exactly dividing m put v_{m,p}=(m/p^a)^{-1} modulo p^a. Define
\[
\tau_p=\sum_{u=1}^{p-1}\chi_p(u)e^{2\pi iu/p},\qquad \tau_p^2=-p.
\]
Then the eigenvalues of S[a,b]=σ(b−a), using Fourier vectors x_b=e^{2πikb/N}, are
\[
\boxed{\lambda_k=\sum_{m\mid N,m>1}\sum_{E\text{ odd}}
\widehat F_m(E)\prod_{p^a\parallel m}L_{p^a,\mathbf1_{p\in E}}(kv_{m,p})}.
\]
Here
\[
L_{p^a,0}(h)=\begin{cases}
p^{a-1}(p-1)&p^a\mid h,\\
-p^{a-1}&p^{a-1}\mid h,\ p^a\nmid h,\\
0& p^{a-1}\nmid h,
\end{cases}
\]
and
\[
L_{p^a,1}(h)=\begin{cases}
p^{a-1}\chi_p(h/p^{a-1})\tau_p&v_p(h)=a-1,\\
0&\text{otherwise}.
\end{cases}
\]
The v_p(0)=∞ convention applies.

**Proof.** Fourier diagonalization gives λ_k=Σ_dσ(d)e^{2πikd/N}; partition d by additive order m. Insert the Boolean expansion. CRT writes t=Σ_p t_p(m/p^a)v_{m,p}, separating the sum into the displayed local factors. For exponent0, subtract the sum over multiples of p from the unrestricted additive character sum. For exponent1, write each unit as u+pj, 1≤u≤p−1. Summing over j vanishes unless p^{a−1}|h; it then leaves p^{a−1} times the quadratic Gauss sum at h/p^{a−1}. At a nonzero residue the substitution u↦u/h gives χ_p(h)τ_p; at zero the character sum is0. Finally conjugation gives bar(τ_p)=χ_p(−1)τ_p=−τ_p. To compute its norm, substitute x=ty in Σ_{x,y≠0}χ_p(x)χ_p(y)ζ_p^{x−y}: the inner sum over y is p−1 at t=1 and −1 otherwise. Since Σ_tχ_p(t)=0, |τ_p|²=p. Thus τ_p²=−p. ∎

This avoids the more delicate classical evaluation of the *sign* of the Gauss sum. Exact arithmetic needs only the actual Gauss sums and their squares.

## Theorem 2: one frequency per divisor suffices for every trace moment

In K=Q(τ_p:p|N), the products τ_E=∏_{p∈E}τ_p form a rational basis; write [z]_∅ for its constant coefficient. Distinct prime valuations show the square classes −p are independent, hence each independent sign change of the generators is a field automorphism. For every integer j≥1,
\[
\boxed{\operatorname{tr}(S^{2j})=
\sum_{d\mid N,d>1}\varphi(d)\,[\lambda_{N/d}^{\,2j}]_\varnothing.}
\]
In fact the same formula is valid for positive odd powers (both sides vanish).

**Proof.** λ_0=0 by antipodality. At k=N/d, a nonzero local factor involving τ_p requires v_p(k)=a−1, with a≤A_p. Therefore p divides d: λ_{N/d} belongs to Q(τ_p:p|d). All frequencies of additive order d are (N/d)u, u∈U(d). The automorphism ζ_d↦ζ_d^u sends the first eigenvalue to this frequency's eigenvalue and acts on τ_p as multiplication by χ_p(u). CRT gives every sign pattern equally often. Averaging any power over these patterns kills every nonconstant basis product, proving the formula. ∎

Exact multiplication is particularly simple:
\[
\tau_E\tau_F=\Big(\prod_{p\in E\cap F}(-p)\Big)\tau_{E\triangle F}.
\]
Thus all moments can be computed using rational arithmetic, without sampling complex eigenvalues. This is divisor-orbit compression, not a polynomial-time claim in logN: the Boolean response tables themselves can be large.

## Corollary: exact fourth-moment objective at N=231

For λ=aτ3+bτ7+cτ11+dτ3τ7τ11 put
\[
C=-(3a^2+7b^2+11c^2+231d^2).
\]
Direct multiplication gives
\[
[\lambda^4]_\varnothing=C^2+
21(2ab-22cd)^2+33(2ac-14bd)^2+77(2bc-6ad)^2.
\]
All valuation strata m=3,7,11,21,33,77,231 are included. Their response counts are2,2,2,4,4,4,16, giving8192 symbols. The first six strata allow only signed odd characters. At m=231, eight responses are signed odd characters and eight are nonlinear majority-type responses. The exhaustive result is:

| Unit response | Symbols | Minimum tr(S⁴) | Minimizers | Maximum | Distinct moments |
|---|---:|---:|---:|---:|---:|
| Signed odd character |4096|17,225,670|16|645,693,510|468|
| Nonlinear |4096|29,082,438|8|175,030,086|429|

**Conclusion:** the first nonlinear enlargement does not improve the best fourth-moment value. This is a useful obstruction to a proposed spectral route, not an obstruction to every Ramsey objective. In particular neither tr(S⁴) nor its ordering determines the largest transitive subtournament, and no Ramsey number follows here.

The receipt includes complete moment histograms, not only optima. Its minimizer indices refer to the fully specified `responses` enumeration: representatives start with+1, remaining signs use lexicographic order −1,+1, and representative values use the same order. The two minimizing tuples, in ascending m order, are (0,0,0,1,3,1,6) and (0,0,0,3,2,2,11). The code reconstructs the actual arc masks from these tables.

For direct mathematical reconstruction, these witnesses have the following responses (s_p denotes the Legendre coordinate at p):

| m | Character minimizer F_m | Nonlinear minimizer F_m |
|---|---|---|
|3|−s₃|−s₃|
|7|−s₇|−s₇|
|11|−s₁₁|−s₁₁|
|21|s₇|s₃|
|33|s₃|−s₁₁|
|77|s₁₁|−s₁₁|
|231|−s₃s₇s₁₁|(s₃+s₇−s₁₁+s₃s₇s₁₁)/2|

## Independent verification

The Gauss/compiler algorithm computes every moment using the closed quartic above and rational local factors (scaled by2 at N231). Independently, **all8192 actual symbols** are built as positive-arc bit masks B. For h≠0, oddness implies
\[
C(h)=\sum_x\sigma(x)\sigma(x+h)=2-N+4|B\cap(B+h)|,
\quad C(0)=N-1.
\]
Since SSᵀ has these correlation entries and Sᵀ=−S,
\[
\operatorname{tr}(S^4)=N\sum_h C(h)^2.
\]
Bit rotation/population count evaluates this expression with no spectral formulas. Each of8192 values is checked against the compiler. Both minimizers are additionally checked by a naive scalar double loop. All52 symbols at N=3,7,9,21,27,33,49 are checked by the generic rational compiler against naive scalar autocorrelation. This includes repeated prime powers and mixed-prime strata; N21 recovers the four Session119 moments15876,19908,23940,44100, each four times.

Run:
```
python check_spectral.py > receipt.json
python -O check_spectral.py > receipt_optimized.json
cmp receipt.json receipt_optimized.json
```
All checks use explicit exceptions rather than removable `assert` statements. Finite verification supports the formulas but does not replace their all-parameter proofs above.

## Next boundary

A distinct, genuinely Ramsey-sensitive question remains: whether nonlinear responses improve the minimum possible maximum transitive-subtournament size. The current negative fourth-moment result neither answers nor licenses an answer to that question. The compiler now makes higher-moment constraints available uniformly without erasing the nonunit strata, and the complete fourth-moment census can prioritize—but cannot certify—such a search.
