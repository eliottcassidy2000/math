# Independent review: fixed-symmetry-slack renewal transfer

Date: 2026-10-21. Verdict: **PASS** for the transfer argument, conditional only on the stated fixed-slack compiler. The parent reports the compiler's independent review has also passed; this review does not duplicate that proof audit. No correctness revision to the transfer is required. No repository writes or Git operations were performed.

Audited source: `session23-root/SLACK_RENEWAL_TRANSFER.md`, SHA-256 `6e391cb643664a8055a9adc6a5741b45f9b36ba46d69d25961660bbcd29dfa30`.

Compiler at the time of these controls: `session23-connections/SYMMETRY_SLACK_RATIONALITY.md`, SHA-256 `22cb6aca899d9a485b2432d37a90d13040dfc6353ee41529a243d8489aa974a1`.

## 1. Moving-pole proof

Write p(z,u)=1-z-uz³ and let R(u) be its root continuing rho=R(1). The triangle inequality, together with positive terms at exponents 1 and 3, makes rho the unique smallest-modulus zero at u=1. The zero is simple, since p_z(rho,1)=−1−3rho². Consequently an isolated analytic R(u) and a uniform separation from all other roots exist on a complex neighborhood of u=1.

Let m=max ell over actual finite configurations. A nonempty slack layer has a witness and all its singleton prefixes, so it has infinitely many orders and m≥1. For L=A−D and

c(u)=1/[R(u)(1+3uR(u)²)],

the leading singular term of L is c(u)/(1-z/R(u)). The finite polynomial D changes lower Laurent terms but not this principal part. Thus the coefficient amplitude is

K(u)=c(u)^m/(m−1)! · sum_{ell_c=m} R(u)^{a_c} u^{b_c}.

It is analytic and strictly positive at u=1. It stays nonzero after shrinking the neighborhood. Finite partial fractions, including smaller pole orders, give uniformly

[z^n]F_d(z,u)=K(u)n^{m−1}R(u)^(−n)(1+O(1/n)).

Cauchy estimates justify differentiating the analytic relative error on a smaller neighborhood. The moment generating function consequently has logarithm

n log(rho/R(e^t)) + log(K(e^t)/K(1)) + O(1/n).

The source's implicit differentiations give mu=rho²/(1+3rho²) and sigma²=rho²(1+rho²)/(1+3rho²)³. They are correct. This establishes the claimed O(1) mean and variance errors and the characteristic-function CLT, not only weak convergence. The d=0 case follows directly with m=1 and D=0.

The exact equality s=3f−d immediately transfers the mean, variance, law of large numbers, and CLT to moved support. In particular Var(s)=9sigma² n+O_d(1).

## 2. Supported parity and variance doubling

The potential ±sqrt(rho) competition must not be treated as a single dominant z-pole. The source correctly restricts to each supported vertex parity.

Here is an explicit expansion that makes that step completely formal. Set x=z² and define D_q=0 for q<0. The palindrome tail is the sum, over (h,e)=(0,0),(1,0),(3,1), of

z^h y^e [A(x,y²)−D_floor((d−h)/2)(x,y²)].

Indeed h is the center's size, e its triangle contribution, and 2|v|+h>d is exactly the stated half-word cutoff. Each paired long slot contributes A(x,y²)−D_d(x,y²). After selecting a fixed decoration, optional center, and the matching vertex parity, every unbounded component family is therefore a finite positive monomial times a product of such word tails. Its only relevant moving x-pole is R(y²), with a positive leading coefficient. A configuration with no long half-slots contributes only finitely many sizes and is harmless.

For n=2N+r, the exponential rate in the normalized generating function is −N log R(e^{2t}). Dividing by n yields −(1/2)log R(e^{2t}) up to an O(1) analytic term. The first and second derivatives are mu and 2sigma². Thus the correct variance is doubled, not unchanged.

There may be different maximal pole orders on the two parities; neither the universal mean nor the variance constant changes. On every parity that has unbounded support, a nonempty positive tail family supplies all sufficiently large orders of that parity. Unsupported parities require no distributional assertion. Finite mixtures of families have O(1) differences between their means and add only O(1) between-family variance.

## 3. Actual automorphism order, not an orbit proxy

For a canonical decorated quotient, the internal kernel is the direct product of all triangle rotations, C3^t where t=sum_i t(w_i). Projection onto the actual word-preserving quotient stabilizer K0 is surjective. Choosing the same canonical model for each repeated word gives compatible lifts, so the extension splits:

Aut(T) ≅ C3^t ⋊ K0.

The semidirect action may permute triangle-rotation factors in short blocks. It does not change the group order, which is exactly 3^t |K0|. It is essential to use the actual stabilizer and canonical blocks, not to count visible triangles as independent factors.

Since f=t+C0, the exact logarithmic correction is

log|Aut(T)|−f log3 = log|K0|−C0 log3.

The compiler's bounds 0≤C0≤d and |K0|≤(4d+1)! make this uniformly bounded. The source uses that only for its law of large numbers and CLT. This is sound. The same bounded correction applies after restriction to self-converse objects.

The actual order-six control C3[11,11,11] has |Aut|=3, f=4, s=6, and delta=6. Its quotient action is diagonal C3, not two independent triangle rotations; 3^f=81 is also not its actual automorphism order.

## 4. The variance pitfall and an optional stronger result

A bounded random correction alone does not imply an O(1) correction to variance. For X_n=±sqrt(n) with equal probabilities and B_n=sign(X_n), one has |B_n|=1 but

Var(X_n+B_n)−Var(X_n)=2sqrt(n)+1.

With Var(X_n)=O(n), the general bounded-error estimate is only O(sqrt(n)). The source's explicit warning is therefore necessary and correct.

There is an optional strengthening available from the full compiler, rather than from boundedness alone. Put kappa_c=log|K0,c|−C0,c log3, a constant on each finite configuration. The exact log-automorphism moment series is obtained by weighting configuration c by exp(t kappa_c) and putting y=exp(t log3). This preserves a nonzero analytic leading amplitude and changes the moving root to R(exp(t log3)). Therefore the same proof gives

E log|Aut|=mu log3 · n+O_d(1),
Var(log|Aut|)=sigma²(log3)² n+O_d(1)

in the ordinary model, with twice that linear variance on supported self-converse parities. These stronger moment statements are not needed to pass the existing text, which does not claim them. They should not be attributed merely to its bounded-perturbation observation.

## 5. Converse-merged distribution

The self-converse denominator p(z²) has all its zeros of modulus at least sqrt(rho)>rho. Its counts are consequently exponentially smaller than the ordinary counts, up to fixed polynomial factors. Non-self-converse orbits contain precisely two ordinary classes; all stated scalar statistics are converse invariant. The uniform merged-node distribution and ordinary distribution consequently differ by exponentially small total variation after projecting to those statistics.

For feedback, the compiler also gives f≤n/3+d, so exponentially small mass contributes an exponentially small amount to first and second moments. Thus the O(1) moment conclusions, not merely the weak limits, pass to the merged model. The same conclusion could be obtained from the crude general polynomial feedback bound.

## 6. Exact finite controls and their limits

`check_slack_transfer.py` passed **1,809 exact assertions** and **96 numerical residual sanity assertions**. Standard and optimized Python outputs are byte-identical (`normal.json`, `optimized.json`). The script does not use assert statements, so optimized mode cannot silently remove checks.

Exact tests include:

- Independent bivariate formal rational division versus direct word-product convolution through order 48 for three explicitly synthetic finite compiler models. These arbitrary configurations are **not** asserted to be any actual slack inventory.
- A synthetic self-converse model with different maximal pole orders on its two parities; another with unbounded even support but only a finite odd exception.
- Literal palindrome enumeration through order 16, including the three different center-specific half-word cutoffs.
- Exact Euler-derivative moment recurrences, checked against complete small bivariate coefficient rows.
- The **actual inherited** delta=0,3,4 rational functions in ordinary and self-converse channels, checked through order 48 against their word descriptions. Their univariate coefficients through order 7 also agree with the existing independently generated tournament census. This review did not rerun or duplicate the compiler's full tournament census.
- Direct brute-force automorphism and feedback enumeration for the actual six-vertex diagonal-group example above.
- The exact bounded-perturbation variance counterexample on square n.

The 96 numerical sanity assertions evaluate exact finite moments at n=120/121, 360/361, and 1080/1081 and compare residuals with deliberately loose fixed thresholds. They are numerical diagnostics, not exact assertions and not proofs of limiting distributions. The proof is the moving-pole argument in Sections 1–5.

## Conclusion

No gap found in the conditional transfer theorem. The ordinary and supported-parity constants, actual group factor, converse merging, and bounded-perturbation caveat are all correct. Once the independently accepted compiler is pinned as a proved dependency, the transfer may be promoted from CONDITIONAL to PROVED. The optional logAut moment strengthening requires its explicit analytic-marker argument, not merely a change of wording.

## Final promotion verification

The final transfer source is now **PROVED**, SHA-256 `18c9bf4e131834de92b6dc72911729cd7cd82b0735638129572b790ec3e365d2`. Its compiler pin remains `22cb6aca899d9a485b2432d37a90d13040dfc6353ee41529a243d8489aa974a1` and matches the compiler file and its independent PASS receipt. Reversing only the status/dependency promotion and the added explicit palindrome-half formula reproduces the initially audited source hash exactly. There are no other changes. The inserted formula is precisely the center-specific tail formula proved and exactly tested above. Final verdict: **PASS**, with no outstanding conditional proof dependency and no required repair. The optional logAut moment strengthening remains a review corollary and is not claimed as an addition to the main text.
