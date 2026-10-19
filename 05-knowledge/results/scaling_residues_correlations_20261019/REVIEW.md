# Independent review: scaling, residues, and arithmetic transfer

19 October 2026. **Verdict: APPROVED as a mathematical packet, with its stated inherited inputs.** No substantive correction is required. The frozen author files were not changed.

Reviewed manuscript: `SCALING_RESIDUES_AND_CORRELATIONS.md`, SHA-256
`7a48b6458f3bd0954eb4c9309093ff3c32f8fbfae61eba8618319cd66e062797`.

## Main findings

1. **Cusp-window transfer is genuinely relative and uniform in the statistic.** For fixed escape radius, the approved threshold phase gives `w_n ~ 2 pi^2 n^-3`. On the stated window, the manuscript's `K_(M,L)` bounds the largest weight divided by the smallest. Thus every nonnegative statistic, even one supported on a single level or depending on M and L, has weighted/unweighted average ratio in `[1/K,K]`. Zero statistics transfer exactly. No lower bound on event density and no independence assumption is needed. Positive integer L may remain bounded; `M -> infinity` and `L=o(M)` suffice. The strict escape convention and the half-open parameter cells agree.

2. **The dimension-four error survives the coordinate change unchanged.** Write `q_m=R(m;J,z)/(L V_m(z))` and let `t_m` be the weighted/unweighted event-density ratio. Then `t_m in [1/K, K]`, and the parameter-space normalized count is exactly `t_m q_m`. If `|q_m-1|<=e`, its error is at most `e+(K-1)(1+e)`. Taking the inherited uniform short-interval limit and `K -> 1` proves (8) with precisely `exp(37-u)`. The theorem's `L=floor(m^theta)` tends to infinity, is `o(M)`, and gives an admissible integer interval. Positivity follows from `u>37`; the lower bound uses only the uniform lower bound for V. The strict condition `r>u/theta` correctly absorbs the factor 6 and yields at most `r-1` prime factors in each neighbor.

   The primary [Thorner–Zaman source, Lemma 6.2](https://msp.org/ant/2019/13-5/ant-v13-n5-p02-s.pdf), pp. 1058–1060, was checked against the inherited specialization. At dimension four its threshold is 37 and its exponential term is `exp(36-s) K^10`; the repository deliberately weakens this to `exp(37-s) K^10`. The fixed pre-sieve, remainder normalization, and order of limits in the inherited short-interval proof retain uniformity in the target and interval location.

3. **The prefix model and its exact total variation are correct.** Uniform input gives independent uniform digit blocks and the stated translated-interval transition kernel. Uniform stationarity holds even when b and W are not coprime. The row total variation is exactly `t(W-t)/(WL)`, including `L<W`, `W|L`, and the initial time `n=0`. The correlation estimate correctly pays the initial-distribution error plus `bar(f)` times the transition error. The normality boundary is stated accurately and is not used to infer anything about pi.

4. **The arbitrary-gcd Fourier statement includes the necessary nilpotent pieces.** On mean-zero functions let `U_0` be the uniform-mean projection of composition by multiplication by b. The identity `bPD=DU_0` and the invertibility of D on that space give `bP = D U_0 D^-1` there. Fourier arrows are `j -> bj mod W`, with arrows into zero killed. Periodic orbits yield roots of unity; the other directions yield the zero-eigenvalue Jordan part. In particular exact uniformity occurs precisely when `W|b^h`, consistent with the total-variation formula. Independent characteristic-polynomial and power-rank checks cover both the spectral and nilpotent assertions.

5. **Affine local factors and the Boolean-gauge iff pass.** The affine forbidden-root union has exactly the three stated cases, including multipliers divisible by p and the collision cases `L=+/-1`. CRT applies to the prescribed squarefree wheel without requiring the multiplier to be invertible. For Boolean gauges, taking logarithms makes every diamond a vanishing mixed difference; the resulting increment is independent of the other coordinates, giving exactly a constant plus singleton contributions. This proves necessity as well as sufficiency, including Boolean ranks zero and one. The antichain mantissa example and its factor-ten repair are exact.

## Independent finite controls

The copied author verifier passes **148,202 active guards** under both ordinary and optimized Python. Both receipts are byte-identical to the frozen receipt.

The separately written `independent_checks.py` passes **178,019 active guards** under both ordinary and optimized Python, with identical receipts. It does not import the author's verifier. Its checks include:

- Integer characteristic polynomials for all `b=2..10`, `W=2..32`; exact power ranks and nilpotence for a stated sub-universe.
- Arbitrary interval lengths `L=1..74`, moduli `W=2..40`, and events attaining the exact total variation.
- Every multiplier and shift modulo each prime from 5 through 101.
- Exact rank/nullity of all Boolean diamond equations through rank six, plus nonmodular hostile controls.
- Nonmonotone positive rational weights satisfying the stated cubic envelope, all statistics in `{0,1,5}` on windows up to length six, and exact relative-error composition.

Reproduce from this directory with `python independent_checks.py` and `python -O independent_checks.py`. These checks verify finite algebra, not the analytic cusp phase, sieve asymptotic onset, or any prime-distribution assertion. The phase-width implication was reviewed using exactly the approved `o(1)` threshold phase, without assuming a quantitative phase rate.

## Optional editorial improvement

For the particular power-sized interval used in Section 3, the inherited full-singular-series tail criterion already holds: `log log(9m^3)/log z -> 0`. Thus the same corollary can also be normalized by `e^(-4 gamma) S_4(m)/(log z)^4`, with unchanged limiting relative-error envelope. Keeping V is correct and more general; the manuscript's tail warning remains appropriate for other window scales. This is an optional strengthening, not a required repair.

The sibling phase hyperlink resolves to the accompanying proof.
