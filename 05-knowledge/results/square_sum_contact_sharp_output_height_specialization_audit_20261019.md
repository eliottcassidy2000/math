# Independent specialization and output-height audit

**Verdict: PASS.** No blocking error was found in the canonical-height argument or its transfer to primitive output length. This audit treats the established cusp-index lower bound and recurrence theorem as inputs; it does not broaden the result beyond the prescribed plus-port orbit.

## Primary sources and hypotheses

1. [Silverman, *A Lehmer-Type Lower Bound for the Canonical Height on Elliptic Curves over Function Fields*, Theorem 2, pp. 1–2; Definition 7, p. 5](https://arxiv.org/pdf/2402.14771). The theorem permits an algebraically closed constant field of characteristic different from 2 and 3, a smooth projective function-field base, extension degree one, and a non-torsion point on a nonisotrivial elliptic curve. These hypotheses hold for `Qbar(t)`: the stated invariant has a genuine pole at `t=1`, and `P` cannot be torsion because its smooth specialization at `t=5/4` is non-torsion. The source's displayed canonical-height convention differs from the note's half-x-height convention by a positive factor, irrelevant to positivity.

2. [Call–Silverman, *Canonical heights on varieties with morphisms*, §§3–4 and Theorem 4.1, p. 183](https://www.numdam.org/item/CM_1993__89_2_163_0.pdf#page=22). The limit is over all good algebraic base points with degree-one base height tending to infinity. No real compactness condition occurs. A direct application uses the smooth projective surface `P^1_X × P^1_t`, the fiber-preserving duplication rational map, the first-factor divisor class, and section `X=0`. Its generic degree is four; the pullback class minus four times that class is fibral. On the smooth good locus the fiber maps are morphisms. The dynamical canonical height is twice the note's elliptic height on both generic and specialized fibers, so the factor cancels. The function-field height is the geometric degree height over `Qbar(t)`. Thus both rational sequences, with heights `log(j+2)` and `log(2j+2)`, are covered even though their real limits are bad fibers.

3. [Silverman's 2006 MSRI lectures, slide 22, PDF p. 27](https://legacy.slmath.org/attachments/workshops/301/HtSurveyMSRIJan06.pdf#page=27) independently display the same unrestricted high-height specialization limit for a fixed section.

The test specialization is exact: `k(5/4)=49/41`, `f(5/4)=720/1681`, and `K` is neither 1 nor −1. On the integral model from the predecessor, the distinct successive doubles have x-coordinate valuations `−2n`; this proves infinite order. A generic torsion relation specializes at that smooth fiber, giving the required contradiction.

## Uniform estimates

- The duplication numerator has coefficient absolute sum at most `16 M²`, and the denominator is at most `16 M² H⁴`. The numerator is nonzero at each denominator root when `K²≠1`; infinity is also valid. Reduction only decreases height. Iterating therefore gives exactly the one-sided comparison (11), including its factor `1/6`.
- The stated bounds for `h(k)`, `h(K)`, and `h(X(Q))` give `hat h(Q)≤(7/3)h(t)+(1/2)log 96≤5h(t)` for `h(t)≥1`.
- The reverse triangle inequality is applied with its positive part before squaring. For every integer `q` above the fixed threshold, it yields `hat h(Q+qP)≥(a/4)q²h(t)`. Combining with the one-sided comparison gives (15) uniformly in `q`; specialization is used only for the fixed section `P`.

## Forward transfer and all-success quantifier

The forward quartic map has the correct direction for the lower bound. In fact, clearing common denominators gives the explicit refinements

`h(z) ≤ 2h(u)+h(k)+(1/2)log 12`,

`h(X) ≤ 2h(k)+h(z)+2h(u)+log 5`.

Hence `h(X)≤4h(u)+6h(t)+O(1)` with absolute constants. These estimates do not require any inverse-map denominator to stay away from zero.

For `v=u` or `v=−1/u`, rational height is unchanged. With `v=n/m>1` reduced, `h(v)=log n`, and division by the primitive normalization factor 1 or 2 gives `|log C′−2h(u)|≤log 2`. Thus the quadratic-in-index lower bound survives reduction, either chart, and leg exchange. The cusp-index bound makes its subtractive `O(h(t)+1)` term absorbable for **every** successful index. It yields `log C′≥c q² log j≥c₁(log j)³`.

The recurrence supplies an `O(log j)` success. Therefore the minimum output length and every least-absolute-index success have matching cubic upper bounds. Larger successful indices are not claimed to have cubic upper bounds.

## Optional clarification for the inherited upper bound

The upper-bound predecessor skips `X=0,K` and infinity. No successful finite quartic point is lost here. Substitution of the forward formula into the quartic gives:

- `X=K` forces `u²(K²−1)=0`, impossible for successful `u≠0`.
- `X=0` forces `u=−1/(2K)`. For `1≤K≤3/2`, its spread is `8K(4K²−1)/(4K²+1)²≥24/25`, hence inadmissible. The inequality follows from
  `f(−1/(2K))−24/25 = −8(K−1)(2K−3)(4K+1)(6K+1)/(25(4K²+1)²)≥0`.
- The forward x-coordinate is finite whenever the accepted `u` is finite and nonzero.

Adding these observations to §7 would make the applicability of the inherited inverse-map upper estimate fully explicit; they are a clarification, not a counterexample or unresolved assumption.

No main proof files were edited and no GitHub changes were made during this audit.
