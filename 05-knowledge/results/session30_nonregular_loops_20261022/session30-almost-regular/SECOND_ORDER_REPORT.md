# Session30 continuation: prospective second-order tests

Status: FINITE-EXACT H targets, double-precision predictions. Candidate decorated-loop engine was supplied by a separate lane; this is an independent-target validation, not an independent implementation of that engine. Its separate formal/core audit and the open uniform remainder must remain distinguished.

## Protocol

Fresh samples selected before outcomes: n12,16,20, three regular triangle-reversal-chain samples at each order n+1, fixed seed30600+n,1000 accepted reversals per successive sample, then delete vertex0. Exact H recomputed using this packet's integer endpointDP. The nine cases are held out from first baseline and formula construction. Eight additional original Paley/rotational/nonvertextransitive-block controls are explicitly retrospective. No coefficients fitted.

All 17 cases improved absolute y error from fullfirst order to fullsecond order. Fresh n20 samples:

|case|exacty|firstprediction|secondprediction|
|---|---:|---:|---:|
|0|1.385242548|1.446382833|1.385894360|
|1|1.410467647|1.471301313|1.412089311|
|2|1.462190649|1.520827594|1.463112773|

Thus errors~0.059–0.061 become0.00065–0.00162. At rotational deletion n20,error0.10547 becomes0.01489. Nonvertextransitive block examples atn20 improve0.07964→0.00203 and0.09367→0.00767.

All exact H, bit strings, intermediate predictor ingredients and timings are in second_order.jsonl. No universal rate is inferred from these finite samples. In y units, first order targets residualO(1/n); second order targetsO(1/n²), conditional on the still-unproved asymptotic expansion. Differences of signs in scaled second-order residuals already occur, so no fixed empirical scalar was inferred.

## Numerical stability control

For fresh n12 case0, rotational n20 andnonvertextransitiveblockn20, recompute the producer's Cauchy derivatives with half radius0.0125 and twicepoint count32, versusradius0.025/16. Results in derivative_stability.json. This is a numerical stability check, not interval certification.

## Exact invariant ladder and hostile rank example

n6: identical spectrum and identical ordered score vector can coexist with H41versusH45; SWITCHING_OBSTRUCTION.md proves an entire2^m-member almost-regular switching family on2mvertices.

n8: the witnesses in invariant_ladder.json have H653and651, equal characteristic polynomial and equal exact tree penalty11169728/10230017. But the figure-eight,theta,dumbbell andlollipop invariants all differ (exact rational arithmetic). Thus spectrumplusone tree-resolvent scalar still cannot determine H.

The richer predictors distinguish them but rank them wrongly:
- H653:exacty1.199192006,first1.271823510,second1.198275672.
- H651:exacty1.220021433,first1.253367447,second1.191368547.

Both firstandsecond predict the reversed y ordering. Secondorder still improves each absolute error, but restores neither guaranteed ranking nor exact reconstruction. Any proposed universal ranking theorem must exclude or overcome this explicit small-order obstruction.

## Genuine leading lollipop obstruction

The held-out non-vertex-transitive regular7 base, lexicographically blown up by rotational oddk, yields n=7k−1 after deletion. Direct numerical nL evaluations extend throughn440; no H computedbeyond20. Applying the theory lane's block-limit formula2mG_Bdiag(H_B−G_B) to this packet's independently reconstructed7x7base gives exactly:

- Exceptionalblock:12607/496782.
- Each otherblock:−12607/2980692.

The matrix and rational results are in lollipop_exact_limits.json. This independently checks the algebraic limit value; its limiting argument is supplied by the theory lane. It demonstrates why the old path-only baseline cannot be regarded as a complete first order formula: the omitted lollipop can have nonzero nL limit.

## Reproduction and dependencies

Run test_second_order.py, invariant_ladder.py, derivative_stability.py with the sibling session30-decorated-engine directory present; these explicitly import its predict.py. The primary exact target code never imports producer count code. switching_family.py andlollipop_scaling.py use only this packet's independently implemented base algorithms.
