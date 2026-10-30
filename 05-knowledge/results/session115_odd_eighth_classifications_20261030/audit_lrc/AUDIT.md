# Independent audit: harmonic completion and odd-eighth weights

**PASS for the stated scoped results.** Independent algebraic review and a separately written exact checker. No unrestricted LRC, lower-chamber classification, or literature-priority claim.

## Proof review

The inherited core has a complete top-window atlas above1/(q+2). In each direction the owner e gives r≤(1−ey)/q for y=q|x|. The first tail bounds r≤my, while the harmonic block bounds r≤1/(k+1). Choosing y=r/m in a least-owner direction proves the claimed minimum formula whenever both candidate bounds lie in the inherited chamber. If the first bound is strict, equality forces exactly the two least-owner endpoints. The restriction2m>s and1≤k≤q is essential to this argument.

For a full q-tail block, attaining its optimum1/(q+1) requires its base phase to be a primitive multiple of1/(q+1). Pulling these phases back to each named core direction gives precisely the positive indices j with ej≤m and gcd(j,q+1)=1. The smallest possible displacement proves the necessity m≥s, and e=s,j=1 proves sufficiency. Distinct core windows and the separation of left/right sides prove that the reflected list has no duplicate points. Thus the census and its Möbius floor form are exact.

At a listed maximizer the two tail owners have indices inverse(j) and q+1−inverse(j). Any additional core owner is slower than both tail owners, including at ej=m. Therefore the literal maximum-owner rule gives the displayed reciprocal-weight formula. At the collision m=s this correctly gives2(q+1)/(q²s), rather than continuing the strict-branch two-owner formula. In particular q=7,m=2 gives8/49.

The affine-period laws follow because L/e is an integer multiple of q+1 for every e when L=(q+1)lcm(E_q). A full unit period contributes φ(q+1) to the count and twice the reciprocal-unit sum to the weight. The asymptotic constants and bounded/reciprocal error terms follow without an equidistribution assumption.

For q=7, the inverse-residue weights are8/7 and8/15 before dividing by7m. Their exact character decomposition is88/105+(32/105)χ₂(j). The odd-prefix χ₂ sums have partial values0,±1 and period8 in the endpoint. Hence the aggregate correction is bounded by6 and has period960. The stated count-weight identity and the real eighth-root cosine identity are valid on the actual four odd phase labels. They are not a transport of LRC dynamics to another mathematical system.

## Independent exact tests

`independent_harmonic.py` imports no producer implementation. It enumerates the entire exact maximizing grid of the harmonic block, then filters with literal core inequalities and computes actual positive/negative endpoint owners from all speeds. This checks390 rows with q=3..15,m=1..30, including both sides of the preservation threshold.

It also checks complete maxima from every opposite-slope rational candidate for94 small high-chamber rows, exact long-period count/weight identities for21 rows, and99 odd-eighth character/period identities. The receipt records the actual gate count. Normal and optimized outputs agree byte for byte; gates use explicit exceptions.

Run `python independent_harmonic.py` and `python -O independent_harmonic.py` from this directory. Python3 standard library only. Finite controls corroborate the all-parameter written proof; they do not extend its low-clearance scope.
