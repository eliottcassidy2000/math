# A finite arithmetic certificate for all periodic failure itineraries

**Status: PROVED, all periods.** This excludes periodic and eventually periodic integer exceptions in the five-word retry language. It does not exclude aperiodic positive integer exceptions.

Let the five failure charts be `f_i(x)=lambda_i x+beta_i`, with the exact source guards from `FIVE_EXCEPTIONAL_CYLINDERS.md`. Every lambda_i>1 and beta_i>0. Their individual real fixed points are

    alpha = (−1, −179/17, −16177/139, −12289/139, −9697/139).

## Convex-hull fixed-point lemma

For any finite composition of affine maps with lambda_i>1, its fixed point is a convex combination of their individual fixed points. Indeed write

    f_i(x)=lambda_i x+(1−lambda_i)alpha_i.

For a chronological word i_1,...,i_k, the fixed point is

    sum_j [(lambda_ij−1) product_{l>j}lambda_il] alpha_ij
    --------------------------------------------------.
                    product_j lambda_ij−1

All weights are positive, and their sum is the denominator by telescoping. Thus every fixed point lies in

    I=[−16177/139,−1].

Every point of a periodic orbit is the fixed point of a rotated word, so the entire orbit lies in I. Consequently an integer periodic orbit uses only the58 odd integers−115,−113,...,−1. This is an all-period reduction, not a bounded-period search.

## Exact guard certificate on that interval

On these58 integers:

- The first chart is n≡3 mod4, with map(3n+1)/2. It strictly decreases every n<−1 and fixes−1.
- The second chart n≡125 mod128 occurs only at n=−3, and sends−3 to−1.
- None of the final three charts occurs: their largest negative representatives are−2739,−1827,−1219, all below I.

Therefore **the only integer periodic cycle in the entire failure language is the fixed point−1**. The finite exception−3 enters it but is not itself periodic. This proof uses source guards as well as affine fixed-point locations; slopes or carries alone would leave additional rational candidates.

## Eventually periodic words also force actual periodic states

An infinite repeated finite failure word cannot describe a nonperiodic integer orbit masquerading as a periodic label sequence. If its total affine map is `(A x+C)/2^D`, its rational fixed point is

    alpha=−C/(A−2^D),

with odd denominator and odd numerator. The fixed point realizes that exact finite valuation word by the source-cylinder formula. Realizing k repeats forces x and alpha into the same residue modulo2^(kD+1). Letting k grow in the2-adics gives x=alpha. Thus an eventually periodic failure itinerary reaches an actual rational periodic point. An integer itinerary must reach−1 by the preceding theorem.

Checking integer preimages through the five exact affine maps gives

    failure-preimages(−1)={−1,−3},
    failure-preimages(−3)=empty.

These are finite divisibility calculations; both surviving preimages satisfy their guards. Hence the only integer sources having an eventually periodic infinite failure itinerary are **−1 and−3**.

In particular every positive integer that stayed forever in the failure language would have to produce an **aperiodic, strictly increasing, unbounded** sequence of first-return values. The periodic class has been arithmetically eliminated, but this final aperiodic possibility remains OPEN. Haar measure zero does not eliminate it.

## Reproduction and transfer mechanism

`check_periodic_certificate.py` verifies the interval and all58 integer guard candidates, the integer preimages, and independently checks all3905 block words of lengths1 through5. The bounded word check is a hostile probe; the proof covers every period through the convex-hull argument. Normal and optimized Python outputs agree in `periodic_checks.json`.

The transferable mechanism is: affine expansion plus positive intercept places all periodic fixed points in a finite convex hull; integer discreteness then turns that hull into a finite arithmetic guard certificate. It can apply to other finite affine languages when their slopes share the same side of1 and their exact source guards remain available. It does not apply unchanged to mixed expanding/contracting alphabets or to unspecified quotient edges.
