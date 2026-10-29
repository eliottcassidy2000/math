# A parameter-selective valuation-clock target and its integer-source guard

**Status: PROVED elementary identities and conditional/finite certificates.**
No proof of Collatz, no claim that Haar-almost-every code covers a positive integer.
These are new packaging/bridges in this session, not claims of literature novelty.

## Inheritance and live concepts

Inherited: Session99 `transfer/INTERCEPT_MOMENT_BOUNDARY.md` keeps original-source pullback and distinguishes bounded slopes from heavy intercepts; Session97 `controllers/SOURCE_PRESERVING_RETRIES.md` keeps `(n0,n,J,D,C)` and warns that local descent is not original-source descent. The hostile controls are the positive 5x+1 cycles and arbitrarily long valuation-one runs. The present concepts are valuation coding, real/2-adic simultaneous sums, fractional slope moments, and source-marked clock barriers.

## 1. Exact code and simultaneous real/2-adic target

Fix an odd integer b>=3. On positive odd integers set

    T_b(n)=(bn+1)/2^a, a=v_2(bn+1)>=1.

For its valuation word `(a_1,a_2,...)`, put S_0=0, S_k=sum_{i=1}^k a_i and

    R_k=sum_{j=1}^k 2^{S_{j-1}}/b^j,
    U_k=b^k/2^{S_k}.

Exact iteration gives

    n_k=U_k(n_0+R_k).                                      (1)

The same code has 2-adic source

    Phi_b(a)=-sum_{j>=1}2^{S_{j-1}}/b^j in Z_2^odd.         (2)

The series converges because S_j>=j. It is a bijection from infinite positive-integer valuation words to odd 2-adic points whose iterates never encounter `-1/b` (where valuation is infinite). Proof: the first term makes Phi odd, and

    b Phi_b(a)+1=2^{a_1} Phi_b(shift a),

with the right-hand residual odd. Thus the exact valuation is a_1 and the map intertwines T_b with shift. Conversely iterating the affine identity determines the series. Excluded finite-termination points are exactly the countable backward orbit of -1/b. Positive integer orbits are never excluded.

Keep BOTH (2) and the ordinary real partial sums R_k. This is a faithful coding only while the original integer-source guard `Phi_b(a)=n_0 in Z_2` is retained. Reading the real value of (2) as that same integer would be invalid: the two completions have different notions of limit.

If S_k/k -> ell, then:

* 2^ell>b implies R_k -> infinity and U_k ->0 exponentially;
* 2^ell<b implies R_infinity<infinity and U_k ->infinity exponentially.

This follows by bounding the individual positive summands by geometric sequences with exponents ell +/- epsilon. In the second case (1) gives n_k ->infinity for every positive integer source with this code. In the first case U_k R_k need not tend to zero; multiplying a vanishing slope by a diverging carry is precisely the lost-information hazard.

## 2. Why the random target distinguishes 3 from 5

Under normalized odd Haar measure, valuations are independent with

    P(a=j)=2^{-j}, j>=1, E a=2.

Indeed a specified prefix of total S_k selects one residue class among the odd classes modulo 2^{S_k+1}, hence has normalized mass 2^{-S_k}. The prefix probabilities factor and determine the product law. The strong law therefore gives ell=2 almost surely. The real boundary in section1 is b=4, between b=3 and b=5.

An equivalent fractional-slope test is

    E[(b/2^a)^r]=b^r/(2^{r+1}-1), r>0.                   (3)

For 1<b<4 there is a sufficiently small r>0 making (3)<1: its logarithmic derivative at r=0 is log(b/4)<0. For b>=4 there is no such r: strict convexity of the log moment-generating function gives log E[M^r]>r E log M>=0. In particular for b=3 all 0<r<1 work (strict convexity between its zeros at0 and1); b=5 fails for every r>0. This is slope contraction in a random model, not an integer convergence result.

## 3. The source-marked clock barrier

**Finite certificate.** Let n_0 be positive odd and N<=n_0 positive. If n_i>=N for i=0,...,k-1, then

    n_k/n_0 <= (b+1/N)^k / 2^{S_k}.                       (4)

Proof: multiply n_{i+1}/n_i=(b+1/n_i)/2^{a_{i+1}}. Consequently, if

    2^{S_k} > (b+1/n_0)^k,

there is some j<=k with n_j<n_0. Otherwise every n_i>=n_0 and (4) contradicts n_k>=n_0. This certifies an earlier source descent, NOT necessarily descent of the displayed endpoint if an earlier descent already happened.

More generally, if an orbit stays >=N forever, then

    limsup S_k/k <= log_2(b+1/N).                         (5)

Use (4) and n_k>=N, divide logarithms by k, and let k grow.

For b=3 this gives the following EXACT reformulation for each positive odd source:

    orbit reaches1
      iff limsup S_k/k > log_2(10/3).                     (6)

Forward: after reaching1 all valuations are2, so the limiting mean is2. Reverse: if1 is never reached, all odd states are>=3, and (5) contradicts the strict inequality. Thus a code-based proof only needs a source-restricted clock lower bound above ~1.737, not full random genericity (mean2). The strict boundary is not asserted sharp.

For b=5, in contrast, a source code with mean2 diverges by section1. The known 5-cycle 1->3->1 has word(1,4), mean5/2; the other cycle 13->33->83->13 has word(1,1,5), mean7/3. These are explicit hostile controls against transporting the 3-parameter conclusion merely from a shared shift space.

## 4. A finite distributional certificate, and the exact missing bridge

Fix odd N>=3 and b=3. Put c=log_2(3+1/N)<2. Any source n>=N with no descent below its own value during k steps has S_k<=ck, by section3. For every t>0, normalized Haar mass of this clock event is at most

    [ 2^{t c}/(2^{t+1}-1) ]^k.                           (7)

This is Markov's inequality applied to 2^{-t S_k}; the bracket is <1 for sufficiently small t since its log derivative at0 is (c-2)log2<0. Each finite clock event is a finite union of exact source cylinders because all its words have S_k<=floor(ck), so (7) also bounds the natural density of positive integer sources >=N failing descent for k steps. For b=5, c>2 and this particular exponentially decaying bound is unavailable.

This is a density statement, not an empty-set statement. Every specific infinite code is Haar-null, including every positive-integer source code. Arbitrarily small finite-prefix failure mass can contain a particular nested sequence of cylinders forever. Even a source-restricted exceptional set of natural density zero could contain an integer. The missing theorem is precisely exclusion of the positive-integer image from the bad-clock set in (6), not additional mixing of the ambient shift. Dropping the exact source guard in (2) changes the theorem.

## 5. Concrete next obligation

The target which genuinely discriminates 3/5 is the integer-source clock predicate in (6), with (1) as its carry guard. A possible noncircular advance would prove that a class of 3-adic/2-adic arithmetic source constraints forces limsup S_k/k above log_2(10/3), or produces a finite source descent certificate. It must identify that class explicitly and exclude the 5-cycle controls. The bijection itself does not provide this arithmetic assertion.

### Certified-basin refinement (same proof, no heuristic)

If a finite exact certificate establishes that every positive odd m<N reaches1, replace3 by N in (6):

    reaches1 iff limsup S_k/k > log_2(3+1/N),

provided N is odd and N>=3. The forward direction still has limiting mean2; a nonconvergent orbit avoids the certified basin and remains>=N. This can move the threshold arbitrarily close to log_2(3) only conditional on correspondingly large finite basin certificates; it never turns a finite certificate into the missing infinite one. The verifier's finite orbit samples are not asserted to be such a complete basin certificate.

Reproduction: `python3 -B check_drift.py` and `python3 -B -O check_drift.py`. These perform 180180 exact affine/source-prefix checks, 63174 triggered source-barrier checks,18 unique-cylinder checks, and the two 5-cycle hostile controls. All pass. The elementary infinite proofs are above; finite tests do not prove their quantifiers.
