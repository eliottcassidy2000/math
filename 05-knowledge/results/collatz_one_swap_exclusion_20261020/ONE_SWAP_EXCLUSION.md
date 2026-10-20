# Excluding positive integer cycles in the primitive one-swap mechanical family

20 October 2026. **PROVED for the stated primitive one-swap family, using the cited Matveev theorem and independently audited exact certificates.** The statement concerns a specified family of parity words, not the general Collatz conjecture. The sole imported analytic input is the explicit real Matveev bound in the version cited below. No published Collatz convergence verification or minimum-cycle computation is required.

## 1. Precise theorem

Let q>=4, 0<p<q, gcd(p,q)=1, and D=2^q-3^p>1. Set

    c_i=floor((i+1)p/q)-floor(ip/q),  0<=i<q.

Let w be any word obtained by swapping cyclically adjacent unequal bits of c. With chronological carry

    B(w)=sum_(i:w_i=1) 2^i 3^(number of ones after i),

we prove

    gcd(D,B(w))<D.

Rotations, swaps across the original cut, repeated descriptions of one output, and outputs returning to the balanced necklace are all included. Coprimality and the one-swap restriction are part of the theorem’s stated scope. This closes the precise one-swap positive-cycle exclusion formulated in the predecessor; it does not exclude arbitrary cycles, repeated nonprimitive shapes, multiple-edit families, or divergent trajectories.

For p=1, B is a power of two, so gcd(D,B)=1<D directly. In the rest assume p>=2. If D divides B, the associated rational fixed cycle is integral and positive: all rotations have the same cancellation denominator, and the parity of each integer entry is forced by integrality of its next shortcut branch. Coprimality excludes a repetition of the ordinary one-odd-step cycle.

## 2. Inputs and what each supplies

1. [Sharp one-swap content bound](../collatz_one_swap_sparse_norms_20261020/SHARP_ONE_SWAP_CONTENT_BOUND.md): g=gcd(D,B)<=3^p-2. Hence a hypothetical positive integer cycle must have q=ceil(p log_2 3).
2. The separately audited [legal-edit refinement](../collatz_one_swap_critical_gap_20261020/CRITICAL_CLOCK_SQRT7_GAP.md) proves the stronger critical-shape gate

       g < 7^(p/2).

   Its mechanism is a sparse norm and the fact that legal edits occupy restricted exponent ranges. It uses root-of-unity orthogonality and the arithmetic-geometric mean inequality, not logarithmic Diophantine approximation.
3. The explicit Matveev theorem, as stated and derived in Bugeaud, Mignotte and Siksek, *Classical and modular approaches to exponential Diophantine equations I*, arXiv:math/0403046, Theorem 9.4, printed page 16: for real positive algebraic numbers, the nonzero multiplicative form has the bound described in Section 3. This is a **CITED** theorem, not proved or claimed new here.
4. Exact rational interval arithmetic and a complete finite check of the remaining small shapes, supplied in `certify_one_swap_closure.py`.

The full argument retains the same p,q and edited word throughout. A resultant cap alone would lose that common witness; the clock and analytic approximation restore the needed constraints.

## 3. Explicit finite upper bound from Matveev

Put

    Lambda = 2^q 3^(-p)-1 = D/3^p >0,
    lambda = q log 2-p log 3 = log(1+Lambda)>0.

The critical-shape content gate and D=g give

    Lambda < (7/9)^(p/2) < exp(-p/9),
    0<lambda<Lambda<exp(-p/9).                         (3.1)

The strict second inequality follows from log(1-2/9)<-2/9.

Use the real-field form of Theorem 9.4 with n=2, number-field degree 1, alpha_1=2, alpha_2=3, b_1=q, b_2=-p, B=q, A_1=7/10 and A_2=11/10. These A values exceed log2, log3 and 0.16, and equal the corresponding required height upper bounds. The theorem yields

    log Lambda > -C(1+log q),
    C=1.4 * 30^5 * 2^(9/2) * (7/10)*(11/10) < 10^9.     (3.2)

Indeed replacing 2^(9/2) by23 gives the exact rational upper value 602,494,200. The choices of heights and the inequalities are checked by the rational certificate.

Every positive integer cycle has q<=2p, by multiplying its odd-step ratios 3+1/n_i<=4. Here gcd(p,q)=1 and p>=2 exclude equality, so q<2p. Combining (3.1) and (3.2),

    p < 9*10^9 (1+log(2p)).                             (3.3)

The function p/(1+log(2p)) increases for p>=1. At p=10^12, log(2p)<30 and 9*10^9*31<10^12. Thus

    p<10^12.                                           (3.4)

This conservative bound is deliberate; no optimized or recalled transcendence constant is needed.

## 4. Continued fractions exclude every p>=100

Let alpha=log3/log2. For p>=100, (3.1) gives

    0<q/p-alpha<exp(-p/9)/(p log2)<1/(2p^2).            (4.1)

For the last inequality it is enough that 2p exp(-p/9)<log2. The left side decreases for p>9, and at 100 the inequality is certified by a positive rational partial sum for exp(100/9) and a rational lower bound for log2.

Since gcd(p,q)=1, Legendre's elementary continued-fraction criterion makes q/p a convergent of alpha. The certificate enumerates every such convergent with100<=p<10^12. There are 18. Those below alpha have negative lambda and cannot be positive cycles. Every remaining convergent has a certified positive lower bound

    lambda > 2^(-min(floor(p/9),256)) >= exp(-p/9),      (4.2)

contradicting (3.1). The first possible upper convergent in this range is 485/306; no assumption about a numerical continued-fraction prefix is used without interval certification.

The terminal convergent is 8573543875303/5409303924479, whose denominator exceeds10^12. Convergent denominators increase thereafter, so the range is exhausted. Intermediate fractions are not omitted candidates: inequality (4.1) specifically forces ordinary convergents, rather than merely good rational approximations.

## 5. Direct finite closure for 2<=p<100

Only the unique critical length q=ceil(p alpha) needs testing. Exact rational log bounds certify that integer for every p. Noncoprime shapes are excluded by the theorem's hypothesis. There are 59 remaining shapes, and the program examines every one of their 3,548 cyclic unequal-adjacent-swap occurrences, retaining duplicates and cut-crossing swaps.

For each actual word it computes the carry twice: by direct affine composition and by the positional sum. Both values agree, and every gcd is strictly smaller than D. This is a complete finite check for the stated residual universe, not a random sample or a check of arbitrary parity necklaces.

There is no overlooked D=1 case at p>=2. If 2^q-3^p=1, reduction modulo 3 makes q even. Then the coprime odd factors 2^(q/2)-1 and2^(q/2)+1 multiply to 3^p; they must be 1 and 3, forcing (q,p)=(2,1). Thus every positive critical clock in this residual primitive universe has D>1, as the verifier checks directly.

Sections 3-5 exhaust all p>=2, while the p=1 case was elementary. The separately proved and audited critical-shape gate applies, so the displayed theorem follows.

## 6. Exact certificate, not floating-point evidence

For t=1/3 and1/2 the checker uses

    log((1+t)/(1-t))=2 sum_(j>=0) t^(2j+1)/(2j+1).

After N=220 terms the omitted tail is positive and bounded above by

    2t^(2N+1)/((2N+1)(1-t^2)).

Every bound is a rational number. The ratio interval for alpha is formed using the lower numerator/upper denominator and vice versa. A continued-fraction digit is accepted only when both interval endpoints have the same floor; reciprocation reverses the endpoints. Gap signs and all exclusion inequalities are rational comparisons. No decimal approximation is used for a decision.

Normal and optimized Python runs produce identical JSON, with 7,302 explicit checks. The script imports no repository mathematical code. The finite proof records the exact intervals, all 18 large convergents, all 98 small p-values, which 59 shapes survive coprimality, and every shape's swap count and maximum content.

Reproduce:

    python3 certify_one_swap_closure.py > closure.replay.json
    python3 -O certify_one_swap_closure.py > closure.optimized.replay.json
    cmp closure.json closure.replay.json
    cmp closure.json closure.optimized.replay.json

The earlier rate 1/2500 draft and its optional Barina route are not dependencies and are excluded from the final publication inventory. The present proof needs neither a general finite Collatz convergence bound nor a numerical claim about irrationality measures of log2(3).

## 7. Boundaries and next question

The statement excludes a sharply described near-balanced cycle family. It neither implies that every hypothetical Collatz cycle belongs to that family nor controls nonperiodic trajectories. The next arithmetic question is whether several edits still admit an exponential norm base below 3 after retaining their legal relative positions. Merely multiplying independent one-edit bounds would not preserve the common modular root or prove such a bound.

No historical novelty is claimed for the classical logarithmic tools, the product identity, Legendre's criterion, or the numerical continued-fraction prefix. The contribution is the explicit connection from the one-edit sparse norm to a complete, reproducible exclusion of this particular family.

## 8. General mechanical shapes and independent audit

The [repetition-factor obstruction](REPETITION_FACTOR_OBSTRUCTION.md) shows that one adjacent unequal swap of any repeated nonconstant binary word destroys the entire repetition multiplier of its carry. In particular, every noncoprime mechanical shape is excluded for either clock sign. Combined with the theorem above and the elementary p=1 calculation, this gives:

**General positive-family corollary.** For every q>=2 and 0<p<q, one cyclic adjacent unequal swap of the zero-intercept mechanical word cannot encode a positive integer Collatz cycle, except for the ordinary one-odd-step case (q,p)=(2,1). No coprimality restriction is needed in this corollary.

This is still a special word-family theorem, not a claim that arbitrary Collatz words are near-balanced. The repeated-word lemma does not require the mechanical assumption and is proved independently in its linked note.

The [independent closure audit](audit/DIOPHANTINE_CLOSURE_AUDIT.md) checked the primary Matveev statement, its exact application, the proof's quantifiers, and the complete range of the continued-fraction certificate. A separately written verifier uses a different logarithm series and a right-to-left carry calculation. It reproduces all 18 large convergents (nine below and nine above the limiting slope), the terminal convergent, and all 59 small shapes / 3,548 swaps. The proof-status text, the now-published dependency link, and this separately audited repetition corollary were finalized after that audit; the audited closure script and its exact receipts are unchanged.
