# Independent audit: amicable deficiency bridge

2026-11-04. Mathematical result: PASS.

Checked the Euler deficiency equivalence in both directions, including D>0, positivity of the two factors, d<e, prime and coprimality guards, and distinctness N<M. Checked the mod-12 sieve under its stated p,q>3 hypothesis. No mathematical error found.

Ran verify_amicable.py normally and with -O. Outputs match byte for byte: 10 certificates, 4 odd kernels, and 9,975 small prime-pair cases. Independently recomputed sigma by paired divisor enumeration, rather than the script's prime-factorization implementation, for all ten emitted certificates. All pass. The three worked examples have common divisor sums 504, 262080, and 157248 respectively.

Primary-source checks passed: Euler's translated §27 explicitly has (cx-b)(cy-b)=b² with b/c=a/(2a-sigma(a)); the Pomerance manuscript is dated March 30, 2014 and states infinitude remains unproved; Pollack's October 2026 introduction states the same and gives the stronger upper bound. The cited Chen–Tang–Zhan version really does state the false odd–odd nonexistence conjecture in §3.8. The explicit odd certificates disprove that assertion without relying on another authority.

Code robustness recommendation delivered to author: replace assert-based validation by explicit checks, since Python -O disables assertions. Matching optimized output by itself is not an independent validation. The small-prime comparison uses an independent enumeration but the same sigma engine; characterize its independence accordingly or add a separate divisor-enumeration engine. These do not invalidate the mathematical proof or independently replayed certificates.

Sources inspected:
- https://arxiv.org/html/2409.08783v2 (§27)
- https://math.dartmouth.edu/~carlp/amicable032914.pdf (introduction)
- https://arxiv.org/html/2610.05661v1 (introduction)
- https://arxiv.org/html/2601.07444v1 (§3.8)

Final clearance: author replaced assertion statements with an explicit require function and added independent divisor-pair sigma checks for every certificate. The hardened script was re-inspected and rerun both normally and with -O; outputs still match. The robustness recommendation is resolved.
