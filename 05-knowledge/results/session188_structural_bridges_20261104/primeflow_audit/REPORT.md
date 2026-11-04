# Independent review of weighted prime-factor closure

**PASS, in the stated local scope.** The mathematical proof, exact examples, source checker, and independent bounded checker agree. Two wording clarifications were requested and incorporated before publication: the closure domain explicitly excludes positive demands at role-disallowed targets, and the Lucas recurrence now displays its minus sign. No theorem-level correction was needed.

## Proof review

For target r different from the source p, the order d divides r-1. Therefore r does not divide d. When d=1, cancellation of p-1 gives valuation v_r(n), not v_r(p-1)+v_r(n). When d>1, order and odd-prime LTE give zero off multiples of d and c+v_r(n) on them. This proves the exact rank divisibility formula for every positive length.

Role compatibility is preserved by the LCM: ordinary allowed orders are odd; Euler allowed orders have 2-adic valuation at most one and the base length supplies exactly one factor of two. Hence L(u)|n is equivalent to the lower demand for role-compatible n. Divisibility of geometric sums proves the least-completion property, with no claim of independent edge choice. Extensivity and this least property give idempotence and the exact fixed-point criterion.

The ordinary L=1 case needs a physical positive exponent, and the proof supplies it. Multiplying a realizing length by any sufficiently large odd prime changes neither the target-prime valuations of the length nor which target orders divide it. This establishes unbounded local realization without a numerical exponent bound.

For the smallest zero-demand box witness, every physical feasible length has an odd prime divisor whose geometric sum divides the original sum. Thus there is a prime witness no larger than any candidate length. At prime length only the order-one target equal to that prime, or a target whose order equals that prime, can contribute. The stated forbidden set is consequently complete. The nonzero-demand branch uses the least LCM length itself. A length ceiling must use the repaired physical prime witness in the L=1 branch, as the proof explicitly says.

The closure-fixed profiles are closed under coordinatewise minima. Closing the coordinatewise maximum gives their least fixed upper bound. The explicit source3 examples correctly show that raw maximum can fail and that the realizable arithmetic lattice contains M3. These are finite joins and meets, not a claim of a complete bounded lattice of finite profiles.

The OPN consequences use only necessary column balance and conditional overflow. The {3,13} example correctly separates local row realization and role-compatible column sums from globally compatible exponents. Untested exterior primes remain an additional obligation. Nothing proves OPN nonexistence or a new general OPN bound.

## Independent finite evidence

`check.py` imports no source implementation. It builds geometric sums by the direct integer recurrence and compares them with a separate order/rank implementation. The complete audit universe is stated in `receipt.json`: sources among the odd primes through43; two-target subsets of odd primes through29, excluding the source; both admissible source roles; direct lengths through240; every lower/upper box with endpoints0..2.

Both normal and optimized runs agree exactly. They give570,240 rank equivalences,580 source/role/target configurations,20,880 box comparisons,2,102 least-completion checks, and4,204 fresh-prime extension checks. Separate integer checks establish primality and orders for the M3 targets4561,368089,71; the actual profiles at lengths3,15,21,35,105; the weighted join at55; the Euler overflow at18; and the local/global hostile.

The author's verifier was independently rerun normally and with `-O`. Both outputs match one another and the supplied normal receipt byte for byte:190 configurations,41,040 bounded boxes,76,230 adjunction comparisons,1,198 closures and extensions, and all named controls pass. Explicit runtime checks remain active under `-O` in both implementations.

## Source and reproduction checks

The two cited repository prerequisite blob hashes were independently fetched at commit8880585ffd98072e93aa39d97476c7a1d4c9129d and match the proof's pins. The [journal landing page](https://journals.tubitak.gov.tr/math/vol42/iss4/3/) confirms the cited authors, title, year, volume, pages, and rank-of-apparition subject. The article is background only, not an imported proof dependency or a novelty claim.

Run `python verify.py` and `python -O verify.py` in the adjacent primeflow directory. Run `python check.py` and `python -O check.py` here. The audit receipt records the reviewed source SHA-256 values; `manifest.sha256` binds the final source and review files. Finite evidence supports the implementation and examples; the preceding mathematical proof review is the basis for the unrestricted local theorem.
