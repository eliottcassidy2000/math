# Independent review: absolute-error and dyadic recovery

Status: PASS for the scoped claims below. Review date: 12 November 2026.
This review covers new recovery mathematics and executable boundaries, not the
historical repository, Chan–Pak completeness, or the entire inherited cache.

## Mathematical review

1. **Every-sector lower bound is valid.** With r odd cycles, internal orbits
   contribute (n-r)/2 and cross-cycle orbits contribute at least r(r-1)/2.
   The centralizer is at most `(n/r)^r r! <= n^r`. Completing the square gives
   the displayed bound independent of r. At n=64 its exponent is 15/2.
   To justify monotonicity for every real n>=64, the auxiliary function
   n/8-1-log2(n) is positive at64 and has positive derivative there onward.
   Consequently the derivative bound used in the proof is valid. Since all
   omitted sectors are positive, the >181 error obstruction applies to every
   proper partial sum, not only fixed support cutoffs. It does not apply to an
   exact aggregate remainder, signed resummation, or modular completion.

2. **The interval is rigorous.** No unordered pair with a moved endpoint is
   fixed by an odd-order permutation: a fixed unordered pair would fix both
   endpoints or swap them, and swapping cannot occur at odd order. Its orbit
   therefore has length at least3. Counting these pairs proves the stated
   integer-rounded deficit. The odd-derangement coefficient recursion counts
   supported permutations, so multiplying its count by the largest possible
   fixed-tournament count produces an upper bound. Lower endpoint ceiling and
   upper endpoint floor preserve all integer candidates. Zero support and
   cutoffs beyond n are correctly handled.

3. **The dyadic arithmetic is exact.** A sector's power-of-two valuation is
   nonnegative; its odd denominator divides the odd part of n!, since the
   original centralizer divides n!. The common-denominator implementation can
   therefore multiply every retained term by an integer and use one final
   inverse. The added valuation at a multiplicity-group boundary is
   nonnegative, including the singleton terminal group, so a branch already
   vanishing modulo2^b cannot later reappear. The response backend aggregates
   this same arithmetic with no lossy right shifts. Its packed digit base
   exceeds n; every response coordinate is at most n, so no digit carry
   corrupts a neighboring coordinate.

4. **Integer uniqueness is strict and correct.** The condition is width < M,
   not width <= M. `lo+((res-lo)%M)` is the first congruent integer in the
   inclusive interval, including negative lower endpoints. The upper check
   detects inconsistent inputs. Width zero correctly needs zero residue bits.

5. **Precision-debt lemma is sound.** Assuming each division on an actual edge
   is exact, a congruence modulo2^B transports across an odd-unit-times-2^s
   edge to a congruence modulo2^(B+s). The backward maximum ensures this
   precision suffices at every consumer; addition preserves the congruence.
   Zero-bit terminal requests and negative shifts are included by the maximum
   with zero. This is a sufficient allocation, not an optimality theorem.

## Independent executable controls

`review_recovery.py` constructs Burnside values through order24 from an
independent ascending multiplicity-vector recursion, without reading the
Session270 cache or importing the production sector formula. It checks both
residue backends at precisions0,1,2,3,7,17,61,151; ten support cutoff choices per
order, including otherwise absent supports and cutoffs beyond n; and hostile
interval boundaries. All907 recorded checks pass. `RECOVERY_CONTROLS.json`
records exact source hashes for the reviewed implementations and proof.

This is a new-code control, not an enlarged numerical Chan–Pak census. The
all-order statements are established by the arguments above, not extrapolated
from order24.

## Scope and cost boundary

The algorithm is exact whenever run to completion with these certified inputs.
It is not an instantaneous enumeration theorem. Large-order approximations
leave a potentially enormous residue modulus; fixed support correction does
not remove the quadratic output-bit burden. The positive-sector obstruction
is not a lower bound for every conceivable future exact algorithm. This review
does not assert a speedup over the inherited optimized exact GMP engine, nor
complete classification of the infinite prime alphabet.

## Final compiled backend and census bridge review

The GMP backend applies the same nonnegative-shift response recurrence,
reducing weights modulo2^b after multiplication and addition and applying one
final odd-factorial inverse. Its integer response coordinates remain exact.
This review found only a CLI input-validation ordering issue: zero requested
bits returned before checking the advertised n range. The review repair moved the
range check before that early return. No arithmetic result for valid inputs
was affected. Independently compiled controls compare all n0..24 at five
precisions against the separate multiplicity enumerator, and check three
invalid inputs (128 checks).

The fixed-H census bridge correctly uses the inherited finite strong-core
formula to obtain polynomial growth in n, whereas every omitted Burnside
sector eventually exceeds any fixed polynomial. One draft sentence inferred
indistinguishability solely from interval width; a wide interval can exclude a
smaller comparison population if placed above it. This was repaired to the
actual lower-sum interval [L,U], where a(n)-L exceeds the removed fixed-fiber
mass. The corrected implication is valid and is not a prohibition against
all future proofs using enumeration.

The order250 GMP benchmark receipt is read as the author's bounded run, not an
independently replayed benchmark. It does not establish constant-time or
polynomial-time enumeration. The Python resource failure remains disclosed.
