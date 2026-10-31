# Independent audit: constructive interval conductor

2026-10-31. Scope: rational inputs, positive integer speeds, strict bad sets / closed safe sets. This audits `../lrc/CONSTRUCTIVE_CONDUCTOR.md`, not unrestricted LRC or a claimed new lower bound.

## Proof findings

The proposed extractor is correct under its stated measure-budget hypotheses. No counterexample was found in the proof audit.

1. **Equality repair is genuinely constructive.** If the left endpoint a is bad in component (c,e), then e>a. The strict individual bound μ(B_i)<W forces e<b: otherwise the component covers [a,b) and has measure W. If e is unsafe, a different bad component contains e in its interior. Its intersection with (c,e) and I has positive rational length. Adding its indicator to 1−Σ1_B never exceeds the safe indicator, because that overlap has bad multiplicity at least two. The corrected integral is strictly positive.
2. **Isolated solutions are not silently discarded.** At equality, if only measure-zero safe times exist, the corrected-positive-integral branch is impossible. The first component endpoint must already be safe. This follows from the same pointwise inequality, without assuming a positive safe margin.
3. **Terminal endpoint test is sound.** All bad event points have denominators at most QV, hence distinct events are separated by at least 1/(Q²V²). A smaller interval contains at most one interior event. Positive safe measure then yields a safe open cell adjacent to an interval endpoint. Closedness of the safe set makes that endpoint safe. Coincident events, an event on an endpoint, and overlapping components do not invalidate this argument.
4. **Bit complexity is polynomial.** log Q is at most the sum of denominator bit lengths; log V is input size. The needed dyadic depth is polynomial in input length, including when interval endpoints are large or very close. Exact floor handles arbitrarily many traversed periods in one operation. A sum of k rational measures has denominator bit length bounded by a sum of polynomial input/intermediate bit lengths; it does not require enumeration of teeth. Ordinary integer arithmetic and gcd/lcm are polynomial-bit computations. The stated arithmetic count excludes preprocessing details but its overall polynomial-bit conclusion is sound.
5. **Phases and negative times work.** The antiderivative with floor, rather than truncation toward zero, is valid on the full real line. Arbitrary rational phase magnitudes require only polynomial-sized integer parts.
6. **Optional extension.** For positive reduced rational speeds v_i=p_i/q_i, replace V in the separation bound by P=max p_i. Event denominators divide Qp_i. The same proof and bit complexity apply. Using max rational speed itself would be incorrect for the denominator bound.

## Adjacent categorical statements

For a fixed rational ambient interval (or circle), finite unions of rational intervals with all endpoint flags form a Boolean algebra. Intersection with ambient unit is symmetric monoidal closed, with A⇒B=(ambient\\A)∪B. This is a thin cartesian closed category, not thereby a realizability topos. Closed-only subsets fail closure under this implication. Nonemptiness does not preserve tensor: two disjoint nonempty sets have empty intersection.

For m≥2 and 0<r≤1/2, the safe relation ||mt+φ||≥r on two circles has no continuous section φ↦t(φ). A section would make g(φ)=mt(φ)+φ land in the safe arc, which avoids0 and so has degree0. But deg(g)=m deg(t)+1, impossible. Each fibre is nonempty; this is a topological obstruction to continuous selection, not an obstruction to discontinuous rational computation. For m=1 the obstruction correctly disappears.

## Independent finite checks

See `check_independent.py` and its normal/optimized receipts. Tests use independently enumerated rational event cells for exact measures and safe geometry; this is deliberately finite audit code, not the production polynomial algorithm. The extractor below uses a separately implemented periodic primitive. All checks raise explicit exceptions rather than removable assertions.

### Recorded result

Both normal Python and `python -O` passed the deterministic 1,203-input corpus: 804 admissible inputs, 399 rejected measure-budget controls, 4 positive-overlap repairs, 490 halvings, 579 inputs with negative left endpoints, and 4,218 direct measure-oracle comparisons. This is finite executable evidence, not a replacement for the arbitrary-input proof. The corpus also exercises reduced rational speeds using the numerator separation bound.

A useful complexity detail is that no lower bound on the positive density integral is needed for termination. Exact sign tests preserve positivity; the final stopping criterion depends only on geometric event separation. In particular, tiny positive overlap lengths do not trigger a hidden precision search.

### Producer implementation cross-check

`witness.py` was inspected separately against the proof. A further 500 deterministic integer-speed cases were classified using the independent full-event geometry and passed through the actual producer: 325 accepted with verified common witnesses, 175 correctly rejected. Normal and optimized runs both passed; see `check_producer.py`, `producer_normal.txt`, and `producer_optimized.txt`. This implementation cross-check deliberately imports the producer, whereas `check_independent.py` does not. The latter includes the rational-speed extension; the producer's public contract remains integer speeds.

## Full-family selector follow-on audit

Inspected `FAMILY_SELECTOR.md` and `family_selector.py`, including the completed dispatcher. Their use of the Session122 universal eight-core/five-tail theorem is explicit and essential; this audit does not substitute finite testing for that inherited theorem.

* Promotion/padding is sound. Before each removal of a bounded tail, at least one bounded core slot remains; an already-present tail can be replaced by the smallest missing bounded speed. This strengthens the original row, and the enlarged core with remaining tails is still an instance of the eight-core theorem. There is no assumption that arbitrary tail speeds are distinct.
* A nonempty point-only rational prefix cannot admit every completion when a tail slot remains: a speed divisible by all point denominators, repeated in the remaining slots, kills every safe point. Thus a positive-width prefix component exists. At zero remaining slots, this argument is unavailable and the implementation correctly permits point-only covers.
* With assigned integer speeds bounded by M, endpoint differences have the form integer/(14uv), not merely a generic denominator-product bound of 1/(196M²). Hence every positive component width is at least 1/(14M²). The upper bound ceil(c_k/W)≤14c_kM² is valid because 14c_k is one of the integers 2,4,6,8,20. At most five expanded finite heads then lie in an input-independent finite universe. This gives a polynomial-bit theorem with potentially huge fixed constants, not a practical uniform bound as the core universe grows.
* `attach` enumerates exactly those safe strips whose closed endpoints meet the current closed interval: ceil(v lo+r−1)≤n≤floor(v hi−r). Inclusive inequalities preserve isolated terminal points. The sorted-tail cutoff makes every remaining speed satisfy the conductor cutoff. The k=1 formula correctly includes both safe boundaries and the equality-width case.
* The returned value is checked on all original core and tail speeds. No quotient clock or phase is substituted.

Final reading of `SAME_CLOCK_CATEGORIES.md` passes, including the no-subobject-classifier argument for the nontrivial posetal category, the failure of a largest closed disjoint object, both exact disjoint phase witness sets, and the local-section/global-degree distinction. The constructive certificate-compilation corollary has the necessary fixed finite, well-founded, witness-preserving hypotheses; it does not establish an unrestricted family theorem.

### Full-family executable evidence

The separate `check_family_independent.py` imports only the public family selector, not its comb intersection routine. It constructs complete safe geometry independently by sorting every rational bad-boundary event on [0,1], checking each event and each intervening midpoint with integer modular arithmetic, and retaining safe cells/points.

Normal and optimized runs pass: all 1,287 eight-cores against complete geometry for tails (14,15,16,17,18); 1,287 repeated bounded-tail rows; all 1,093 larger-core choices; eight mixed huge-input controls through 4,096 bits; and seven invalid-input rejections. The full-geometry rows exercise 841 conductor returns, 80 finite leaves, and 366 one-tail returns. The completed bounded pool's safe set consists of six isolated points, and the algorithm preserves a witness. An additional exact one-tail boundary control uses core {1,...,12}, tail 24 and maximal width 1/168, so vW=1/7 exactly.

Status: **PASS on proof and implementation audit**, conditional on the explicitly inherited Session122 universal family theorem. Finite observations do not establish a universal small cutoff, and no unrestricted LRC result is asserted.
