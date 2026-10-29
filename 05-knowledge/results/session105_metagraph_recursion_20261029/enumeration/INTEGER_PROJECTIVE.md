# Integer realization of the fixed-degree projective response recurrence

**Status: PROVED arithmetic realization; FINITE-EXACT verification; measured implementation experiment.** No new Burnside formula, new projective quotient, universal speedup, bounded-state theorem, or classification of tournament representatives is claimed.

## Inheritance and actual baseline

This lane recovers and combines, rather than replaces:

1. `05-knowledge/results/tournament_future_boundary_20261020/README.md`, finite-horizon gcd-response pruning and insertion.
2. Session 98, recovered locally as `session98_growth/recurrence/PROJECTIVE_RESPONSE.md` (published packet `05-knowledge/results/session98_tournament_growth_20261029/recurrence/PROJECTIVE_RESPONSE.md`), fixed-degree projective normalization `s_j=h_j-jh_1` and its minimality for individual exact-degree continuations. The exact inherited checker is copied unchanged as `baseline98.py`.
3. Session 87, `session87-four-vertices/counting/fast_counts.py`, specifically **degree_bucket_dp**, not its slower direct-partition routine: packed responses, shared degree/multiplicity metadata, odd-factorial integer arithmetic, and terminal unit-tail fold. Copied unchanged as `baseline87.py`. This baseline returns every coefficient through N, whereas our projective engine returns only coefficient N.
4. Session 88, `archive/poset_future_state_theorem.md`, §2: full gcd/divisor response faithfulness, with an explicit boundary between a future-operation state and aggregate algorithms. Projective fixed-degree merging does not contradict its all-future faithfulness result.
5. THM-514, `the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`, and THM-516, `burnside-core-kernel-phi-reframe.md`: analytic/formal tail elimination and gcd energy; these are inherited mechanisms, not new discoveries.

The source hashes are retained in the verification receipt. Source copies deliberately preserve their original function and module documentation. No claim is made to be faster than unbenchmarked C/GMP implementations.

**Concept board.** Fixed-degree observability; exact integer lattices under a projective gauge; carry-safe packed response deficits; cached transfer metadata; soluble tails. Hostile controls are nonexact right shifts, composite odd lengths, and the faster all-coefficients Session87 baseline.

## 1. Full-factorial integer gauge

Fix N. A history mu contains already processed odd lengths, total mass d=N-r, t cycles, energy e and centralizer z. All future lengths are strictly smaller than processed lengths. The inherited normalized response is

`s_j = sum_(a in mu) gcd(a,j) - jt`, with `s_1=0`.

Session98 accumulates `W = sum_mu 2^(e+rt)/z` in a state `(r,s)`. Define instead

`V = (N!/r!) W`.

**Theorem.** Every individual history contribution to V is a nonnegative integer. Inserting k cycles of the current odd length p, with r'=r-kp, updates it by

`V' = V * [r!/(r'! p^k k!)] * 2^E`,

`E = k(p-1)/2 + p*k(k-1)/2 + k*s_p + r'*k`.

The bracket is an integer. If E<0, division by `2^(-E)` is exact **after multiplication by that bracket**, including for a merged state. The successor response is `s'_j=s_j+k(gcd(p,j)-j)`, restricted to j<=r' and j<p.

**Proof.** `N!/(r! z)` counts permutations with prescribed processed cycles on N labels and r unused labels, so is integral. The exponent e+rt is nonnegative. The bracket counts selections of k disjoint p-cycles on r available labels. Algebraically divide `(N!/r'!) W'` by `(N!/r!) W` using the inherited Session98 insertion identity. The power E can be negative, but each old history maps to an individual new integer contribution by the preceding combinatorial argument. All histories in one state share E and the bracket, so their sum is divisible as well. No cancellation, approximation, or post-hoc rounding is involved. QED.

Unit cycles can be filled without materializing their stage: add `V*2^(r(r-1)/2)` to the final numerator, then divide the total by N!. This follows because the unit-tail contribution is `2^(rt+r(r-1)/2)/r!` before normalization. For N=0 the initial empty history directly gives 1.

This is the new arithmetic bridge: Session98's sometimes negative rational transitions admit an exact integer implementation **without abandoning the projective merge**. Clearing denominators itself was already present in Sessions22/87; the new point is compatibility with this coarser projective state and its negative shifts.

## 2. Smaller inherited odd-local scale also survives the gauge

Let `O(d)=odd(d!)`, `v(d)=v_2(d!)=d-popcount(d)`, and define

`U=O(d) W`.

For histories made entirely of processed p>=3, U is integral: `O(d)/odd(z)` is integral, while `e>=sum_p v(m_p!)=v_2(z)`. Indeed each p>=3 contributes at least m_p to e and `v(m_p!)<=m_p`. The extra rt is nonnegative. Thus the same state merging can use smaller weights.

For d'=d+kp, its transfer is

`U' = U * [O(d')/(O(d) p^k O(k))] * 2^(E-v(k))`.

The bracket is integral: it is the odd part of the integer `d'!/(d! p^k k!)`. Negative shifts are again exact by individual-history integrality. To finish the unit tail, multiply a state by

`[O(N)/(O(d) O(r))] * 2^(r(r-1)/2-v(r))`,

sum, and divide by O(N). The bracket is the odd part of binomial(N,r). This merges Session87's smaller integer scale with Session98's projective state; neither ingredient alone is claimed as new.

## 3. Carry-safe packing and reusable transitions

Put `D_j=-s_j=sum_(a in mu)(j-gcd(a,j))`. In descending elimination, all processed a exceed every active j. Hence `0<=D_j<=jt<=sum a=d<=N`. This uniform bound also holds immediately after insertion for each retained j<p. Packing each D_j in `max(1,bit_length(N))` bits is therefore carry-safe. Pack only j=3,5,..., since D_1=0. The top digit gives `-s_p`; addition of `k*(j-gcd(p,j))` gives the next response, and a mask expires inaccessible coordinates. This is a proof for the actual descending order; it is not an unrestricted bound for arbitrary operation orders.

`count_packed` buckets by r and shares coefficients, offsets and masks across all states with that r. The tuple reference implementation `count` has an optional `(p,r)` transition-template cache, allowing the cache alone to be measured. Tail=3 explicitly sums all three-cycle multiplicities and then fills units, using the inherited soluble-tail formula; it is not a new quotient. Both tails retain every composite odd length.

## 4. Quotient and structural boundary

No strictly coarser **termwise exact-r projective response** is obtained here. Session98's proof is already decisive: compare the r-unit continuation with one j-cycle plus r-j units. Any history-independent ratio forces `h'_j-h_j=j(h'_1-h_1)` for every available j<=r. That obstruction is tied to individual continuations; it does not exclude identities among sums, signed algorithms, or other representations.

The map from histories to `(r,s,U)` preserves the weighted coefficient under all remaining insertions. It loses the individual cycle multiset, its centralizer, and all information about individual tournament automorphism groups, extension distributions, or metagraph edges. Recovering these requires an additional structural carrier, not a faster Burnside evaluation. All automorphism weights remain in the aggregate arithmetic; no regularity, prime-length, rigidity, or pair-free restriction is introduced. This result is not a finite structural prime atlas.

## 5. Reproduction and limitations

Run from this directory:

```
python3 verify.py > verification.json
python3 -O verify.py > verification.optimized.json
cmp verification.json verification.optimized.json
python3 benchmark.py > benchmarks.json
```

Verification uses exceptions, never removable assertions. It compares direct odd partitions through 30; full-factorial and odd-local packed scales, cached/uncached tuple implementations and both tail folds through 30; then inherited Fraction98 and all-coefficients integer87 at 40,60,100,150,200. Exact projective stage counts agree with Session98 after omitting only the analytically folded unit stage. Invalid inputs and a deliberately nondivisible right shift are rejected. The packed carry bound is checked on all odd partitions through 30.

Timing is a measured implementation property, not a theorem. In particular the faster inherited Session87 implementation returns a strictly richer output, and must remain in the benchmark. See `benchmarks.json` for repeats, machine/Python details and all tested variants. No historical priority claim is made beyond the recovered repository comparison.

### Measured comparison

Median seconds from three runs; N=200 uses the separate replay after this lane's main checker had finished. Other agents may still share the machine. These are wall-clock observations, not isolated CPU benchmarks.

| N | Session87 all coefficients | Fraction98 single target | Packed odd-local single target |
|---:|---:|---:|---:|
| 40 | 0.0015 | 0.0084 | 0.0008 |
| 60 | 0.0048 | 0.0382 | 0.0031 |
| 100 | 0.0324 | 0.3454 | 0.0266 |
| 150 | 0.2311 | 3.1490 | 0.2635 |
| 200 | 1.3817 | 19.9700 | 1.9507 |

At N=200 the new packed realization is about **10.2× faster than the inherited Fraction projective engine**, but **41% slower than Session87**, which moreover produces the whole prefix. At N=40/60/100 it is faster than that richer-output baseline in these measurements. **Retain Session87 for production prefix queries; no uniform performance improvement is established.**

The N=200 tuple-cache ablation falls from 8.20 s uncached to 5.98 s cached; the three-tail variant takes 4.35 s. Tail folding and metadata sharing help this implementation, but are inherited ideas and do not improve the state quotient. At N=100 three-tail folding was slightly slower than cached one-tail (0.0633 vs 0.0620 s); the extra elimination is not an unconditional win. The packed engine uses 95,709 peak projective states at N=200 versus 96,836 raw-response states in Session87: fewer states alone do not guarantee faster arithmetic.
