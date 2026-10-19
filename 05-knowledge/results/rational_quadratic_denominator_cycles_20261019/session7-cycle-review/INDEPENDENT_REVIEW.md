# Independent review: rational quadratic denominator signatures

Date: 2026-10-19.

## Verdict

**ACCEPT the elementary reduction, the complete fixed-denominator algorithm, and the stated finite exclusion through reduced parameter denominator 1,000,000, relative to the explicitly cited Eliahou–Fares theorems.** I found no omitted sign branch, denominator case, period length, or cyclic-closure condition.

Independent implementations reproduce all 62 period-three cycles in the 101 searched denominators through `b=1000`, with no cycle of any other length at least three. The symbolic table is a proved *superset*, not an exact realizability classification. The author now makes this distinction explicitly.

**ACCEPT the added 43-adic obstruction.** Its Taylor argument proves more than the finite computation: no genuine period-six point exists anywhere in the specified three-cycle residue tube over `Q_43`. However, modulo `43³` there are still 21 period-six cycles; these are new, closer-to-center cycles reducing to period three modulo `43²`. Only lifts of the old period-six cycles grow to period 258. Preserve this distinction in summaries.

No general Flynn–Poonen–Schaefer conjecture, unbounded-denominator exclusion, priority claim, or proof of rationality from local compatibility is established.

## Material reviewed and separation of implementations

Reviewed the author's `denominator_signature_theorem.md`, particularly Sections 3–6 and 8, and the complete source of `rational_cycle_signature_probe.py`. The exact source hashes used at the final comparison appear in `comparison_receipt.json`. I also inspected the primary outputs through `b=400` and `b=1000` and the finite-lift script/output.

All new review files were written outside the repository, in this review directory. No Git writes were performed. The independent scripts import no author code. The projected-pair script imports only the reviewer's own elementary factorization function.

## Primary-source verification

Directly read the currently listed arXiv version, [Eliahou–Fares, arXiv:1903.08865v1](https://arxiv.org/pdf/1903.08865). The arXiv history lists only v1. The [publisher record](https://doi.org/10.4153/S0008439521000588) identifies the published article in *Canadian Mathematical Bulletin* 65 (2022), 598–611; the publisher's accessible abstract confirms the denominator and two-prime conclusions. The full theorem numbering used here is explicitly arXiv-v1 numbering.

- Proposition 2.1 supplies the common exact point denominator `b` and reduced parameter denominator `b²`.
- Theorem 2.12 supplies `4|b` for cycles of length at least three.
- Corollary 3.11 assumes reduced `den(c)=b²` and `4|b`; it bounds the total number of rational periodic points by `2^ω(b)+2`.
- Theorem 3.16 excludes periods greater than three when `ω(b)≤2`.

These are exactly the required hypotheses. Replacing the paper's `x²−c` by `x²+c` changes the parameter sign only. The implementation does not silently import an unconditional period-six theorem.

## Proof audit

### 1. Exact denominators and pair-sum support

The reduction requires `c=C/b²`, `b>0`, `gcd(C,b)=1`, and every numerator coprime to `b`. A reduced denominator which is not a square admits no rational periodic point. A cycle of length at least three forces `4|b`, so all numerators are odd and the three half-sums are integers.

Distinct cycle points cannot be negatives: they would have equal successors, whereas the map permutes its cycle. Thus every pair sum is nonzero. For any nonzero lag, multiplying

`b(a_(i+1)−a_(j+1))=(a_i−a_j)(a_i+a_j)`

around the exact cycle cancels nonzero factors and gives a product of sums equal to `b^n`. Every sum is an integer, so it has no prime divisor outside `supp(b)`. There is no hidden cancellation by rational denominators in this step.

### 2. Every individual bad-prime word has the full global period

At an odd bad prime, unit numerators lie in two nonzero opposite classes modulo that prime. Equal symbols make their sum a unit; unequal symbols make their difference a unit. Equal-symbol iteration lowers the difference valuation by `v_p(b)>0`.

At two, the two symbols are the odd classes modulo four. Equal symbols give sum valuation one, and unequal symbols give difference valuation one. Equal-symbol iteration lowers the difference valuation by `v_2(b)−1>0`.

If two distinct cyclic shifts had identical infinite symbols at even one bad prime, this positive drop would continue forever, contradicting integrality. Therefore each individual local word is primitive of length exactly `n`. It is not merely necessary that the least common multiple of several shorter local word periods equal `n`.

Stopping at the first differing symbol proves

- odd `p`: `v_p(a_i−a_j)=e_p L(i,j)`;
- `p=2`: `v_2(a_i−a_j)=1+(e_2−1)L(i,j)`.

One further use of the difference identity proves the author's pair-sum formulas. Dividing sums by two removes exactly the additive one in the 2-adic formula, yielding bases `p^e_p` for odd primes and `2^(e_2−1)` for two. Both the offset and the smaller 2-adic base are necessary.

### 3. Three-edge signature bound

The four initial binary types `000`, `001`, `010`, `011`, up to complement, give precisely the claimed containing families. For the `011` case, the first run of ones produces adjacent and two-step disagreement times `k` and `k−1`. For `001`, the next bit agrees with exactly one of the two previous unlike bits, so exactly one nonzero waiting time is one.

Two distinct rotations have the same Hamming weight. Their disagreement positions cannot consist of a single position. Thus, after a disagreement, another occurs in at most `n−1` steps. This proves the upper exponent bound, including wraparound cases, and the containing set has `4N−6` elements.

The table does contain nonrealizable entries at a fixed length. For example, at exact length ten there are 31 realized triples, while the stated superset has 34. The union over all lengths three through ten has 32. These extra entries are harmless because the procedure only claims completeness followed by exact acceptance tests.

The cycle-product valuation identities do not supply a further obstruction: for any primitive binary word and any lag, the positive disagreement gaps partition the cyclic index set, so their sum is `n`. This validates the paper's stated limitation of circulation-only arguments.

### 4. Signs, reconstruction, and closure

With `u=(a_0+a_1)/2`, `v=(a_0+a_2)/2`, `w=(a_1+a_2)/2`, inversion is exactly

`(a_0,a_1,a_2)=(u+v−w,u+w−v,v+w−u)`.

Subtracting the first two recurrence equations gives exactly

`u(v−w)=(b/2)(u−v)`.

The primary implementation tries both signs of `u` and both signs of `v`, and retains both possible signs of the prescribed magnitude of `w`. Zero is correctly absent, and collisions among the initial three vertices are rejected. Primitive numerator checks establish the exact denominator. The additional `gcd(C,b)=1` check is redundant after `C=b a_1−a_0²` and `gcd(a_0,b)=1`, but safe.

The loop accepts only a first return to the starting point, after at least three distinct points. It rejects a tail entering a different cycle. The `range(N+1)` placement checks for return before the `step==N` break, so a genuine cycle of length exactly `N` is not lost. Nonintegral or nonprimitive subsequent numerators are rejected. No reduction modulo a good prime is used to infer a global period.

There is no free parameter-numerator search bound: every actual cycle determines one of the finitely many signed half-sum candidates, from which `C` is reconstructed. Exact recurrence checking provides the converse.

### 5. All denominator and period cases in the finite conclusion

Through `b=400`, only 28 multiples of four have at least three distinct prime factors. They all have three prime factors, since the least eligible four-prime denominator is `420`. The external theorem bounds their periods by ten, and the implementation checks all lengths three through ten, not merely length six.

Through `b=1000`, there are 96 three-prime cases and five four-prime cases (`420,660,780,840,924`). The latter have period bound eighteen. Five distinct prime factors would require at least `4·3·5·7·11=4620`. Cases with at most two primes are covered by the cited theorem; nonmultiples of four cannot have period at least three; nonsquare reduced parameter denominators admit no periodic points. Therefore the finite exclusion really does cover every rational parameter of reduced denominator at most 1,000,000, without bounding the numerator.

## Independent exact controls

### A. Exhaustive primitive binary words

`check_words.py` independently enumerates every primitive binary cyclic word of lengths three through eighteen, up to complement. It computes each exponent as the first subsequent disagreement, without using the author's signature formula. No computed triple is omitted by the containing table. Detailed counts and surplus triples are in `word_results.json`.

### B. Broader smooth-pair enumeration through b=400

`independent_search.py` uses no binary-word or triple-signature table. It enumerates every pair of signed half-sums whose bad-prime exponents are independently between zero and `N−1`, derives `w` from the recurrence, and only then tests its bounded smooth support and the exact orbit.

It tests 2,678 primitive initial triangles versus 1,653 in the author's three-edge search. Both methods recover exactly the same 28 cycles, all of period three, for all 28 denominators. The independent control cases `b=4,12,24,60,84,120` also agree.

### C. Independently generated projected-word pairs through b=1000

`word_projected_search.py` enumerates primitive words directly for every length `3..N` and retains only the pair `(q(0,1),q(0,2))`. It does not use the closed-form triple table. It enumerates the resulting two-half-sum possibilities, calculates `w`, and checks only bounded smooth support before exact recurrence closure.

At `N=18`, the independently generated pair table has 34 entries. This requires `34^4=1,336,336` unsigned pairs per four-prime denominator, rather than `66^4=18,974,736` unsigned triangles. All five four-prime cases are checked, not skipped.

The complete 101-denominator run produces 8,906 primitive initial triangles and exactly the same 62 period-three cycles as the primary output. No longer cycle is found. Both this run and the broader `b≤400` run were repeated under `python -O`; mathematical outputs agree. Explicit comparison code raises exceptions on mismatches rather than relying on assertions. See `comparison_receipt.json`.

### D. Complete finite graphs and Taylor coefficients at 43

`check_43_lift.py` independently decomposes the complete functional graphs modulo `43²` and `43³`; it does not start from the primary list. It finds 126 period-six points modulo `43²`, and all their 5,418 lifts lie in 21 cycles of exact length 258 modulo `43³`.

It also detects 21 different period-six cycles modulo `43³`. They are not lifts of the old period-six points. This is the essential finite-level qualifier.

Exact polynomial arithmetic at the three centers gives `λ=35/8`, `λ²−1=1161/64`, and quadratic coefficients of `T³` equal to `1107/64`, `−297/64`, and `−81/64` at centers `−7`, `5`, and `−1`. All higher coefficients are 43-adically integral. The composition's quadratic coefficient is `Aλ(λ+1)`, divisible by 43. Thus, for displacement valuation `s≥1`, the linear term has valuation `s+1`, the quadratic term at least `2s+1`, and the cubic remainder at least `3s`. The linear term is uniquely minimal, proving the uniform local obstruction.

## Reproduction

Run from the workspace root:

```sh
python session7-cycle-review/check_words.py
python session7-cycle-review/independent_search.py
python -O session7-cycle-review/independent_search.py --output session7-cycle-review/independent_results_optimized.json
python session7-cycle-review/word_projected_search.py
python -O session7-cycle-review/word_projected_search.py --output session7-cycle-review/projected_results_optimized.json
python session7-cycle-review/check_43_lift.py
python session7-cycle-review/compare_receipts.py
```

The files ending in `_optimized.json` include runtimes, so full-file byte identity is neither expected nor claimed. The comparator checks denominator lists and every exact parameter, numerator cycle, and period.

## Addendum: two-edge compression and arbitrary-order resonance

These refinements were reviewed after the initial census and do not alter its scope.

### Universal two-edge candidate theorem

Projecting the proved containing table onto its first two coordinates yields

`P_N={(0,k):0≤k≤N−1} ∪ {(k,k−1):1≤k≤N−1}`,

with exactly `2N−1` elements. Every pair satisfies `q_01≤q_02+1`. Because `h=b/2` is exactly the product of the half-sum prime bases, any signed `u,v` assembled from these pairs automatically satisfies `u | h v`. Thus `w=v−h+h v/u` is integral without an extra divisibility search.

There are at most `4(2N−1)^ω(b)` signed candidates before reconstructing vertices and testing recurrence closure. This is a candidate-count bound, not a bit-complexity or total-runtime claim. The bounded smooth-support test on `w` is a safe acceleration; logically the full exact orbit test is sufficient to reject an inadmissible third edge. The argument remains complete because every actual cycle's first two half-sums occur in the projected table.

`universal_pair_search.py` independently replays this larger universal pair table rather than the smaller word-realizable projection. Its output `universal_pair_results.json` again matches all 101 denominators and all 62 period-three cycles exactly. It imports the reviewer's own reconstruction/closure code, not the author's implementation.

### Stronger simple-resonance theorem: exclusion of all nontrivial periodic points

Let `p` be odd, `g∈Z_p[t]`, `g(α)=α`, `α∈Z_p`, `λ=g'(α)`, and let `r≥1` be an integer. Assume

`λ−1 ∈ Z_p^×` and `v_p(λ^r−1)=1`.

Then for every nonzero `t∈p Z_p` and every integer `N≥1`,

`v_p((g^r)^N(α+t)−α−t)=v_p(t)+1+v_p(N)`.

Consequently `α` is the only periodic point of `g` in the ball `α+p Z_p`.

**Proof.** Write `g(α+t)=α+λt+A t²+O(t³)`. The quadratic coefficient of `F=g^r` is

`A λ^(r−1)(1+λ+⋯+λ^(r−1)) = A λ^(r−1)(λ^r−1)/(λ−1)`.

It is divisible by `p`, while all higher coefficients are integral. Put `s=v_p(t)≥1`, `t=p^s y`, where `y` is a unit. The scaled map

`H_s(y)=(F(α+p^s y)−α)/p^s`

has the polynomial form

`H_s(y)=y+p u y+p²V(y)`, where `u=(λ^r−1)/p` is a unit and `V∈Z_p[y]`.

Use the following elementary iteration identity. If `H(y)=y+p^a B(y)` with `a≥1`, then for odd `p`,

`H^p(y) ≡ y+p^(a+1)B(y) (mod p^(a+2))`

as polynomials over `Z_p`. To check it, expand successive compositions to second order: the additional quadratic-in-`p^a` term has factor `binom(p,2)`, divisible by `p`, and all remaining terms have valuation at least `3a≥a+2`. Induction gives

`H_s^(p^k)(y)=y+p^(k+1)B_k(y)`, with `B_k(y)≡u y (mod p)`.

For a positive integer `m` prime to `p`, another first-order composition gives

`(H_s^(p^k))^m(y)−y ≡ m p^(k+1)B_k(y) (mod p^(k+2))`.

At a unit `y`, the coefficient on the right is a unit times `p^(k+1)`. Taking `N=p^k m` and undoing the scaling proves the valuation formula. It is finite for every `N`, so a nonzero displacement cannot be periodic for `F`, and therefore cannot be periodic for `g`. QED.

For the original map `T(a)=(a²−29)/4`, take `g=T³`, `p=43`, and `r=2`. Here `λ=35/8`, `λ−1=27/8` is a 43-adic unit, and `λ²−1=1161/64=43·27/64`. The union of the three balls centered at `−7,5,−1` contains no `Q_43`-periodic points other than those three centers. More precisely, any nonzero displacement `t` from a center satisfies

`v_43(T^(6N)(α+t)−α−t)=v_43(t)+1+v_43(N)`.

This is an elementary tangent-to-identity iteration argument; no novelty claim is made for that general mechanism.

**Boundary caution.** The `λ−1`-unit condition cannot simply be dropped. For `p=3`, `g(t)=4t−t²` fixes `α=0`, has `v_3(λ−1)=1`, and nevertheless has the second fixed point `t=3` in the same ball. The odd-prime restriction is also retained throughout the iteration argument; no 2-adic analogue is asserted.

`check_general_resonance.py` performs 7,875 exact finite-modulus controls for primes `3,5,7,11,13`, resonance orders `2,4,3,10,12`, 25 quadratic/cubic coefficient choices per prime, three displacement depths, and several iteration counts including prime squares. All satisfy the stated valuation identity. The controls and boundary hostile appear in `general_resonance_controls.json`; the proof above supplies the unbounded statement.

### All-depth finite-shell classification at 43

For modulus `43^k`, the three-center tube consists of the central three-cycle and, for every `1≤s≤k−1`, the shell of points with displacement valuation exactly `s`. Every point on that shell has exact `T`-period

`6·43^(k−s−1)`.

Indeed, reduction modulo 43 forces any return time to be divisible by three. An odd number of `g=T³` iterates changes the displacement by a quantity of valuation exactly `s`, because its linear multiplier reduces to `−1` and all nonlinear terms have valuation at least `2s`. Such an iterate cannot return modulo `43^k`. For even `g`-iterate counts, the exact sixth-return valuation formula determines the least required power of 43.

The shell has `3·42·43^(k−s−1)` points across the three disks. Dividing by the exact period gives precisely 21 cycles on every shell. The central three-cycle together with these shells contains

`3 + sum_(s=1)^(k−1) 3·42·43^(k−s−1) = 3·43^(k−1)`

points, which is the entire tube. This explains simultaneously why old period-six cycles grow and why new period-six cycles continue to appear closer to the three centers at every finite depth `k≥2`.

## Final bounded addendum: complexity and critical-cycle controls

**ACCEPT the factorization-given fixed-parameter enumeration theorem in Section 6 and the critical-cycle controls in Section 9A.** The new universal-pair implementation also passes replay. One minor wording correction was requested: additional shell cycles occur at every finite depth `k≥2`, while depth one has only the central three-cycle. This does not change any theorem or computation.

### Factorization-given bit complexity

The supplied factorization is an explicit precondition, so the theorem makes no claim about deterministic polynomial-time factoring. With `s=ω(b)` and `N=2^s+2`, the number of signed seeds is at most `4(2N−1)^s`. Each prime-base magnitude has at most `(N−1)log_2(b/2)` bits. The stated bounds on `w`, the reconstructed vertices, and `C` therefore imply `O(N log b)` bits uniformly.

The escape guard is mathematically safe. For `A=|a|` and `S=sqrt(|C|)`, the condition `A>b+floor(S)+1` implies `A>b+S`. Hence

`a²+C ≥ A²−S²=(A−S)(A+S)>bA`,

so the next iterate is positive and strictly exceeds `A`, and the same inequality continues forever. Such a point cannot be periodic. All nonrejected orbit values are bounded by `b+floor(S)+1`, and computing one next numerator only doubles the bit bound by a constant factor. The period bound permits at most `N` iterations per seed. Standard deterministic arithmetic gives a computable function of `s` times a polynomial in `log b`, whose exponent is independent of `s`. The precomputed bounded-smooth set has `N^s` elements, also absorbed into that parameter-only function.

This proves fixed-parameter tractability of *enumeration at an input denominator*, allowing all parameter numerators. It does not uniformly exclude longer cycles for every denominator. Restricting the enumerated periods to at least three is essential; fixed points and two-cycles vary through infinitely many parameter numerators at fixed denominator.

The new `rational_cycle_pair_probe.py` correctly implements the universal pair table, sign choices, exact reconstruction, denominator tests, and real-escape guard. Acceptance does not depend on Python assertions. The supplied-factorization CLI expressly assumes prime factors rather than silently claiming primality certification. Its bounded-denominator mode uses elementary trial division. I replayed all 101 denominator cases and recovered the same 62 cycles; see `author_pair_full_replay.json`.

For `b=4·3^20·5^20=1330102692031860351562500`, I independently ran my own universal-pair reconstruction and exact recurrence code **without** the new escape guard. It gives 6,859 unsigned pairs, 16 primitive triangles, and no cycle of period at least three, agreeing with the author's guarded result. The author's script reports 27,436 signed seeds and six escape rejections. Fresh normal and optimized executions produce byte-identical JSON files, `author_large_exponent_replay.json` and `author_large_exponent_replay_optimized.json`.

### Critical-cycle contraction theorem

For `c∈Z_p`, a nonintegral point escapes because its negative valuation doubles at each iterate. Every integral periodic point reduces to a finite-field cycle. If the finite-field map has a unique cycle of length `m`, every periodic point must lie in its residue disks.

The map `f_c^m` preserves each such disk. The zero derivative product implies its derivative vanishes modulo `p` throughout that disk. For two points in the same disk, the polynomial difference quotient equals the derivative at one point modulo `p`; thus the quotient is divisible by `p`, giving Lipschitz constant at most `1/p`. Completeness of the disk and strict contraction give one fixed point and exclude longer periodic orbits. The original map permutes these fixed points with exact period `m`, since their reductions are distinct. This proves the stated unique affine periodic cycle over `Q_p`.

The residue graphs are correct: `c≡1 mod5` gives the sole cycle `0→1→2→0`; `c≡2 mod3` gives the sole cycle `0→2→0`. Both include the critical point zero and have multiplier zero. The examples `−29/16` and `−7/4` belong to their respective classes. The qualification about possible absence of rational points is necessary and is present: `c=1` has the local three-cycle but no real periodic point, since `x²+1>x` for every real `x`.

My `check_final_addenda.py` independently decomposes every complete finite graph for every compatible parameter residue modulo `p^k`, for `1≤k≤4`. All 40 parameter cases for `p=3` and all 156 for `p=5` have exactly one cycle of the claimed length. Its output `final_addendum_controls.json` also contains the independent large-exponent search. These finite checks support the separately proved infinite-family contraction statement; they are not substituted for that proof.

Final reviewed author-note SHA256: `a9ba51d9487c2c31bdc0fabe58b6eb9175938fefcbd22884ce5f14b59265aadf`. The requested depth-one wording correction is confirmed in both the author note and this audit. Final verdict remains **ACCEPT**.
