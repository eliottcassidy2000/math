# The asymptotic unit-distance problem: a positive count with a cancelling spectral kernel

**death-star-2026-07-22-S108** (HYP-8945). Owner: work the unit distance problem; creative new routes; mine
unrelated threads. **The asymptotic Erdős unit-distance problem is open (upper `O(n^{4/3})`, lower
`n^{1+c/loglog n}`, a 40-year gap); this is exploration, not a solution.** The repo has a real unit-distance
cluster, but it is *small-`n` exact* (THM-431 `u(21)=57`, THM-432 Moser-lattice optimum, citing
Alexeev–Mixon–Parshall) and *lattice density* (THM-412 density-quantization of the unit-distance layer, THM-421
the `3N` floor) — plus the Moser/chromatic (Hadwiger–Nelson) side. It has **nothing** on the asymptotic gap or its
tools (Szemerédi–Trotter, Guth–Katz, sums-of-two-squares: 0 files each). This note opens that problem in the repo:
verified anchors, and — the point of the session — a precise **cross-thread placement** that says *why* the repo's
signature move (the positive-definite manoeuvre that unlocked GMC) does **not** transfer here, and which of the
repo's threads unit distances actually belongs with.

## 1. The problem and state of the art

`u(n)` = max number of unit distances among `n` points in `ℝ²`. Erdős (1946): `u(n) = n^{1+o(1)}`, conjecturally
`u(n) ≤ n^{1+c/loglog n}`, matched by the `√n×√n` grid. Proven: lower `n^{1+c/loglog n}` (grid), upper
`O(n^{4/3})` (Spencer–Szemerédi–Trotter 1984, unchanged in 40 years). The gap between `n^{4/3}` and the believed
`n^{1+o(1)}` is the whole problem.

## 2. Verified anchors (`unit_distance_asymptotic_probe_deathstar_S108.py`)

- **Lower bound is number-theoretic.** In the `m×m` grid (`n=m²`), the count at squared-distance `k` is
  `≈ ½·r₂(k)·n`, where `r₂(k)` = #`{(a,b)∈ℤ²: a²+b²=k}`. Verified: the maximizing `k*` and the count `U/n` both
  grow — `U/n` ran `1.9, 2.9, 4.4, 6.1, 8.4, 9.8` for `n` up to `10⁴`. **The extremal squared-distance is always a
  product of primes `≡ 1 (mod 4)`**: `k* = 25=5², 65=5·13, 325=5²·13, 1105=5·13·17, …` (verified through
  `N=20000`: `argmax r₂` runs `25, 65, 325, 4225=5²13², 5525=5²·13·17`). `max_{k≤N} r₂(k) = N^{(log2+o(1))/loglog N}`
  — the `n^{c/loglog n}` growth, with `c` from doubling `r₂` per new prime `≡1 (mod 4)`.
- **Upper bound is incidence-geometric.** Unit distances = point–unit-circle incidences. The bipartite incidence is
  `K_{2,3}`-free (two unit circles meet in ≤2 points ⟹ two centers share ≤2 unit-neighbors; verified: max common
  unit-neighbors `= 2` on a `15×15` grid at distance `√5`). Kővári–Sós–Turán ⟹ `O(n^{3/2})`; Szemerédi–Trotter
  (or Székely's crossing-number proof) refines to `O(n^{4/3})`.
- **The repo's small-`n`/lattice cluster is the same object at finite scale**: THM-412's lattice unit-distance
  layer *is* the `r₂` count above; THM-431/432's `u(21)=57` is its extremal small-`n` value. This note is their
  `n→∞` continuation.

## 3. Cross-thread placement (the point of the session)

**(a) The extremal distance is a smooth-number / "clock" object.** `r₂` is multiplicative with
`r₂(p^a)=4(a+1)` for `p≡1`, `r₂(p^a)=0` for odd-power `p≡3`; its maximizer over `k≤N` is the product of the
smallest primes `≡1 (mod 4)` (5,13,17,29,…). This is the *same shape* as the repo's LRC extremal objects — the
"divisor-complete cores" and clock towers are products of primes under residue constraints, maximizing a
multiplicative function over smooth numbers. **Honest status: a structural analogy** (both are
"maximize-a-multiplicative-function-over-smooth-numbers" problems), not a transfer — the multiplicative functions
(`r₂` vs the LRC clock density) differ. But it says the unit-distance lower bound lives in the repo's native
prime-product/residue-class world, and its extremal search is a smooth-number optimization the LRC machinery is
built for.

**(b) Unit distances belong on the CANCELLATION side of the repo's dichotomy — which is why the GMC move fails
here.** The repo's signature unlock (memory `gmc-lrc-same-positivity-manoeuvre`) is to re-express a hard count as a
*positive-definite* pairing (Hankel/autocorrelation), killing the cancellation wall. Unit distances *look* like a
candidate: `U(P) = ½·(autocorrelation of `P` against the circle measure `σ`)`. **But the manoeuvre does not
transfer, for a precise reason.** The spectral kernel is `σ̂(ξ) = J₀(2π|ξ|)`, the Bessel function, which
**oscillates and changes sign** (decay `|ξ|^{-1/2}`). So although the *direct* count `U(P)` is a sum of
non-negative indicators (an edge count, trivially `≥0`), the *spectral* bound — the only route to `n^{1+o(1)}` for
structured sets — pairs `|P̂|²` against a **sign-indefinite** kernel. That is exactly the repo's cancellation
regime: LRC covering (`|G_δ| = Σ ∏ĝ`, signed sinc weights, all-cancellation, no positive regime — boxeph
S228/S229) and DvdK's coincident-cycle stratum (signed `c^r`, cancellation — S101/S106). **Unit distances are a
positive count with a cancelling spectral certificate** — the count is combinatorially positive, but the tool that
would give the sharp bound faces the same oscillating-kernel wall that forced the LRC `χ`/topology route and the
DvdK Galois route. This is the honest reason the positive-definite manoeuvre — and hence a clean spectral
`n^{1+o(1)}` — is blocked, and it predicts that the sharp upper bound (if it follows the repo's pattern) needs a
*topological/arithmetic* certificate that survives the sign cancellation, not a positivity argument.

**(c) The upper-bound obstruction is unexploited translation-congruence.** All `n` circles are *congruent*
(radius 1) — translates of one circle, an abelian (translation-only) 2-parameter family. Szemerédi–Trotter treats
them as a generic 2-parameter curve family and cannot see the congruence; the believed truth `n^{1+o(1)} ≪ n^{4/3}`
is exactly the statement that congruence forces a much smaller bound. This is the dual of the distinct-distances
problem, where Guth–Katz (2015) *did* exploit the rigid-motion group `SE(2)` via the polynomial method; the reason
their method cracks distinct distances but not unit distances is precisely (b) — distinct distances is a
`max/counting` extremal with a favorable algebraic parametrization, unit distances a same-value count whose
spectral certificate cancels.

## 4. A concrete route-shaped statement

Combining (b)+(c): the sharp upper bound should come from an **arithmetic/topological certificate for the
translation-congruent incidence that is stable under the `J₀` sign cancellation** — the unit-distance analogue of
codex's LRC `χ`-criterion (THM-2047, cancellation-blind volume replaced by a topological count) and of THM-2067's
Galois orbit-product (manufactured valuation surviving unit-coefficient cancellation, S106). Concretely: for
structured `P` (lattice or algebraic), `U(P) = n^{1+o(1)}` is *provable* via the `r₂` number theory (§2), far
below `n^{4/3}`; the open content is the **structure-vs-randomness** step — that no set beats the structured
optimum. This is exactly the regime the repo's additive-combinatorics and cancellation-certificate threads target;
the missing object is a certificate that (i) is translation-covariant, and (ii) survives `J₀`'s sign changes.

## 5. Honest scope

- **Not a new bound.** §2 anchors are classical (`r₂` lower bound, Szemerédi–Trotter upper), here verified in-repo
  and connected to the existing THM-412/431/432 cluster. §3–4 are a *placement and route*, not a theorem.
- **What is genuinely useful:** opening the asymptotic problem in the repo; recording *why* the positive-definite
  manoeuvre fails (the `J₀` sign-indefinite spectral kernel), which prevents wasted effort and correctly files
  unit distances with LRC-covering and DvdK-resonant on the *cancellation* side; and identifying that the sharp
  bound wants a cancellation-surviving arithmetic/topological certificate (the repo's `χ`/orbit-product template),
  not a positivity argument.

Cross-links: THM-412 (lattice unit-distance density = the `r₂` layer), THM-431/432 (`u(21)=57`, Moser-lattice
optimum — the small-`n` values), memory `gmc-lrc-same-positivity-manoeuvre` (the manoeuvre that does *not*
transfer, and why), boxeph S228/S229 (LRC covering = all-cancellation, the sibling regime), S106 (Galois
orbit-product = a cancellation-surviving certificate template), S107 (the resonance dictionary — unit distances add
a fourth entry on the *cancellation* pole), codex THM-2047 (`χ`-criterion, cancellation-blind volume replaced by
topology — the analogue certificate). Script `04-computation/unit_distance_asymptotic_probe_deathstar_S108.py`
(+ `.out`). HYP-8945.
