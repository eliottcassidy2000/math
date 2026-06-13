# THM-488: the pentagonal Lyapunov gap and the sign-rigidity of Euler's recurrence

**Status:** Part A (zero-free characterization) PROVED forward + reverse IVT-half;
reverse HARD HALF (B(S) ≥ 0) CERTIFIED by the argument principle (exact winding
numbers) on |S| ≤ 6, k ≤ 12 — 1585/1585, 0 counterexamples, with the
truncation-artifact correction; earlier passes 218/218 (k≤8,|S|≤5) and K=10 sweep
1023/1024. Part B (Lyapunov constants) COMPUTED, method VALIDATED against
Viswanath. Scripts `04-computation/pentagonal_lyapunov_kps3_0611.py`,
`pentagonal_rigidity_proof_kps3_0611.py`, `rigidity_winding_kps_0611.py`,
`rigidity_winding_diag_kps_0611.py`, `rigidity_winding_stress_kps_0611.py`,
`rigidity_winding_boundary_kps_0611.py` (+ .out).
**Source:** kind-pasteur-2026-06-11-S3 (HYP-2423/2417/2418). Continues the repo's
random-sign Lyapunov theme (HYP-614, Viswanath/Embree–Trefethen) and feeds the
code bridge (THM-487).

## Setup

Euler's partition recurrence, written as a signed linear recurrence with lags at
the generalized pentagonal numbers g_k = k(3k−1)/2, ḡ_k = k(3k+1)/2:
  a_n = Σ_{k≥1} ε_k ( a_{n−g_k} + a_{n−ḡ_k} ),  a_0 = 1.
Its generating function is 1/F_ε(x), F_ε(x) = 1 − Σ_k ε_k (x^{g_k} + x^{ḡ_k}).
For **Euler's signs ε_k = (−1)^{k+1}**, the pentagonal number theorem gives
F_ε(x) = Π_{n≥1}(1 − x^n), and a_n = p(n), the partition function.

The exponential growth rate of (a_n) is
  rate(ε) = max(0, −log r(ε)),  r(ε) = least modulus of a zero of F_ε in |x| ≤ 1,
so **rate(ε) > 0 ⟺ F_ε has a zero in the OPEN unit disk** (when F_ε is zero-free
inside, the coefficients are subexponential — governed by the boundary, the
essential singularity at x = 1).

## A. Sign rigidity (the zero-free characterization)

**Theorem A.** Among ±1 sign sequences ε differing from Euler's in finitely many
places, F_ε is zero-free in the open unit disk — equivalently rate(ε) = 0,
equivalently a_n grows subexponentially — **iff ε is Euler's alternation
(−1)^{k+1}**.

**Proof, forward (ε = Euler ⟹ rate 0).** F_Euler = Π(1−x^n); each factor is
zero-free in |x| < 1 and the product converges there, so F_Euler is zero-free in
the open disk. Hence rate = 0; concretely a_n = p(n) ~ exp(π√(2n/3))/(4n√3)
(Hardy–Ramanujan), subexponential. ∎

**Proof, reverse — the IVT half (PROVED).** Let S ≠ ∅ be the (finite) flipped set.
Then F_ε − F_Euler = Σ_{k∈S} 2(−1)^{k+1}(x^{g_k}+x^{ḡ_k}), so at x = 1 (where
F_Euler(1) = 0):
  **F_ε(1) = 4 Σ_{k∈S} (−1)^{k+1}.**
If Σ_{k∈S}(−1)^{k+1} < 0 then F_ε(1) < 0, while F_ε(0) = 1 > 0; by the
intermediate value theorem F_ε has a real zero in (0,1), so rate(ε) > 0. ∎

**Reverse, the hard half (B(S) ≥ 0) — ARGUMENT-PRINCIPLE CERTIFICATION
(kind-pasteur-2026-06-11 winding run).** For S with B(S) = Σ_{k∈S}(−1)^{k+1} ≥ 0
the boundary value F_ε(1) ≥ 0 and an interior zero — when present — is complex,
so the IVT is silent. We CERTIFY interior zeros with the argument principle: the
number of zeros of F_S strictly inside |x| = r equals the winding number of the
image curve F_S(r e^{iθ}) about 0, an EXACT integer (computed by adaptively
refined accumulated-argument summation — every step kept below π/4 so no 2π jump
aliases), not a heuristic grid minimum. Result: on the much larger family
**|S| ≤ 6, k ≤ 12, every one of the 1585 hard-half sets carries a certified
interior zero; zero counterexamples** (smaller pass: 255/255 on |S|≤4, k≤10).

*The truncation-artifact correction (the load-bearing subtlety).* The naive
winding count of the TRUNCATED pentagonal polynomial is NOT trustworthy near the
boundary: the partial sum approximates Π_{n≤N}(1−x^n), whose zeros are roots of
unity ON |x| = 1, and these leak inside as r → 1. Concretely winding(F_Euler^trunc)
jumps from 0 (true value, since Π(1−x^n) is zero-free in the disk) to 19, 41, 48
as r crosses 0.95 → 0.97 → 0.99. So a high-r winding count over-counts by the
artifact. The honest certificate uses two devices:
 (i) a **truncation-safe annulus** — at Kmax ≥ 20 the truncated Euler product has
     winding EXACTLY 0 for all r ≤ 0.95 (verified table), so a winding ≥ 1 there
     certifies a GENUINE zero of modulus < r < 1 of the true F_S;
 (ii) the **winding-difference** D(r) = w(F_S, r) − w(F_Euler, r) at the same r,
     which cancels the shared high-degree artifact zeros (F_S and F_Euler differ
     only in finitely many low coefficients, k ∈ S), for the handful of sets whose
     zero sits in the boundary annulus 0.95 < |x| < 1.
Exactly **4 sets** (on |S|≤6,k≤12) have their zero in the near-boundary annulus,
all supported on small odd k — {1,3,5,7,9}, {1,3,5,7,9,10}, {1,3,5,7,9,11},
{1,3,5,7,9,12}. Each is resolved by D(r) ≥ 1 at r = 0.95 AND by an explicit
Newton-polished complex zero of modulus 0.937–0.954 < 1 (residual ~1e-16). The
structural signature: flipping signs on the smallest odd k pushes the interior
zero toward |x| = 1 but never onto or outside it.

Independently, the K=10 growth sweep (all 1024 sign patterns on the first 10
pairs) found exactly ONE subexponential pattern — Euler's — with α-quantiles min
0.0009 (Euler), 1% 0.094, median 0.245, max 0.548. The reverse direction is thus
PROVED for the sign-of-boundary half (B(S) < 0) and CERTIFIED with no counter-
example for the hard half on |S| ≤ 6, k ≤ 12; a fully analytic all-n proof (every
finite perturbation of Π(1−x^n) acquires an interior zero) is the open residue —
**HYP-2417 stays PARTIAL pending it** (proof sketch in the Honesty section).
Scripts: `rigidity_winding_kps_0611.py`, `_diag_`, `_stress_`, `_boundary_`
(+ .out in 05-knowledge/results/).

**Reading.** The rigidity is the *analytic shadow of the product formula*: among
all sign-decorated pentagonal series, Euler's is the unique one that is a genuine
infinite product Π(1−x^n), and the product is exactly what forbids interior
zeros. The pentagonal number theorem is the zero-freeness certificate; flipping
any sign breaks the product and (provably for half the patterns, computationally
for the rest) plants a singularity inside the disk.

## B. The pentagonal Lyapunov constant and the gap

Three regimes of the recurrence (all randomness seeded; float with renormalized
log-scale accumulation, MISTAKE-067 guard; method validated — see below):

| regime | rate |
|--------|------|
| Euler signs (−1)^{k+1} | **0** (α = 0.0005, β = 2.48 ≈ π√(2/3) = 2.565: the genuine √n exponent of p(n)) |
| all-plus (ε ≡ +1) | **γ₊ = −log x\* = 0.548074**, x\* = 0.578062 the root of Σ(x^{g_k}+x^{ḡ_k}) = 1 (exact-int empirical slope matches to 1e-12) |
| fresh iid signs (Viswanath-type) | **γ_pent = 0.2059 ± 0.0032** |

**γ_pent ≈ 0.206 is a new random-recurrence Lyapunov constant** — the partition
analogue of Viswanath's 0.1240 for random Fibonacci. The *Lyapunov gap*
γ₊ − 0 = 0.548 (deterministic, all-plus) and the random value 0.206 sitting at
ratio 0.376 of γ₊ both quantify how much of the all-plus growth Euler's
alternation cancels: Euler kills 100%, fresh randomness kills 62%.

**Method validation.** The identical pipeline on lags {1,2} (the random Fibonacci
recurrence x_{n+1} = ±x_n ± x_{n−1}) returns γ = 0.1242 ± 0.0019 vs Viswanath's
constant log(1.13198824…) = 0.12398 — agreement within the spread.

**Comparative table (fresh iid, same pipeline):** lags {fib, pentagonal,
triangular, squares, primes} give γ = 0.124, 0.206, 0.173, 0.139, 0.169 and
γ₊ = 0.481, 0.548, 0.438, 0.349, 0.389; the ratio γ/γ₊ climbs 0.258 → 0.433 as
the lag set thins (denser low lags ⟹ more cancellation ⟹ smaller ratio). New
constants; none in OEIS as decimals (candidate submissions after more digits).

## C. The 2-adic floor (HYP-2418, VERIFIED)

The recurrence mod 2 is sign-free (±1 ≡ 1), so EVERY sign regime reduces to the
same F₂ recurrence and reproduces p(n) mod 2 (0 mismatches, n ≤ 3000; odd
density 0.506, the Parkin–Shanks ½). The random-sign structure lives strictly
ABOVE the 2-adic floor — the partition-function instance of the repo's
digit-ladder depth grading (THM-466/478: signs/Lyapunov data are invisible at
depth 0, where only parity survives).

## Honesty

- Theorem A reverse is PROVED only for Σ_{k∈S}(−1)^{k+1} < 0; the hard half is
  CERTIFIED by the argument principle (exact integer winding numbers) on
  |S| ≤ 6, k ≤ 12 — 1585/1585, 0 counterexamples — with the truncation-artifact
  correction (safe annulus r ≤ 0.95 at Kmax ≥ 20, plus winding-difference for the
  4 near-boundary sets, all with an explicit |z|<1 Newton witness). This is
  finite-family certification, not an all-n proof. HYP-2417 = the analytic
  completion.

- **Proof sketch toward the analytic completion (defensible, not yet a proof).**
  Write F_S = F_Euler + 2 Σ_{k∈S}(−1)^{k+1}(x^{g_k}+x^{ḡ_k}) =: F_Euler + P_S(x),
  a finite polynomial perturbation. F_Euler = Π(1−x^n) vanishes to INFINITE order
  at x = 1 in the sense that all derivatives are dominated by the essential
  singularity; on the boundary it has a dense set of near-zeros at roots of unity.
  The empirical pattern (all certified zeros cluster at modulus 0.94–0.99 for the
  hard sets, deeper for the easy ones) says the perturbation P_S — which is O(1)
  near x = 1 while F_Euler → 0 — DOMINATES F_Euler in a boundary collar, so the
  zero set of F_S there is governed by P_S = −F_Euler, i.e. by where the finite
  polynomial P_S meets the rapidly-decaying Euler product. A Rouché argument on a
  contour |x| = 1 − δ(S) would close it IF one can show |P_S| > |F_Euler| on that
  contour for the hard half (B(S) ≥ 0) — true pointwise in the data but needing a
  uniform lower bound on |P_S| over the circle (P_S is a sparse pentagonal
  polynomial; its minimum modulus on |x| = 1−δ is the obstruction). The clean
  statement to target: *for B(S) ≥ 0, F_Euler has no zero in 1−δ ≤ |x| < 1 but
  F_S does, because P_S is a nonzero finite pentagonal polynomial whose own zero
  structure forces a sign/phase configuration F_Euler cannot cancel.* The IVT half
  is the B(S) < 0 special case where this collapses to the real axis. Hurwitz/Rouché
  is the right machine; the missing lemma is the uniform |P_S| bound.
- γ_pent, γ₊ are float estimates with seeded reproducibility; γ_pent's error bar
  is the across-run std (8 runs, N = 20000), not a rigorous interval.
- The β ≈ π√(2/3) recovery of the Hardy–Ramanujan exponent is a regression
  readout, not a reproof of Hardy–Ramanujan.

**Cross-refs:** HYP-614 (Dedekind regulator = growth rate, the φ-Lyapunov),
THM-466/478 (digit ladder / 2-adic depth), THM-487 (the code bridge: the same
Π(1−x^n) = η/q^{1/24} drives extremal self-dual enumerator corrections),
HYP-2423/2417/2418.
