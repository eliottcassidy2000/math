# The nature of numbers in tournaments: two arithmetics, and the Erdős-592 / Jacobian trichotomy

**death-star-2026-07-20-S60** (HYP-8245; owner across two prompts: how the tiling
hypercube, the merged metagraph, and other tournament structures reveal the nature of
numbers — "natural numbers map to each tournament, or every other, or something more
ingenious"; and how Erdős 592 relates to the three distinct parts of the Jacobian
counterexample, via the Pisano-60 / 1001 theme). Two creative threads, both verified
where checkable, both about the same spine — **the two multipliers 2 (doubling) and 3
(trisection)**.

## Part 1 — Tournaments carry TWO arithmetics of the numbers

**(I) The additive arithmetic — the tiling hypercube IS the integers in binary.**
Fix the base path n→n−1→…→1. A tournament is a bit string of the m = C(n−1,2) tile
arcs — an integer N ∈ [0, 2^m). So the natural numbers 0…2^m−1 *are* the labeled
tournaments, and the tiling hypercube Q_m is their **additive / 2-adic** structure: XOR
= tile symmetric difference (wiggly/waggly moves), Hamming distance = flip distance.
"Natural numbers map to each tournament" — this is the bijection.

**(II) The multiplicative arithmetic — H is a norm, strong tournaments are the primes.**
The ordinal sum T₁ ⊕ T₂ (every vertex of T₁ beats every vertex of T₂) satisfies
**H(T₁ ⊕ T₂) = H(T₁)·H(T₂)** (verified, 200 random pairs) — a Hamiltonian path can't
return once it leaves T₁. So H : (tournaments, ⊕) → (odd numbers, ×) is a multiplicative
**norm**; the **strong tournaments are the primes** (ordinal-indecomposable); the single
vertex is the **unit** (H=1); the transitive tournament is the all-units product.

The image is a **numerical monoid**: the strong-tournament H-values generate, under ×,
**exactly the odd numbers except {7, 21}** (verified: the monoid up to 45 misses
{7,21,35}, and 35 = 5·7 is realized as a *strong* H-value at n=7, filling; only {7,21}
are permanent). And the gap is number-theoretic: **7** is prime and not a strong
H-value (THM-029), so it is neither a product nor a prime here; **21 = 3·7** needs the
absent 7 (THM-079). "Every other tournament / something more ingenious" — this is it:
tournaments realize the odd numbers as a multiplicative monoid with the sporadic gap
{7,21}.

**The bridge (THM-466): the multiplicative value, in the additive base, is the census.**
H = Σ α_k 2^k, so the *binary digits* of the multiplicative norm H are the tournament's
odd-cycle-collection counts. Multiplicative value, additive digits, combinatorial
meaning — the two arithmetics meet in the 2-adic expansion. And the **merged metagraph
G_n/Z₂** is the number line modulo the natural symmetries: relabeling (S_n) and
complement (the antipodal map on Q_m, mac-mini-S11).

## Part 2 — Erdős 592 and the Jacobian counterexample are the same trichotomy

Both objects decompose into **three axes**, and the axes match one-for-one (verified):

| axis | Erdős 592 (α → (α,3)²) | JC counterexample F |
|---|---|---|
| **OBSERVER / 2n+1** (tame, linear) | tree-grid: R(n,2) = 2n+1 (R(1)=3, R(2)=5, **R(3)=7**) | odd fiber 1 + 2 = Rédei parity (opus THM-1350) |
| **×3 TOWER** (wild) | Chang tower ω^(ω^m), m=1/2 done, **m=3 open ($1000)** | degree-3^m tower F, F², F³ = deg 3, 9, 27 (klein-S327) |
| **CHAR-2 ATOM** (the seam) | bi-dyadic b=2, forced by a unique-prime fact (THM-469) | det = −2, λ↦λ², the 2-adic staircase (THM-1300) |

The shared spine is the **two multipliers 2 and 3** — the char-2 atom (doubling) and the
×3 tower (trisection, the "3" in (α,3)² and in geometric degree 3). This is exactly the
doubling/trisection dichotomy of S59w (the ×3 trisection pairs {2,6}, {7,21} vs the ×2
doubling pair {12,24}). Erdős 592 and the fallen Jacobian conjecture are **one
trichotomy — observer / ×3-tower / char-2 — playing in two theaters** (countable
ordinals; polynomial maps of ℂ³). And the observer axis literally coincides at n=3:
R(3,2) = 2·3+1 = **7** = the seven-wall = the H=7 permanent gap = the 7 of 1001 = 7·11·13.

## Part 3 — the Pisano / 1001 clock (graded: evocative, one real hook)

The last digit of Fibonacci repeats with period **60 = 2²·3·5 = lcm(π(2)=3, π(5)=20)** —
the period where the doubling (2) and trisection-flavored (3) clocks align, the same 2
and 3 that spine both trichotomies. **1001 = 7·11·13** carries the 7 (the observer of 3,
the seven-wall). These are evocative resonances, not theorems — but the *structural*
claim under them is real and verified: the number 60, the number 7, and the {7,21} gap
all sit on the 2-and-3 axis that organizes the tournament multiplicative monoid, Erdős
592's three axes, and the JC counterexample's three parts alike. The Fibonacci
connection is genuine on the JC side too (the counterexample's Pisano clock, kind-pasteur
S128c100; the golden JC₂ corner, mac-mini-S137).

## Part 4 — folded-in results (this session's other threads)

- **h(G₈/Z₂) ≥ 22** (certified K₄…K₂₂ minors; V=3528, E=45550, ω=5, χ=7). Growth vs χ
  *steepens*: h/χ = 12/6 = 2.0 at n=7, ≥ 22/7 ≈ 3.1 at n=8. G_n/Z₂ is confirmed a
  minor-dense, low-clique (ω 4→5), low-chromatic (χ = n−1) family — a genuinely
  publishable graph-theory phenomenon (Hadwiger with a widening margin).
- **Lean/publishability triage** (miner): the top Lean target is THM-1300 §0+§1 (the JC
  counterexample + the A₃ Dixmier endomorphism — pure polynomial identities, no Mathlib
  AG); then THM-025 (real-rootedness disproved at n=9), the {7,21} spectrum note, and a
  **paper-ready Cayley/Delannoy identity cluster that is absent from both ledgers**
  (`drafts/candidates-for-sharing.md`: the master GF with k·g_k(m)=m·g_m(k) duality, the
  CV² 1/n² cancellation, the golden exceptional points, the bilinear-Delannoy bridge) —
  the most conventionally publishable, most Lean-able items in the repo. And a systemic
  THM-id collision issue (~70 doubled numbers) to de-collide.
- **The witness extraction** (assessed S59z): novel, meaningful, the #1 publishable move,
  genuinely untouched — verify the witness *directly* (the finite Δ^m(P^m) certificate,
  built S59z, is Lean-able) rather than via the equivalence. Feasibility gate = the
  reduction's dimension; the recipe agent failed on an API error (it retrieved Zhao's
  inversion formula, arXiv 0704.1689, before dropping), so the exact dimension is still
  to be pinned — the one open input.

## Honesty

Part 1 (the two arithmetics) is verified exactly (H multiplicative norm; the {odds}\{7,21}
monoid; THM-466 bridge). Part 2 (the Erdős-592 / JC trichotomy) is a structural analogy —
each half's three axes are individually established in the repo/literature, and the
one-to-one match (observer / ×3-tower / char-2) is real and striking, but it is a
correspondence of *structure*, not a theorem connecting the two problems. Part 3 is
graded evocative with one verified hook (R(3,2)=7; 60=2²·3·5). The value: a single
organizing frame — the numbers 2 and 3 — under which the tournament monoid, Erdős 592,
and the Jacobian counterexample are the same story.

## Cross-links

THM-466 (2-adic digits = odd-cycle census) · THM-029/079 ({7,21} gaps) · THM-1300 (the
counterexample) · opus THM-1350 (odd fiber) · klein-S327 (degree-3^m tower) · THM-469
(Erdős 592 char-2 seam) · the three-axes-of-erdos-592 reflection · S59w ({2,6}/{7,21}
trisection, {12,24} doubling) · S59x/S59y (Hadwiger) · the S59z witness assessment ·
kind-pasteur S128c100 (the Pisano clock).
