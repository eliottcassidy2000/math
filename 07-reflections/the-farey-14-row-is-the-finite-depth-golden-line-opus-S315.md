# The Farey-14 row is the finite-depth golden line — the six speeds, the primes, and the Fibonacci rate in the radius-7 resonance framework

*opus-2026-07-15-S315; owner directive: "see how 6 speeds relate to prime
numbers and the rate of the fibonacci sequence (it feels 6 periodic, according
to past work in the repo) and zeckendorf."*

## The decode

The repo's proven 6-periodicity is the **mediant-attainer trichotomy**
(mac-mini HYP-4572 + opus HYP-4516, LRCBinderInfeasible.lean): the canonical
family is a gap member iff N ≡ 1 (mod 6), through the mod-30 binder gate. The
Fibonacci lineage is klein-S124 ("Fibonacci is the covering min's FOIL, not
its lever"), klein-S269 (LRC = 2-term Ostrowski), and the Farey–Fibonacci
bridge reflections (S169). The new frame that ties them to the CURRENT
frontier is THM-863/864's resonance functional

> Y*(x₁, x₂) = min₁≤q,p≤13 |q·x₂ − p·x₁|,

which governs the E-restriction error of Hunter tree edges at the seven-comb
wall (THM-864: err ≤ 13κρ/(y(p+q−1)) + dust).

## Fact 1 — the Fibonacci edge law (exact, one-line proof)

**Y*(F_n, F_{n+1}) = F_{n−7}, attained at (q,p) = (F₇, F₈) = (13, 21).**
*Proof:* d'Ocagne, F_k F_{n+1} − F_{k+1} F_n = (−1)ⁿ F_{n−k}; the best
approximants of F_{n+1}/F_n = [1;1,…,1] are exactly the Fibonacci ratios, so
the minimizing q ≤ 13 is the largest Fibonacci ≤ 13, i.e. F₇. ∎ (Verified
exactly n = 9..25.) Golden-ratio speed pairs are thus PROVABLY non-resonant
at the exponential rate φⁿ — "the rate of the Fibonacci sequence" enters the
radius-7 machinery as an edge-quality guarantee.

## Fact 2 — but the FINITE-depth extremal is the Farey-14 row, six per block

The landscape scan (x = 2003, all ratios in (1,13)) puts EVERY one of the
top-10 Y*-classes at p/14 with gcd(p, 14) = 1, achieving Y* = ⌊x/14⌋ = 143
exactly — while the golden line (φ+1, φ+2, 2φ) reaches only ≈ 0.0554·x.

- *Upper bound:* by Dirichlet/three-distance, min_{q≤13} ||q·α|| ≤ 1/14 for
  every α, so Y* ≤ x/14 + O(1) always.
- *Attainment:* near α = p/14 with p a unit mod 14, ||q·p/14|| ≥ 1/14 for all
  1 ≤ q ≤ 13 — the just-out-of-range rationals saturate the ceiling.

> **At resonance depth Q = 13 the extremal non-resonant directions are the
> (Q+1) = 14th Farey row's units — and |（Z/14)*| = φ(14) = 6: SIX extremal
> ratio classes per block, {1,3,5,9,11,13} mod 14.** The owner's "6 speeds"
> at this frontier: six optimal edge-ratio classes per 14-block, recurring
> 14-periodically up the ratio ladder, i.e. 6-per-period. The primes: every
> prime > 7 is a unit mod 14 (odd, not 7), so the extremal classes are
> exactly where the primes live; and mod 6 those units read {1,3,5,3,5,1} —
> the ±1 mod 6 prime skeleton appears with the two 3's contributed by the
> 7-column (3·? — the classes 3, 9 ≡ 3 mod 6 are the ones divisible by 3;
> honest note below).

**The foil/lever duality sharpened (klein-S124 vindicated and completed):**
Fibonacci/φ is the extremal of the DEPTH-∞ functional (Hurwitz/Markov:
lim inf q||qα|| maximized on the noble line — Zeckendorf/Ostrowski
numeration is its bookkeeping). The DEPTH-13 functional — the one the
LRC(14) machinery actually uses — is maximized instead on the problem's OWN
denominator row p/14. Fibonacci is the foil (0.0554x); the 14th Farey row is
the lever (0.0714x = x/14). **LRC(14) is self-tuned: its optimal Hunter
edges are the ratios with its own denominator.**

## Fact 3 — the golden sawtooth is Pisano-28

ρ(F_n, F_{n+1})'s sawtooth kernel is a function of (F_{n+2} mod 13,
F_{n−1} mod 13), hence 28-periodic in n (π(13) = 28); exact zeros occur at
n = 13, 14 (s ≡ 12). The "6-periodic Fibonacci rate" is NOT this: the
relevant Fibonacci sixes are φ⁶ ≈ L₆ = 18 (near-integer) and π(4) = 6 —
neither is load-bearing here (checked).

## Honest ledger

- The TWO sixes — the mediant law's N ≡ 1 (mod 6) and the totient six
  φ(14) = 6 — are DIFFERENT structures (moduli 6 = 2·3 vs 14 = 2·7). Both
  are real; their identification is NOT established and should not be
  assumed. What they share: both are unit-group counts (φ(6)·… the mediant
  gate is a mod-30 unit condition; the edge extremals are mod-14 units) —
  "units of small even moduli" is the common shape, logged as a lead, not a
  law.
- Numerology killed: φ⁻⁶ = 0.0557281 vs the AP floor 2833/50700 = 0.0558777
  — not equal (Δ = 1.5e-4). Reported per the BH discipline: 1 candidate
  identity scanned, 0 hits.
- Zeckendorf's role here is organizational: Ostrowski numeration is the
  bookkeeping of the resonance ladder (each ratio's CF = its resonance
  signature); Zeckendorf is its golden instance and governs the DEPTH-∞
  regime only. klein-S269's "LRC = 2-term Ostrowski" lives on the other side
  (the time/return dynamics), untouched by this.

## What to do with it

1. Hunter tree construction should PREFER edges with ratio near units/14 —
   a concrete, checkable optimization of the THM-863 tree choice (the
   per-edge error budget gains the full 1/14-rate).
2. The p/14-plateau width and its interaction with the sheet canon (exact
   ratio p/14 means 14x₁ ≈ qx₂-relations — the 14-sheet!) is the next exact
   question: the extremal edges sit just off the 14-dilate resonance, the
   problem's own modulus appearing as both the wall and the lever.
3. Cross-reference: THM-863/864, klein-S124/S269, HYP-4572/4516 (the mod-6
   mediant gate), the S169 Farey–Fibonacci bridge.
