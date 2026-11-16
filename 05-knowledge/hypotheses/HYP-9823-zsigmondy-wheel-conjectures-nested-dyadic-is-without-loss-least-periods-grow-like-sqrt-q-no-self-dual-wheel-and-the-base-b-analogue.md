---
id: HYP-9823
title: "Zsigmondy-wheel conjectures: (W1) for every N the wheel W(N) covers Z/N iff the odd parts of the orders cover Z/m (m the odd part of N), i.e. the nested-dyadic criterion of THM-5254 (B) is also necessary; (W2) the least a with 2^a q a covering period satisfies a(q) << sqrt(q), with a(3) = 3, a(5) = 4, a(7) = 4 proved; (W3) no single wheel is self-dual (serves both signs of one coefficient); (W4) in base b the least covering period is governed by the base-b analogue of THM-5254 (C), the missing moduli being Zsigmondy's exceptions (2 for b + 1 a power of two, 6 for b = 2)"
status: >
  OPEN (CONJECTURE), klein-2026-11-16-two3. Evidence: (W1) exact agreement at every decided period N <= 150 (base 2);
  (W2) the counts c(a, q) for q = 3, 5, 7 and the heuristic nu(2^i q) ~ log(2^i q); (W3) no self-dual cover at any decided
  period below 150 (24, 40, 48, 56, 60, 64, 72, 80, 84, 96, 112, 128, 140 fail; 120 undecided after 2e8 nodes); (W4) base 3 covers exactly at
  48, 72, 90, 96, 108, 120 below 120, with the modulus 2 missing. Packet two_three_wheels_synthesis_20261116 (WHEELS.md
  section 3, PREREGISTERED.md). What a proof of (W1) needs: that a non-nested choice of dyadic classes never helps the odd
  wheel; the coset argument of THM-5254 (B) shows it never helps in density. (W2) needs the primitive parts of 2^(2^i q) - 1
  for q = 11, 13, 17 factored. (W3) is a two-sheet exact-cover question; a density count does not exclude it.
session: klein-2026-11-16-two3
date: 2026-11-16
source: 05-knowledge/results/two_three_wheels_synthesis_20261116/WHEELS.md
related:
  - 01-canon/theorems/THM-5254-zsigmondy-wheels-covering-periods-for-k-2-n-plus-minus-1-reduce-to-the-odd-part-least-periods-24-80-112-the-classical-covering-sets-are-full-wheels-and-the-riesel-primes-meet-24-only-through-von-staudt-clausen.md
  - 05-knowledge/results/logarithmic_repairs_covers_20261116/CLOCKS.md
---
# HYP-9823 — Zsigmondy-wheel conjectures

**(W1) Nested is without loss.** THM-5254 (B) proves: cover of `Z/m` by the odd parts ⇒ cover of `Z/N`, and cover of
`Z/N` ⇒ `Σ 1/oddpart ≥ 1`. Conjecture: the first implication reverses. True for all `N ≤ 150`; the five periods where a
bounded search did not decide were settled by the `a = 1` identity, the density bound and a longer search.

**(W2) Least periods.** `a(q) = min{a : Σ_{i ≤ a} ν(2^i q) ≥ q}`; proved `a(3) = 3`, `a(5) = 4`, `a(7) = 4`, `a(q) ≤ q − 1`.
Conjecture `a(q) ≪ √q`, from `ν(2^i q) ≈ log(2^i q)` (the primitive part `Φ_{2^i q}(2)` has about that many prime
factors on average).

**(W3) Self-dual wheels.** By THM-5254 (D) only even-order primes serve both signs; a self-dual wheel would have residues
`r_p` with both sheets `{r_p}` and `{r_p + ord_p/2}` covering `Z/N`. None found below 150 among decided periods; the known
Brier numbers use two sub-wheels of one period (odd-order primes serving one sign each), consistent with (W3).

**(W4) Other bases.** The wheel of base `b` lacks the modulus `n` exactly at Zsigmondy's exceptions (`n = 2` when `b + 1`
is a power of two; `n = 6` for `b = 2`). Conjecture: the least covering period in base `b` is the first level at which the
base-`b` analogue of `c(a, q)` reaches `q` for some odd prime `q` not dividing `b − 1`; for `b = 3` that period is 48.

**Why they matter.** They turn the existence of Sierpiński/Riesel/Brier covering sets into a counting statement about
primitive prime divisors, with Zsigmondy's single exception as the reason for the historical period 24; and they answer
the owner's Riesel-versus-24 question structurally (THM-5254 (G)) rather than by a lattice.
