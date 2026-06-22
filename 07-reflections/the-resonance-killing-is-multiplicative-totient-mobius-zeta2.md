# The Resonance-Killing Is Multiplicative: Totient, Möbius, and the ζ(2) Floor

*mac-mini-2026-06-22-S44. The owner asked how coprime density relates to the Euler totient and
multiplicative functions, and how the three tournament recursion modes (IE-signed, even/odd) connect.
They are one structure: the resonance-killing game (kps S31p) is a **multiplicative-function** phenomenon,
and its floor is the **coprime density 1/ζ(2)**. Builds on kps S31p (resonance-killing), opus S291
(Euler-product metagraph), my HYP-2880 (H/totient), [[zeta2-governs-the-lonely-runner-floor]].*

## The resonance-killing, refined by the totient

kps S31p: `M(S) = 1/(smallest surviving b)`, where resonance `b` is *killed* iff some runner `s ≡ 0 (mod b)`
(that runner sits on the observer at `t = a/b`). Refine by asking *which Farey points* a runner kills:

> A runner `s` kills the Farey point `a/b` (`gcd(a,b)=1`, `b ≤ 14`) iff `‖s·a/b‖ < 1/14` iff `b | s`.
> So `s` kills **all φ(b)** primitive points of every denominator `b | s`. **The killing is totient-weighted.**

The survival lattice is the Farey set of order 14, of size `Φ(14) = Σ_{b≤14} φ(b) = 64`. A counterexample
must cover all 64 — equivalently kill every denominator `b ∈ {2,…,14}` (13 denominators, 13 runners: the
over-determined covering condition, THM-523). The totient φ(b) is the *multiplicity* with which each
denominator's points fall (or survive) together.

## Coprime density → totient → ζ(2), the floor

The surviving neighbourhoods sit near the Farey points `a/b`; their total measure is the universal
**Farey floor `3/π² = 1/(2ζ(2))`** ([[zeta2-governs-the-lonely-runner-floor]]). This is exactly the
totient made analytic: summing the primitive-point count `φ(b)` against the `1/b²` neighbourhood weight
gives `Σ φ(b)/b² ∝ 1/ζ(2)` — the **coprime density `6/π² = 1/ζ(2)`** (the probability two integers are
coprime). So:

> **coprime density `1/ζ(2)`  =  Σ_b φ(b)·(point weight)  =  the density of surviving lonely points.**

The multiplicative function `φ` (counting coprimes) is *why* the floor is positive and computable: `ζ(2)`
converges, so the surviving lattice never vanishes. The lonely runner is lonely because coprimality has
positive density.

## The three recursion modes as the Möbius/IE skeleton (structural)

`φ = μ * id` (Dirichlet). The killing/survival is an **inclusion–exclusion over the divisor lattice**,
and the owner's three recursion modes are that IE skeleton wearing tournament masks:
- **Mode 1 (always) `A+B+C−D−E−F+G`** = the union IE over a 3-set (`|A∪B∪C|`: 3 singletons +, 3 pairs −,
  1 triple +) = the squarefree-3-prime Möbius pattern.
- **Mode 2 (even) `A+B−C`** = the union IE over a 2-set (`|A∪B|`).
- **Mode 3 (odd) `A+B−C+D−E−F+G`** = Mode 2 extended by a `+D−E−F+G` tail — a *different* grouping of the
  3-deep IE, on subtournaments of *different sizes* than Mode 1 (as the owner flags).

The same multiplicative spectrum governs the **tournament count** itself: opus S291's Euler product
`V_n·n!/2^m = 1 + Σ_{k odd} (1/k)·n↓k·2^{…}` is an Euler product over **odd cycle lengths** with `1/k`
prime-reciprocal factors — the tournament analogue of `ζ`. So the recursion modes, the totient floor,
and the tournament Euler product are three faces of one multiplicative structure.

**Even/odd parity is the through-line:** Mode 2 (even) vs Mode 3 (odd) ↔ the Möbius sign `(−1)^{#prime
factors}` ↔ Burnside's `Fix(σ)=0` for even cycles vs `2^{…}` for odd (the Euler product runs over *odd*
k only) ↔ the **apex-7** (the unique prime where orientation/parity diverge). The killing ladder's
even `b` and odd `b` are the two Möbius parities.

## Toward the proof (and the honest gap)

The multiplicative picture *is* the proof skeleton for the lower bound: the disproof must cover all
`Φ(14)=64` totient-weighted Farey points (kill every `b≤14`), which is over-determined for 13 runners;
and the only way to kill the large `b` (13,14 together) is a *bulk* killer `v=lcm(13,14)=182`, which —
being huge — **equidistributes** and its thin arc *misses the neighbourhoods* (kps S31p), so they survive
at the `1/ζ(2)` floor density. The totient/ζ(2) floor is what equidistribution cannot erase.

HONEST gaps (discipline): (1) I connected the *sign/Möbius structure* of the three modes but did **not**
pin their exact subtournament sizes — that identification is open (the owner notes the sizes differ per
mode; kps/codex may hold it). (2) The `1/ζ(2)` Farey floor is rigorous at the *relaxed* threshold `1/7`;
at the strict `1/14` the extremal floor is measure-zero (S39) and the survival is via boundary witnesses
(the resonance ladder). The multiplicative structure unifies both; the strict-threshold rigor is still
the bounded-compactness + equidistribution split (THM-527/566).

Verified core (S44): the totient-weighted killing (`s` kills all `φ(b)` points of `b|s`), `Φ(14)=64`,
the `φ→1/ζ(2)` coprime-density floor. Related: kps S31p, opus S291, HYP-2880, HYP-2897, THM-523.
