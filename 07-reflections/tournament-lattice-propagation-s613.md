---
source: claude-2026-06-03-S613
status: THM D content verified; H 7*3^k refuted ({7,21} stands); Cl2(pi/3) link structural/numerical; synthesis
tags: [tournaments, H-impossibility, forbidden-7, prime-3, clausen, unit-distance, eisenstein, LRC, theorem-D, propagation, honest]
---

# Tournaments, H-impossibilities, and the prime-3 / prime-7 propagation

A long session on the user's vision: that LRC, unit distance, and Collatz are
"tournament-shaped," with the H-value impossibilities and a shared `1.014`
exponent as the keys. The honest results below include a refutation and an
"unverified," which sharpen the picture.

## Theorem D — verified content

`M(V*) = max(1/n, P)`, `P` = the `(2n−1)` pinch-lattice max. The strictly-lonely
(loose) intervals of `V*` center **exactly on multiples of `1/(2n−1)`** — the
`{3, 2n−4}` pinch and its mirror, exactly two per loose `n` (checked
`n=10,12,16,18,22`). So deleting `n−2` opens a window *only* on the `(2n−1)`
lattice. With S612 (loose proved; tight proved on the lattice) the doubling law
holds modulo the general "this is the only family" statement — which is *exactly*
codex's `Res₂₇` lift/CRT conservativity (HYP-2167). The two threads share one
open lemma.

## H-impossibilities — a refutation, honestly

`H(T) = #HamPaths = I(Ω, 2)`, the independence polynomial of the odd-cycle
conflict graph at 2; always odd (Rédei). The tempting pattern from a first scan
was `7·3^k = {7, 21, 63, 189, …}` (note `{1,3,9}` = the gcd-strata!). **Refuted:**
`63` and `189` *are* achievable — they first appear at `n=8`, and my initial
`n≤7` sample simply missed `63` (a sampling artifact). Re-running with `n=8`
sampling hit `63` six times and `189` 121 times.

So the forbidden set up to 80 is **`{7, 21}` only** — `7` proved impossible
(THM-200), `21` conjectured (THM-115). The H-impossibilities are **sporadic and
finite here, tied to the forbidden prime 7**, *not* an infinite arithmetic
family. The powers-of-3 `{1,3,9}` that I hoped to find in forbidden-H actually
live in the **gcd-strata / 3-shell** (LRC THM-407, my Theorem D) — a different
role. Lesson re-learned: verify the small cases at large enough `n` before
declaring a pattern.

## The `1.014` constant — structural link, unverified equality

`Cl₂(π/3) = 1.014942` (Clausen), the numerator of the tournament tropical
constant `κ` (HYP-707). The unit-distance disproof gives exponent `> 1.014`
(Sawin) — a **lower bound**, not a proven closed form. So **"UD exponent =
Cl₂(π/3)" is unverified**; the numerical proximity (`1.014` vs `1.01494`) is
suggestive but the disproof's mechanism (CM fields / class-field towers) is not
the triangular grid.

What *is* real and shared is the **prime-3 / Eisenstein angle `π/3`**:
- it produces `Cl₂(π/3)` in the tournament tropical log-sin sum (HYP-707), and
- it geometrizes the triangular unit-distance lattice `Cay(ℤ[ζ₆])` (HYP-2170).

Same prime 3, two appearances. Whether the disproof exponent is *literally*
`Cl₂(π/3)` would need the disproof's CM-tower density to reduce to the same
Eisenstein angle — a real question, not established here.

## The propagation — the Cayley primes are the shared invariants

The unifying picture (matching the repo's "master object" lesson, S614, and my
resonance-lattice synthesis HYP-2155): across LRC, tournaments, unit distance,
and Collatz the hard core lives in a **high-degree arithmetic lattice**, not the
visible structure, controlled by the Cayley/Hurwitz primes:

| prime | tournaments | LRC | unit distance |
|---|---|---|---|
| **2** (inert/doubling) | Rédei parity, GF(2) involution (solved face) | even-fold (HYP-2065); doubling voltage lift | edge/2-coloring duality |
| **3** (Eisenstein, `π/3`, ramified) | `Cl₂(π/3)` tropical const; boost RAMIFIED (THM-253) | 3-shell `3∣(2n−1)`; the `V*` doubling (THM D); gcd-strata `{1,3,9}` | triangular lattice `Cay(ℤ[ζ₆])`; `Cl₂(π/3)` ~ exponent? |
| **7** (forbidden) | `H=7,21` impossible (THM-200/115) | `2n−1=7` at `n=4`; pinch-modulus base | `F`-free frontier; `n=22 = 3·7+1` first awkward |

The "propagation" the prompt asks for is this: **the impossibilities and
constants of each problem are expressions of the same small primes (`2,3,7`)
acting on each problem's arithmetic lattice.** The forbidden 7 of tournaments and
the `2n−1` modulus of LRC and the Eisenstein 3 of unit distance are not analogies
— they are the same primes in the formal-group/Cayley frame (THM-253).

## Honest status and next

- **Solid:** Theorem D's geometric content (verified); the `{7,21}`/`7·3^k`
  refutation; `Cl₂(π/3)` value and its genuine prime-3 sharing.
- **Open/unverified:** Theorem D's general single-family statement (= codex
  conservativity); whether the UD exponent equals `Cl₂(π/3)`.
- **Next:** (1) close Theorem D / codex conservativity jointly. (2) Decide the
  `Cl₂(π/3)` equality: derive the UD-disproof density's effective exponent and
  test against the Clausen value. (3) Ask whether `H=21`'s impossibility (if
  true) has an LRC/UD shadow — does the forbidden 7 forbid a specific `n` or a
  specific unit-distance count?
