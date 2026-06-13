# The cubic lens: where 1729 is a theorem and where it is a coincidence

*kind-pasteur-2026-06-10-S1. Dispatched by the human: "this repo came across 1729 a few
times but was unable to flesh it out … look at the mathematical progress we have made
through the cubic lens." Builds on THM-434 (Moser rosettes), HYP-2306 (the severed
tournament spine), HYP-2307/THM-438 (R(p) → e), and lands THM-462, THM-463, LEM-003,
HYP-2370, HYP-2371. Every claim below was adversarially verified by an independent
agent this session (36/36 confirmed).*

## The question 1729 actually poses

A number that keeps reappearing poses exactly one question: **is the recurrence a
mechanism or an accident?** This project has now answered it for 1729 three times, and
the three answers are different in kind — which is the real finding.

1. **Tournament lane: ACCIDENT** (HYP-2306, 2026-06-07). `r(11) = H(T_11)/|Aut| = 1729`
   is the last gasp of small-prime smoothness in `H(T_p)`; the splitting law dies at the
   very next genus-1 Paley prime. Tested, refuted, closed.

2. **Taxicab × Moser lane: MECHANISM — now a theorem** (THM-463, this session).
   `x³ + y³ = (x + y) · N(x + yω)` factors through the *same* Eisenstein norm form whose
   representation count `r_E` is THM-434's rosette law `12 + r_E(t)`. Two-cube
   representability collapses to ≤ τ(n) divisor checks (the parity side-condition is
   automatic — a small sharpening of the textbook criterion), and for primitive
   representations every cofactor is `3^{0,1} ×` (split primes): the cofactors are
   FORCED onto the split axis, the exact axis on which THM-434's record rosettes live.
   1729 = 7·13·19 is the product of the three smallest split primes, `B = τ = 8`, so it
   simultaneously heads the taxicab list (Ta(2)) and the rosette record list (60 unit
   vectors) — one mechanism, two appearances. The honest boundary makes it a theorem
   rather than a slogan: taxicab numbers are NOT all completely split (1952 of the 5464
   doubly-primitive ones below 10¹² carry inert primes), but inert content is provably
   confined to `gcd(d₁, d₂)` with full multiplicity — it never touches a cofactor.

3. **Ramanujan's deathbed lane: MEMBERSHIP — also a theorem** (HYP-2370). 1729 is not
   merely "adjacent to" Ramanujan's near-miss family `a³ + b³ = c³ ± 1` (the
   `1 − 82x − 82x² + x³` generating functions): it IS the family's `n = −1` term. The
   recurrence has unit leading and trailing coefficients, so it runs backwards, and the
   backward extension hits `(−9, 12, 10)`: `9³ + 10³ = 12³ + 1`. Proved for all `n ∈ ℤ`
   two independent ways (order-27 Cayley–Hamilton on the Kronecker cube; factor
   `t³−82t²−82t+1 = (t+1)(t²−83t+1)` + Vandermonde). So Hardy's taxi number sits one
   unit above Klein's `1728 = j(i)` *and* one step behind `(135, 138, 172)` — and the
   "unit" is literal: a near-miss says the Eisenstein factorization misses a Fermat cube
   by an element of the order-6 unit group of `ℤ[ζ₆]`, FLT(3) says the gap never closes,
   and Ramanujan's family realizes the closest approach infinitely often.

The meta-lesson sharpens MISTAKE-028/036/055 ("a pattern is not a law"): **a number's
meaning is not in the number — it is in which of its appearances share a mechanism.**
The integer 1729 carries one theorem, one membership, and one accident, and only the
factorization 7·13·19 tells you which is which.

## The cubic degree ladder (why the lens has exactly three rungs)

The session's results arrange themselves by degree of collapse:

- **Two cubes = a curve that factors.** `x³+y³` splits off `(x+y)`; representability is
  a divisor property (THM-463); everything is decidable by factoring. This is the rung
  where 1729 lives.
- **Three cubes = a surface that doesn't.** `x³+y³+z³` is irreducible (the Fermat plane
  cubic is smooth; reducible cubics are singular — one line), so there is no norm-form
  collapse: `x³+y³+z³ = k` is integral points on a log-K3 surface, and no algorithm is
  even known to decide one `k`. The human dispatch's "42 is open" was seven years stale
  — 42 fell to Booker–Sutherland in Sept 2019 (we re-verified the 17-digit identity at
  residual 0) — but the *frontier* is real: {114, 390, 627, 633, 732, 921, 975} remain
  below 1000, and every one sits in the rigid classes `±3 mod 9` where solutions are
  forced to satisfy `x ≡ y ≡ z (mod 3)`. The forbidden classes `±4 mod 9` are the
  cubes-side forbidden-values law; the rigid classes hoard the open cases. Obstruction
  and stubbornness are mod-9 neighbors.
- **The cubic invariant of our own objects = no obstruction at all.** The repo's
  forbidden-values culture (H ∉ {7, 21}) primed us to hunt forbidden 3-cycle counts.
  The opposite is true (THM-462): the c₃ spectrum is gap-free `[0, M(n)]` at every n —
  Kendall–Babington Smith's 1940 identity makes c₃ a score-sequence functional, Landau
  plus a Lagrange four-square perturbation covers the top window, and a dominant vertex
  does the rest. **Impossibility begins strictly above the cubic layer.** The
  degree-3 channel (the FAST/contradiction channel of the trichotomy) has a complete
  integer scale; H's forbidden values are a higher-degree phenomenon.

## The dissolved question (the best kind of answer)

The backlog's last "honest bridge" handle asked: *why* does `|Aut| | H` for Paley
tournaments — is there a `ℚ(√−3)`/QR reason? LEM-003 dissolves it: Aut acts freely on
directed Hamiltonian paths of ANY digraph, because a path's arc set anchors at its
unique source. Zero Eisenstein content; the only QR fact is the SIZE `p(p−1)/2`. The
repo had recorded this fact four times (asserted, assumed, circularly "proved" —
MISTAKE-070 — and special-cased) without the one-paragraph proof. The honest boundary
is again where the content lives: the lemma FAILS for Hamiltonian cycles (RQ5 has both
its Ham cycles rotation-fixed) — freeness is a property of *paths having a first
vertex*, nothing more.

## Forward: the falsifiable edge

The cubic lens ends in a number we do not know yet. From the proven form
`R(p) = e(1 − C/p − …)` (with the proven negative that truncated cluster sums cannot do
better at finite p — the resurgence is violent, `δ₃/p⁴ ≈ +262`), this session predicts

> `R(31) = 2.59599 ± 0.00650`, i.e. `H(T_31) ≈ 1.988 × 10²⁵`, with `H ≡ 465 (mod 930)`,

and hands a validated compute design to whichever node runs it. `H(T_31)/465` — an odd
integer near `4.275 × 10²²` — is the sequence's next term after 1, 9, **1729**,
6857869865, 62293308207033. The number that started as Ramanujan's joke about a taxi
is now the third term of a sequence whose limit law is `e` and whose next term is a
prediction with an error bar. That is what "fleshed out" looks like.
