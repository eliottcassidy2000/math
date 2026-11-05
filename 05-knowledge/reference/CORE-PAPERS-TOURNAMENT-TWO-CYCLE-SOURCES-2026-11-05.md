# Thomassen 1980 / Busch 2006 source intake: Hamiltonian cycles in 2-connected tournaments

**Status:** CITED PUBLISHED SOURCES (klein-2026-11-05). This sidecar keeps the
bounded-startup catalog compact; it imports one classical lower bound and
records exactly how it is used and how it was read.

## Primary sources

- C. Thomassen, "On the number of Hamiltonian cycles in tournaments,"
  *Discrete Mathematics* **31** (1980), 315--323,
  [DOI 10.1016/0012-365X(80)90142-9](https://doi.org/10.1016/0012-365X(80)90142-9).
- A. H. Busch, "A note on the number of Hamiltonian paths in strong
  tournaments," *Electronic Journal of Combinatorics* **13** (2006), #N3,
  [article](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3).
- J. W. Moon, "The minimum number of spanning paths in a strong tournament,"
  *Publ. Math. Debrecen* **19** (1972), 101--104 (the path bound Busch makes exact).

## Exact imported role

Thomassen proves that every 2-connected tournament of order `n` has at least
`α^{c(n)}` Hamiltonian cycles with `α = ⁴√6 ≈ 1.565` and `c(n)` an explicit
linear function of `n`; Busch's Corollary 2 replaces `α` by `β = ∛5 ≈ 1.710`
with the same exponent, as a consequence of his exact value of Moon's minimum
number `hp(n)` of Hamiltonian paths in a strong tournament
(`hp(n) ≥ β^{n−1}` for `n ≡ 1 (mod 3)`, `3β^{n−3}` for `n ≡ 0`, `9β^{n−5}` for `n ≡ 2`).
The survey of Bang-Jensen and Gutin states the result qualitatively: a
constant `α > 1` with at least `α^n` cycles for every 2-strong tournament, and
Thomassen notes `α < 2` is necessary.

**Reading caveat.** The journal typesetting of the exponent in Busch's note
(`β^{(n/3 ... 2 − 1)}`) extracted as "n/32 − 1"; the readings consistent with
the regular tournament `RT₅` (2-connected, two Hamiltonian cycles) are
`n/3 − 1` and `(n − 3)/2`-type exponents, which give at least three cycles from
order 7; the weakest conceivable reading `n/6 − 1` gives three cycles from
order 14. The repository uses only: **every 2-connected tournament of order
`n ≥ 14` has at least three Hamiltonian cycles**, and the two-cycle census of
`two_cycle_tournaments_gluing_20261104` covers orders 6--16 directly.

Consumers: [THM-4877](../../01-canon/theorems/THM-4877-two-cycle-tournaments-complete-generator-census-to-16-gluing-of-two-uniquely-hamiltonian-pieces-and-the-count-law-as-a-proved-lower-bound.md)
(Lemma B of `PROOF_CONVERSE.md`: a tournament with exactly two Hamiltonian
cycles of order `≥ 6` has a cut vertex).

## Scope firewall

Neither source proves anything about the structure of tournaments with
exactly two Hamiltonian cycles, the gluing classification, the count law
`(n − 2)F(2n − 9) − L(2n − 10)`, or the Hamiltonian-path spectrum; the
minimum number of Hamiltonian cycles in 2-connected tournaments is open
(Kühn--Osthus survey, §4.2).
