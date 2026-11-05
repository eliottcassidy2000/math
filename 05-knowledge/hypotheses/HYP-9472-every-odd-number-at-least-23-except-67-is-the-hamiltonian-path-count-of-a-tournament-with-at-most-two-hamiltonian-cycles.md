---
id: HYP-9472
title: "Every odd number N >= 23 other than 67 is the Hamiltonian-path count of a tournament with at most two Hamiltonian cycles; equivalently (THM-4877) of a uniquely Hamiltonian tournament or of a gluing of two rooted uniquely Hamiltonian tournaments. FINITE-EXACT: over orders 6-14 the uniquely Hamiltonian and two-cycle families together attain every odd number in [15, 2001) except 21 and 67 (21 is never a path count, THM-1370; 67 is attained only with three or four cycles at order 7); the uniquely Hamiltonian family alone misses 26 values below 2001, the two-cycle family 38; the two-cycle floor is f(n) + 8 a(n) with f Busch's strong floor and a = 1, 2, 3, 5 a(n-3). Pre-registered: orders 15-16 extend the window to about 4,000 with the same two exceptions."
status: >
  OPEN (klein-2026-11-05). FINITE-EXACT through order 14 from the complete generator of THM-4877 (data/cls_{n}.txt,
  scripts/hval.c, scripts/hc2_spectrum.py and its receipt in two_cycle_tournaments_gluing_20261104); the family is
  completely classified (THM-4877 (2)-(4)), so the conjecture is a numeration statement about pairs of rooted uniquely
  Hamiltonian tournaments, stronger than the Chan-Pak completeness (every odd except 7 and 21 is a path count, OPEN,
  MISTAKE-559) which it implies. No proof route beyond the composition laws (the path count of a gluing is a
  two-path-cover sum over the pieces; uniquely Hamiltonian pieces have Shelor's grammar and THM-4868's roundness) is
  known. The two-cycle floor formula is FINITE-EXACT at orders 6-14.
session: klein-2026-11-05 (codex/session-spectrum-20261104)
source: 05-knowledge/results/two_cycle_tournaments_gluing_20261104/README.md
related:
  - 01-canon/theorems/THM-4877-two-cycle-tournaments-complete-generator-census-to-16-gluing-of-two-uniquely-hamiltonian-pieces-and-the-count-law-as-a-proved-lower-bound.md
  - 01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md
  - 05-knowledge/hypotheses/HYP-9457-single-cycle-floor-band-and-lorentzian-join-crevasses.md
  - 05-knowledge/hypotheses/HYP-9459-few-cycle-tournaments-by-minimum-out-degree-f-2n-9-and-2f-2n-11-and-the-out-degree-1-lemma-for-one-and-two-cycles.md
  - 01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md
---

# HYP-9472 — the path spectrum of the tournaments with at most two Hamiltonian cycles

Let `𝒰` be the uniquely Hamiltonian tournaments and `𝒟` the tournaments with exactly two Hamiltonian cycles (`𝒟` = the vertex and arc gluings of rooted members of `𝒰`, THM-4877). Let `h(T)` be the number of Hamiltonian paths.

**Conjecture.** `{h(T) : T ∈ 𝒰 ∪ 𝒟} = {odd N ≥ 15} − {21, 67}`.

**Evidence (packet README §10).** Orders 6–14, all `121,393 + 361,595` classes down to order 6 through the complete generator: the union attains every odd number in `[15, 2001)` except 21 and 67; 21 is impossible (THM-1370) and 67 is realised only by strong tournaments with `hc ∈ {3, 4}` (order 7), so both exceptions are permanent (`f(15) = 1875`, `f_2(15) = f(15) + 8 · 125 > 67`). The uniquely Hamiltonian family alone has 26 permanent gaps below 1875; the two-cycle family fills all but 67 of them.

**Why it matters.** The Chan–Pak completeness (every odd number other than 7 and 21 is a Hamiltonian-path count; OPEN) would follow from this statement about a *classified* family: `h` of a gluing is an explicit two-path-cover sum over the two uniquely Hamiltonian pieces (PROOF_CONVERSE.md §8), and uniquely Hamiltonian tournaments have Shelor's insertion grammar (THM-4867) and THM-4868's run structure. A proof would be a numeration theorem for pairs of rooted uniquely Hamiltonian tournaments.

**Pre-registered tests.** (P1) Orders 15 and 16 (running at the time of writing) extend the covered window to about `2 · f_2(15)`–`4,000` with no new exception. (P2) The two-cycle floor stays `f(n) + 8a(n)` at orders 15–16 (2875, 4875). (P3) The number of permanent gaps of the uniquely Hamiltonian family grows (26 below 1875); the gaps of the union stay `{21, 67}`.
