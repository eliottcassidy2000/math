---
id: THM-1845
title: "THE LARGEST TRANSITIVE SUBTOURNAMENT IS SANDWICHED n − c₃ ≤ β ≤ smax + 1, both PROVED — and this is the first fruit of a TxGraffiti / Written-on-the-Wall style AUTOMATED CONJECTURE GENERATOR for tournament invariants (the WOWII front-end klein-S395 named as the one ingredient the repo lacked; owner: leverage the WOWII-103 refutation). Let β = the size of the largest transitive (acyclic) subtournament, c₃ = number of 3-cycles = tr(A³)/3, smax = maximum score. (1) LOWER, β ≥ n − c₃ (min feedback vertex set ≤ c₃): a tournament is acyclic iff it has no 3-cycle, so deleting one vertex from each 3-cycle (≤ c₃ vertices) leaves a 3-cycle-free = transitive subtournament of size ≥ n − c₃. TIGHT exactly on the near-transitive family — c₃ = 0 (transitive, β = n) and c₃ = 1 (THM-1830's transitive-skeleton-+-one-3-cycle, β = n − 1) — so the WOWII 'odd-cycle core + pendants' witness template IS my unstable-non-transitive family. (2) UPPER, β ≤ smax + 1: the top vertex of a transitive subtournament of size β has out-degree β − 1 inside it, hence ≥ β − 1 overall, so smax ≥ β − 1. (3) THE GENERATOR, run over all 33864 tournaments n = 3..6, reproduces known structure (H and the arborescence count arb₀ are INCOMPARABLE — mac-mini THM-1580 — both directions auto-refuted) and surfaces further candidates holding to n = 7: c₃ ≤ H, H ≤ 2^{n−2}·c₃ + 1. (4) A WOWII-STYLE OFF-BY-ONE REFUTATION found automatically: srange ≤ β holds for all n ≤ 6 but FAILS at n = 7 (witness c₃ = 4, srange = 6, β = 5) — the same n = 7 phase-transition wall as THM-1825/1830, now surfaced by the conjecture engine rather than by hand. So the repo now HAS the WOWII loop: generate inequalities on the invariant zoo → the tight survivors are candidate theorems (two proved here) → the failures come with an explicit tuned witness, and the witness template is the 3-cycle atom"
status: >
  (1) PROVED (feedback-vertex-set / tournament-acyclic-iff-no-3-cycle, elementary and
  complete).  TIGHTNESS on c₃ = 0, 1 verified (THM-1830 family).
  (2) PROVED (one line).
  (3) The generator is VERIFIED to reproduce THM-1580 incomparability; c₃ ≤ H and
  H ≤ 2^{n−2}c₃+1 are candidate inequalities holding on all n = 3..6 exhaustively and a
  150k-sample at n = 7 — NOT proved, offered as generator output for follow-up.
  (4) srange ≤ β VERIFIED to hold exhaustively n ≤ 6 and to fail at n = 7 with the stated
  exact witness.
  This is a TOOL + two small proved theorems + candidate inequalities.  The value is the
  loop, not any single inequality; β ≥ n − c₃ is likely folklore (it is the FVS bound) and is
  claimed as generator-surfaced-and-proved, not as new.
source: kind-pasteur-2026-07-21-S128c134 (owner: leverage the WOWII-103 counterexample ideas for repo problems / new ones)
depends_on:
  - THM-1830    # the 3-cycle-atom family = the tournament WOWII witness template; tightness of (1)
related: [THM-1580, THM-1780, THM-1805]
external:
  - "DeLaViña, Written on the Wall II (Graffiti.pc); google-deepmind/formal-conjectures PR #4482 (WOWII-103 refutation). klein-S395 reflection."
script: 04-computation/tournament_graffiti_kps_S128c134.py (+ .out)
---

# THM-1845 — the transitive-subtournament sandwich, from a tournament Graffiti engine

The owner asked to leverage the WOWII-103 refutation. klein-S395 named the missing ingredient:
the repo has the *search → exhaustive-check → Lean* discipline but not the **conjecture
generator** that Graffiti/WOWII puts in front. This builds that front end for tournament
invariants and reports its first fruit.

## The sandwich (both proved)

Let `β` = size of the largest transitive (acyclic) subtournament, `c₃` = number of 3-cycles,
`smax` = maximum score.

> **`n − c₃ ≤ β ≤ smax + 1`.**

**Lower (`β ≥ n − c₃`).** A tournament is acyclic iff it contains no 3-cycle (a tournament with
any directed cycle has a directed 3-cycle). Choose one vertex from each 3-cycle — a set `W` with
`|W| ≤ c₃`. Every 3-cycle of the original has a vertex in `W`, so the induced subtournament on
`V∖W` has no 3-cycle, hence is transitive, of size `≥ n − c₃`. ∎ (This is min-feedback-vertex-set
`≤ c₃`.)

**Upper (`β ≤ smax + 1`).** In a transitive subtournament of size `β`, the source vertex has
out-degree `β − 1` *within* it, so its out-degree in the whole tournament is `≥ β − 1`; thus
`smax ≥ β − 1`. ∎

**Tightness of the lower bound is exactly the WOWII witness template.** `β = n − c₃` holds at
`c₃ = 0` (transitive, `β = n`) and `c₃ = 1` (THM-1830's transitive skeleton + one 3-cycle atom,
`β = n − 1`). The 3-cycle atom is the tournament analog of WOWII-103's triangle, and the
transitive singletons are its "leaves."

## The generator, and what it found

`tournament_graffiti_kps_S128c134.py` computes a battery of iso-invariants (`c₃, H, β, dom,
kings, scc, smax, smin, srange, sumC2, arb₀`) over all `33 864` tournaments `n = 3..6` and
machine-generates tight linear inequalities between them.

- **It reproduces known structure.** `H` (Hamiltonian-path count) and `arb₀` (rooted
  arborescence count) are **incomparable** — both `H ≤ arb₀` and `H ≥ arb₀` are auto-refuted at
  `n = 3` — exactly mac-mini THM-1580.
- **Candidate inequalities** holding on all `n ≤ 6` and a `150k` sample at `n = 7` (offered, not
  proved): `c₃ ≤ H`; `H ≤ 2^{n−2}·c₃ + 1`.
- **A WOWII-style off-by-one refutation, found automatically.** `srange ≤ β` holds for every
  `n ≤ 6` but **fails at `n = 7`** — witness `c₃ = 4, srange = 6, β = 5` (`6 ≤ 5` false). This is
  the *same* `n = 7` wall as THM-1825/1830, this time surfaced by the engine rather than by hand.

## The point: the repo now has the WOWII loop

> **generate** inequalities on the invariant zoo → **the tight survivors are candidate
> theorems** (two proved here) → **the failures come with a tuned witness**, and the witness
> template is the **3-cycle atom** (THM-1830).

That is the exact shape of WOWII-103 (an easy invariant bounding a hard one, refuted by an
odd-cycle-core-plus-pendants graph) transported to tournaments. It is also, in miniature, what
the project has been doing by hand for a hundred sessions — "the pattern breaks at `n = 6/7`" —
now with a machine front end that proposes the patterns.

## Named next

- **Prove or refute `c₃ ≤ H` and `H ≤ 2^{n−2}c₃ + 1`** past `n = 7`; the former looks like it
  wants an injection from 3-cycles into Hamiltonian paths.
- **Run the full WOWII list on `G_n`, `E_n`** (klein-S395 started this) and on the tournament
  invariants at `n = 7, 8` — every off-by-one survivor is a candidate theorem, every failure a
  structured witness.
- **Formalize the sandwich in Lean** (`β`, `c₃`, `smax` are all `decide`-able on a fixed
  tournament) — the WOWII PR's `native_decide` discipline applies directly, and the repo's
  `TournamentH7` already has the harness.
