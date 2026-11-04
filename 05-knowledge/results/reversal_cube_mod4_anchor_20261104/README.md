# One-child modulo-four separation of module-reversal cubes (THM-4864, HYP-9456)

kind-pasteur-2026-11-04-S9.
- **Status.** PROVED (THM-4864), given the cited reviewed inputs.
  - The pre-registered confirmations and their outcomes are in `PREREGISTRATION.md`.
  - The independent audit `audit/REVIEW.md` found it SOUND WITH CORRECTIONS. All corrections are applied
    (MISTAKE-703).
- **Question answered.** Session189 shows that, in a full-response collision, a transitive child of order ≥ 5 can
  only be **unchanged or wholly reversed**. This packet separates both cases modulo 4, at every order and in every
  fixed ambient.
  - **The query.** It is uniform on every other child and on the child minus one end u, and named only at u. Here u
    is the source or the sink.
  - **Widened form.** Transitivity is never used. Any child with a vertex that is a source or sink of both its source
    and target interiors is separated, provided m ≥ 2.
- **Also separated by coarse queries alone:** any quotient converse with a **singleton child**, for arbitrary
  same-size target interiors.

## Contents

- `PROOF.md`: the proof.
  - Contraction lemma, coefficient formulas, and the endpoint lemma.
  - **Weighted E1 (Proposition 4).** For every vertex v of a weighted tournament on at least 3 vertices,
    sum_{u≠v} [a_u c_v] delta + sum_{v→u} [a_u a_v c_v] delta = w(v) (mod 2).
  - Theorem A (singleton children) and Theorem B (transitive children, both targets).
  - The collision corollary.
- `PREREGISTRATION.md`: predictions P1–P5, pushed in 92aaef65f8 before any confirmatory run, with outcomes appended.
- `wmodel.py`: the weighted port-cycle parity polynomial g_Z, computed exactly over F2.
- `fcoef.cpp`: exact named response coefficients by Hamilton DP with Boolean Möbius inversion. It supports fixed
  exterior bits.
- `confirm.py`, `confirm.out`, `confirm.json`: the confirmatory runs P1-P5(b). All passed, with 0 failures.
- `p5c.py`, `p5c.out`, `p5c.json`: the P5(c) ambient controls. The guard in `confirm.py` was infeasible, so these
  ran separately with the same prediction. All passed. The disclosure is in `PREREGISTRATION.md`.
- `audit/`: the independent referee. It contains `REVIEW.md`, an engine of its own and drivers. Verdict: SOUND WITH
  CORRECTIONS, all applied.
- `residual_c3/`: EXPLORATORY probe of the residue. Paley-7 with C3 children, kept unchanged, is blind modulo 4 but
  exactly separated modulo 8.
  - 870 coefficient differences, minimum 2-adic valuation 2.
  - Exact run-word substitution DP, validated against direct counts.

## Reproduction

    g++ -O2 -std=c++17 fcoef.cpp -o fcoef.exe     # (or -o fcoef on POSIX; confirm.py picks the name)
    python confirm.py > confirm.out

## What the collision frontier looks like now

Combining S188, S189 and THM-4864: in a surviving original-module collision, every highest changed node satisfies
- its quotient Q is prime of order ≥ 7, with p(Q) = 0;
- none of its children is a singleton;
- none is a transitive tournament of order 2, 3 or ≥ 5;
- every transitive child of order 4 is reordered so that its two middle names become the ends.

The rest is HYP-9456. The Paley-7 C3 probe shows that its non-transitive part needs modulo 8.
