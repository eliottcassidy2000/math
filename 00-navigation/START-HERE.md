# Start Here

**Rolling brief — refreshed 2026-07-21; update whenever a headline status
changes.** This page is a router, not a proof source. The startup packet prints
the exact current commit. Status labels and proof links live in
[`CURRENT-FRONTIER.md`](CURRENT-FRONTIER.md).

## The five-minute orientation

1. **Main prize:** the Lonely Runner Conjecture for 14 total runners (13
   nonzero relative speeds) is **OPEN**. LRC is known through 13 total runners.
2. **Do not chase the old shortcuts:** uniform good period `q <= 25` is false,
   and uniform emptiness of the twelve-speed sporadic tight branch remains open.
3. **Major new closure:** NC2, hence unrestricted GMC(2), is **PROVED in repo
   canon** by THM-2022. Its full proof is not yet formalized in Lean.
4. **Truth discipline:** correction/repaired canon outranks canon, which
   outranks exact computation, hypothesis, synthesis, and historical prose.
5. **Research posture:** use an Anchor / Niche / Wildcard portfolio. Recover
   prior work, explain mechanisms and failures, generate alternate objects and
   quotients, and test connections rather than collecting slogans.

Run the bounded packet before reading large files:

```bash
python3 agents/start_session.py --topic "<target statement or object>"
```

Then read the relevant frontier section, the active guardrails, and the exact
theorem or evidence it links. `SESSION-LOG.md`, `TANGENTS.md`, reflections, old
frontier snapshots, the full mistakes ledger, and the giant hypothesis index
are searchable history—not mandatory warm-up.

## Twelve high-signal facts to reuse, not re-derive

1. **LRC(14) is the first unresolved runner count in the present literature.**
   Sungkawichai–Trakulthongchai prove the cases with 10, 11, and 12 nonzero
   speeds; see the scoped import in
   [`CORE-PAPERS.md`](../05-knowledge/reference/CORE-PAPERS.md).
2. **The `q <= 25` period claim is refuted.** THM-762/764 exhibit
   `26*{1,...,12} union {339}`, whose first good period is `27`.
3. **Good-period existence is a maximum statement, not an average or count.**
   The tight AP at its resonant ruler defeats the tempting mean arguments; see
   MISTAKE-127/129/130.
4. **The twelve-speed tight/sporadic locus is finite but not classified.**
   THM-763 gives a huge finite height; THM-1171 handles the AP locus; THM-1284
   eliminates the first-gap single-far stratum. Non-AP and deeper multi-defect
   branches remain.
5. **Exact LRC(14) computation is sound through height 55.** THM-1290, rerun
   after MISTAKE-194, also empties `(1/14,3/41)` through height 64.
6. **The 13-speed floor is isolated from above.** THM-1289 imports the
   Giri–Kravitz one-sided accumulation theorem. Its gap is ineffective; it does
   not make the whole first window finite without an additional conjecture.
7. **NC2 and GMC(2) are proved.** THM-2022 uses algebraic torus descent, a
   lowest balanced face, the Duistermaat–van der Kallen constant-term seed,
   and good-prime Kummer/Lucas/Frobenius amplification of the *whole* face.
8. **Whole-packet Frobenius is transferable but not an LRC proof.** THM-2041
   preserves exact-order/parity/conductor packets. LRC still needs a nonzero
   safe seed and a pointwise exit.
9. **Tournament structure is operational.** Order-join makes Hamiltonian-path
   count multiplicative and triangle count additive (THM-1862); signed Rédei
   data is join-multiplicative (THM-1936); zeta lives on the strong core
   (THM-1926). The known invariant lattice is exact only through `n <= 6`.
10. **A reciprocal integer sequence has two profiles.** THM-2000/2005 separate
    support from indexed multiplicity by a collision tax; the reusable object
    is the support Dirichlet profile, with Abel–Stieltjes/log-block Dini and the
    full Bertrand boundary at `z=1`.
11. **The rank-two Poisson, Dixmier, and planar Jacobian scopes now separate
    sharply.** THM-2044 proves the two-pair Poisson conjecture false by an
    explicit symplectic suspension. DC(2) and planar JC remain open; THM-2045
    proves only that the factorized family `R=x(a-b*x^r*q^s)` has no planar
    polynomial mate.
12. **Complete period-14 coordinates can still be globally blind.** THM-2043
    proves parity-Hasse completeness for reduced period-14 functions, then
    gives an infinite AP-alias family with the same local packet and a strict
    `17/41` exit. Exact owner height or adaptive resolved phase is indispensable.
    [THM-2047](../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md)
    now proves the labelled phase-height carrier and pair-sum recovery;
    deletion/localization forcing the AP core remains open.

## Where to go by topic

| Topic | Read first | Then use |
|---|---|---|
| LRC(14) proof | [`CURRENT-FRONTIER.md#lrc14`](CURRENT-FRONTIER.md#lrc14) | [`LRC14-PROOF-MAP.md`](LRC14-PROOF-MAP.md), [`HYP-8800`](../05-knowledge/hypotheses/HYP-8800-lrc14-face-carry-frobenius-transfer.md) |
| LRC literature | [`CORE-PAPERS.md#lonely-runner-conjecture`](../05-knowledge/reference/CORE-PAPERS.md#lonely-runner-conjecture) | [`LRC-TECHNIQUE-INDEX.md`](LRC-TECHNIQUE-INDEX.md), source paper |
| NC2 / GMC | [`CURRENT-FRONTIER.md#nc2-and-gaussian-moments`](CURRENT-FRONTIER.md#nc2-and-gaussian-moments) | THM-2022, THM-2040/2041, active guardrails |
| Jacobian / Dixmier / Poisson | [`CURRENT-FRONTIER.md#other-active-portfolio`](CURRENT-FRONTIER.md#other-active-portfolio) | [`CORE-PAPERS.md#jacobian-dixmier-and-poisson`](../05-knowledge/reference/CORE-PAPERS.md#jacobian-dixmier-and-poisson), THM-2044/2045 |
| Tournaments | [`CURRENT-FRONTIER.md#tournaments`](CURRENT-FRONTIER.md#tournaments) | [`TOURNAMENT-INVARIANT-ZOO-ATLAS-klein-S399.md`](TOURNAMENT-INVARIANT-ZOO-ATLAS-klein-S399.md), [`THE-ZOO.md`](THE-ZOO.md) with its freshness banner |
| Integer sequences | [`CURRENT-FRONTIER.md#integer-sequences`](CURRENT-FRONTIER.md#integer-sequences) | THM-2000/2005 and their frozen outputs |
| Broad portfolio | [`PROBLEM-LEDGER.md`](PROBLEM-LEDGER.md) | [`PROBLEM-ATLAS-2026-07-20.md`](PROBLEM-ATLAS-2026-07-20.md) only as historical provenance |
| Formalization | relevant theorem's `formalization` field | `00-navigation/*FORMALIZATION*`, then a clean build/axiom audit |
| Research method | [`RESEARCH-PROTOCOL.md`](RESEARCH-PROTOCOL.md) | [`META-PATTERNS.md`](META-PATTERNS.md), relevant reflections |
| Agent process | [`../agents/README.md`](../agents/README.md) | [`CONCURRENT-SESSIONS.md`](CONCURRENT-SESSIONS.md) |

## Before promoting a result

- Search the exact statement, not only the method vocabulary.
- Read [`ACTIVE-GUARDRAILS.md`](../01-canon/ACTIVE-GUARDRAILS.md) and recent
  matching `MISTAKE-*` entries.
- State scope, quantifiers, logical direction, mechanism, boundary, dependency,
  and what the result does **not** prove.
- For a connection, name the map, preserved predicate, information loss,
  sidecar, and decisive test.
- For computation, freeze universe, filters, controls, command, output, and
  code/output hashes where practical.
- For identifiers, check filename, frontmatter, indexes, and remote history;
  cite ID plus slug because legacy collisions exist.

## Maintaining this brief

Only put session-independent facts here. When a headline changes, update this
page, `CURRENT-FRONTIER.md`, and any stale startup banner in the same commit.
Preserve the old argument in its historical file; make supersession visible
rather than silently rewriting provenance. The startup script prints current
commits, so this page should not duplicate an inevitably stale commit diary.
