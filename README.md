# Math Research Workspace

This repository is a cumulative, multi-agent laboratory for proving,
refuting, computing, formalizing, and connecting mathematics. Its current main
prize is the **14-runner case of the Lonely Runner Conjecture**, but the corpus
also develops tournament invariants, Gaussian moments and Mathieu–Zhao
phenomena, integer-sequence Dirichlet profiles, formal proofs, and a wide
portfolio of adjacent problems.

The organizing principle is simple: current truth should be quick to find;
detailed evidence should be reproducible; and failed ideas should leave behind
more structure than they consumed.

## Start here

- [`00-navigation/START-HERE.md`](00-navigation/START-HERE.md) — the short
  startup brief and topic router.
- [`00-navigation/CURRENT-FRONTIER.md`](00-navigation/CURRENT-FRONTIER.md) —
  rolling mathematical status and exact live residuals.
- [`01-canon/ACTIVE-GUARDRAILS.md`](01-canon/ACTIVE-GUARDRAILS.md) — current
  corrections agents must know before reusing attractive claims.
- [`00-navigation/RESEARCH-PROTOCOL.md`](00-navigation/RESEARCH-PROTOCOL.md) —
  how exploratory sessions recover history, generate perspectives, test
  connections, and explain success or failure.
- [`00-navigation/META-PATTERNS.md`](00-navigation/META-PATTERNS.md) — reusable
  research moves with triggers and counterindications.
- [`05-knowledge/reference/CORE-PAPERS.md`](05-knowledge/reference/CORE-PAPERS.md)
  — exact external imports and what each paper does *not* prove.
- [`00-navigation/PROBLEM-LEDGER.md`](00-navigation/PROBLEM-LEDGER.md) — the
  broader problem portfolio.

Agents should run:

```bash
python3 agents/start_session.py --topic "lonely runner"
```

The packet is intentionally bounded. Huge historical ledgers are searched when
needed, not read wholesale at every startup.

## Current orientation

As of 2026-07-21:

- **LRC(14) is OPEN.** The conjecture is known through 13 total runners; the
  repo's remaining work is structural, not the once-claimed uniform `q <= 25`
  shortcut, which is false. THM-2051 now closes the relation-dissociated
  branch after paying pairs exactly: every hypothetical counterexample has a
  genuine 3--5-term relation of height at most `2^21`. THM-2052 already forces
  eleven independent bounded support-at-most-three relations and a rational
  plane; pointed
  rank-or-Euler transport is the live prize.
- **NC2 and GMC(2) are PROVED in repo canon** by THM-2022's lowest-balanced-face
  Frobenius argument. The Lean work is partial and sorry-free at its completed
  kernel nodes—including complete face-sum Frobenius/noncancellation—but
  algebraic descent, DvdK/Kummer survival and layer assembly, root wiring, and
  the final `nc2 : NC2` theorem remain.
- **GMC is false from dimension 3 onward.** Dimension two is the sharp surviving
  case handled by THM-2022; dimensions one and two must not be conflated with
  the higher-dimensional counterexamples.
- **The two-pair Poisson conjecture is false, while DC(2) and planar JC remain
  open.** THM-2044 is an explicit symplectic suspension; THM-2045 obstructs
  only its factorized first-coordinate family from de-stabilizing to the plane.
  THM-2049 makes the local correction complex acyclic, leaving polynomial
  termination and the coupled `D` column rather than a grade-six obstruction.
- **Tournament work has moved from isolated invariants to operations and
  decomposition:** order-join, strong cores, signed Rédei data, local
  subtournament censuses, and invariant-independence witnesses.
- **For reciprocal integer sequences, the underlying object is the support
  Dirichlet profile**, not one scalar reciprocal sum; repeated values require an
  explicit collision tax.

These are summaries, not proof citations. Follow the links in
`CURRENT-FRONTIER.md` before relying on them.

## Repository layers

| Path | Role |
|---|---|
| `00-navigation/` | Current routes, frontier, ledgers, and historical maps |
| `01-canon/` | Definitions, proved/cited theorems, and correction ledger |
| `02-court/` | Genuine unresolved disputes and adversarial review |
| `03-artifacts/` | Drafts, reports, and application artifacts |
| `04-computation/` | Reproducible searches, exact checks, and Lean projects |
| `05-knowledge/` | Hypotheses, raw results, variables, and references |
| `07-reflections/` | Idea provenance and cross-domain structural synthesis |
| `agents/` | Concurrent-session utilities and communication archives |

Truth precedence is correction/repaired canon, proved canon, reproducible exact
computation, hypothesis, navigation synthesis, then historical log/reflection.
A theorem identifier is not globally unique in legacy material, so cite its ID
and slug or file path together.

## Contributing

Read [`AGENTS.md`](AGENTS.md). Mathematical contributions should state exact
scope and quantifiers, distinguish necessary/sufficient/iff directions, explain
the mechanism and boundary, preserve negative results, freeze computational
evidence, and connect new objects to the other live concepts in the session.

The preferred research portfolio is **Anchor / Niche / Wildcard**. LRC(14) is
the default anchor; it is never a reason to suppress a precise and promising
idea elsewhere.
