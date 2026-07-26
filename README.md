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

- [`START-HERE.md`](00-navigation/START-HERE.md) — startup brief and topic router.
- [`CURRENT-FRONTIER.md`](00-navigation/CURRENT-FRONTIER.md) — current status and live residuals.
- [`ACTIVE-GUARDRAILS.md`](01-canon/ACTIVE-GUARDRAILS.md) — corrections required before reuse.
- [`RESEARCH-PROTOCOL.md`](00-navigation/RESEARCH-PROTOCOL.md) and [`META-PATTERNS.md`](00-navigation/META-PATTERNS.md) — exploratory process and reusable moves.
- [`CORE-PAPERS.md`](05-knowledge/reference/CORE-PAPERS.md) and [`PROBLEM-LEDGER.md`](00-navigation/PROBLEM-LEDGER.md) — scoped external imports and the broader portfolio.

Agents should run:

```bash
python3 agents/start_session.py --topic "lonely runner"
```

The packet is intentionally bounded. Huge historical ledgers are searched when
needed, not read wholesale at every startup.

## Current orientation

These headline states are intentionally evergreen; dated mechanisms and exact
live debts belong to `CURRENT-FRONTIER.md`, which is authoritative before reuse.

- **LRC(14) is OPEN.** The repository has strong structural reductions and exact
  finite certificates, but no full proof.
- **NC2 and GMC(2) are PROVED in repo canon** by THM-2022; GMC is false in every
  dimension at least three.
- **JC(2) and DC(2) remain OPEN.** Repo counterexamples concern JC from dimension
  three and Poisson/DC variants only in their stated higher ranks; they do not
  descend automatically to two variables.
- Tournament, knot, and integer-sequence claims are operation- and
  observable-relative: an analogy or matching invariant is not a bridge without
  an explicit target-preserving map and the needed sidecar.
- Reciprocal integer sequences use the support Dirichlet profile plus collision
  tax, not a single reciprocal sum.
- Entropy, modular, Paley, and cusp language remains heuristic unless a typed map
  preserves the target predicate.

These are routing summaries, not proof citations. Follow exact slugs, statuses,
dependencies, and corrections in `CURRENT-FRONTIER.md` and `ACTIVE-GUARDRAILS.md`.

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
