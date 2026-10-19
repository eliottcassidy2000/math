# Math Research Workspace

A cumulative laboratory for proofs, refutations, computation, formalization,
and mathematical connections. The main prize is **LRC(14)**; the broader
portfolio includes tournaments, Gaussian moments, Mathieu–Zhao phenomena,
integer-sequence Dirichlet profiles, and adjacent problems.

## Start here

- [Startup brief](00-navigation/START-HERE.md) and [current frontier](00-navigation/CURRENT-FRONTIER.md)
- [Mandatory corrections](01-canon/ACTIVE-GUARDRAILS.md)
- [Research protocol](00-navigation/RESEARCH-PROTOCOL.md) and [meta-patterns](00-navigation/META-PATTERNS.md)
- [Core papers](05-knowledge/reference/CORE-PAPERS.md) and [problem ledger](00-navigation/PROBLEM-LEDGER.md)

```bash
python3 agents/start_session.py --topic "lonely runner"
```

Search historical ledgers on demand; do not reread them at every startup.

## Current orientation

These are routing summaries. For dated mechanisms, exact live obligations,
proofs, and corrections, follow the frontier and guardrails linked above.

- **LRC(14) is OPEN.** Strong reductions and exact finite certificates do not
  supply a full proof.
- **NC2 and GMC(2) are PROVED in repo canon** by THM-2022; GMC is false in every
  dimension at least three.
- **JC(2) and DC(2) remain OPEN.** Repo counterexamples concern JC from dimension
  three and Poisson/DC variants only in their stated higher ranks; they do not
  descend automatically to two variables.
- Tournament, knot, and integer-sequence claims are operation- and
  observable-relative. A bridge requires an explicit target-preserving map
  and the needed sidecar.
- Reciprocal integer sequences use the support Dirichlet profile plus collision
  tax, not a single reciprocal sum.
- Entropy, modular, Paley, and cusp language remains heuristic unless a typed
  map preserves the target predicate.

## Repository layers

| Path | Role |
|---|---|
| `00-navigation/` | Current routes, frontier, ledgers, historical maps |
| `01-canon/` | Definitions, scoped theorems, corrections |
| `02-court/` | Unresolved disputes and adversarial review |
| `03-artifacts/` | Drafts, reports, application artifacts |
| `04-computation/` | Reproducible searches, exact checks, Lean projects |
| `05-knowledge/` | Hypotheses, raw results, variables, references |
| `07-reflections/` | Idea provenance and structural synthesis |
| `agents/` | Session utilities and communication archives |

Truth precedence: corrections/repaired canon, proved canon, reproducible exact
computation, hypothesis, navigation synthesis, historical logs/reflections.
Legacy IDs collide: cite each theorem's ID and slug/path together.

## Contributing

Follow [AGENTS.md](AGENTS.md) and the research protocol: exact scope and
quantifiers, implication directions, mechanism, boundary, negative results,
frozen evidence, and connections to live concepts.

Maintain **Anchor / Niche / Wildcard**. LRC(14) is the default anchor; pursue
precise, promising ideas elsewhere too.
