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

As of 2026-07-21:

- **LRC(14) is OPEN.** The conjecture is known through 13 total runners; the
  repo's remaining work is structural, not the once-claimed uniform `q <= 25`
  shortcut, which is false. THM-2051 now closes the relation-dissociated
  branch after paying pairs exactly: every hypothetical counterexample has a
  genuine 3--5-term relation of height at most `2^20`. THM-2052 forces relation
  rank at least eleven: rank twelve is a finite maximal-minor box, while rank
  eleven reduces to finitely many two-anchor, one-projective-parameter stars.
  THM-2053 adds the exact transverse deck `D_N(m)` and a sufficient determinant
  gate whose strict failure carrier is an indexed union of 26 open tangent
  disks; THM-2055/2056 compress that fixed-basis
  determinant side into a signed-hull normal fan and a rational Kelvin/Farey
  certificate. THM-2057 closes two scaled AP one-tail planes: the omitted-12
  core by `12a`/`14a` clocks and an `84a|w` binding phase, and the full
  `{1,...,12}` core by `13a`/`14a` clocks and a `182a|w` deep-well phase. The
  later chain makes the residual sharper: THM-2058/2059 give phase packets and
  CRT joins; THM-2060/2064 isolate a dyadic seam; THM-2061/2062 fold it and add
  a hereditary wheel; THM-2065 makes circuit-free templates finite; THM-2066
  closes divisor-complete cores through `24`; THM-2068 compresses that bounded
  census; and THM-2072/2073 replace the impossible universal fixed bank by an
  adaptive depth-at-most-eight safe-child tower. LRC(14) remains open at its
  hereditarily primitive terminal, persistent circuits, and other atlas cells.
- **NC2 and GMC(2) are PROVED in repo canon** by THM-2022's lowest-balanced-face
  Frobenius argument. THM-2067 supplies its bare one-variable constant-term
  seed internally; full DvdK is only a stronger alternate. Lean checks the
  residue, face transport, contradiction, extractor, and height ingredients,
  while `DvdK1` and `HeightWitnessSupplier` remain explicit interfaces—the
  paper theorem is not yet fully formalized. THM-2070 explains why a cofinite
  support-return semigroup cannot replace coefficient noncancellation.
- **GMC is false from dimension 3 onward.** Dimension two is the sharp surviving
  case handled by THM-2022; dimensions one and two must not be conflated with
  the higher-dimensional counterexamples.
- **JC is false from dimension three; JC(2) and DC(2) remain open.** THM-1300
  gives the exact three-variable certificate. THM-2063 tames one-fiber-linear
  planar pairs, and THM-2071 forces a square leading coefficient with odd
  reduced degree at least three in the first quadratic-fiber survivor. Neither
  is a cover-degree classification or an NC2/GMC bridge (MISTAKE-236/237).
- **Tournament work has moved from isolated invariants to operations and
  decomposition:** order-join, strong cores, signed Rédei data, local
  subtournament censuses, and invariant-independence witnesses. Rank matches
  are not bridges without an explicit map, saturation index, and discriminant;
  small-prime Paley spectra are background, not an LRC periodic table.
- **Entropy and modular-form language is observable-relative, not a universal
  invariant.** MISTAKE-230--233 quarantine S217--S220, including the composite-
  `14` Paley/Frobenius analogy and the `f14`/modular-cusp bridge. No genus,
  period field, or Paley spectrum is presently an LRC obstruction without an
  explicit map preserving the loneliness predicate. HYP-8885's “cusp” frame is
  therefore a typed heuristic, not a universal difficulty theorem.
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
