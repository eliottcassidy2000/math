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
  live prize is a lossless deck/clock/Farey/endpoint discharge for every
  remaining rank-eleven star cell; rank twelve is confined to a separate
  finite box that still needs exact decision. THM-2059 gives an exact
  arbitrary-clock CRT histogram join on one-tail cells; it does not force
  packet overlap.
  A failed gate means uncertified, not unsafe; THM-2058 is an empty reservation.
- **NC2 and GMC(2) are PROVED in repo canon** by THM-2022's lowest-balanced-face
  Frobenius argument. Its modular Lean spine is root-imported and sorry-free at
  the checked nodes, including descent, face construction, DvdK seed/reference
  interfaces, height/gap, Lucas, and whole-face Frobenius/noncancellation. A
  separate `GMC2NC2Capstone` module typechecks the conditional assembly skeleton
  with one `sorry` and is deliberately not root-imported. The remaining work is
  its concrete normalized-Wick-channel wiring and a formal proof of the cited
  DvdK input itself. HYP-8878 removes that citation on any lowest face with a
  unique minimum-mass balanced channel; its `98/116` prevalence is only one
  finite small-support census. MISTAKE-234 shows that support-return
  reachability does not control mixed-sign cancellation.
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
