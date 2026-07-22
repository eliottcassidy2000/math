# Agent Operating Contract

This is a persistent, concurrent mathematical research workspace. `AGENTS.md`
is the small universal policy; current mathematics lives in the documents it
routes to. Never use a chronological log as the current truth source.

## Bounded startup

1. Start from a clean, current worktree. Prefer an isolated worktree made by
   `bash agents/new_session_worktree.sh <tag>`. If a shared checkout is dirty,
   identify and preserve every unrelated change before syncing.
2. Run `python3 agents/start_session.py --topic "<your topic>"`. Then read:
   - `00-navigation/START-HERE.md`;
   - the relevant section of `00-navigation/CURRENT-FRONTIER.md`;
   - `01-canon/ACTIVE-GUARDRAILS.md`; and
   - the topic routes printed by the startup packet.
3. Read `00-navigation/RESEARCH-PROTOCOL.md` before a mathematical session and
   `05-knowledge/reference/CORE-PAPERS.md` before making literature claims.
4. Check targeted agent/human messages with
   `python3 agents/processor.py --status --peek --limit 8`. Do not scan the full
   session log, mistakes ledger, backlog, or hypothesis index unless needed.

The intended startup is minutes and hundreds of lines, not a corpus reread.

## Truth and status discipline

Use this precedence order:

1. explicit correction/retraction and its repaired theorem;
2. current proved canon;
3. reproducible exact computation;
4. named hypothesis;
5. current navigation synthesis;
6. historical session log, reflection, message, or draft.

Reflections and logs preserve idea provenance, not truth. Before using a claim,
search its exact statement, constants, quantifiers, synonyms, theorem ID, and
recent `MISTAKE-*` entries. Cite theorem ID **plus slug/file path**: legacy ID
collisions exist. Mark scope as `PROVED`, `CITED`, `FINITE-EXACT`, `VERIFIED`,
`CONDITIONAL`, `OPEN`, `RESERVED`, `REFUTED`, or `SUPERSEDED`; never blend them.
`RESERVED` is never a proved result or proved dependency. A reserved file may
be an empty namespace stub or an explicitly provisional proof candidate under
audit; only an audited status promotion moves it into the proof graph.

Repair demonstrated errors promptly, preserve the correction lineage, and log
the mechanism in `01-canon/MISTAKES.md`. Use court cases for genuine unresolved
disputes, not as a reason to leave a known false statement live.

## Mathematical-session directive

The main near-term prize is **LRC(14)**, which is open. It is an anchor, not a
tunnel. Every self-directed session keeps an **Anchor / Niche / Wildcard**
portfolio: the assigned or highest-value target, an orthogonal underexplored
thread, and a curiosity lane. A niche or wildcard may overtake the anchor when
it yields a precise object, obstruction, experiment, or theorem.

Always:

- begin with an inheritance pass: name the closest proved mechanism, canonical
  hostile example, corrected near miss, and least-used relevant sidecar;
- recover and connect prior work before deriving anew;
- prefer underexplored operations, duals, scales, boundary cases, and discarded
  objects over another pass through a saturated route;
- determine **why** a claim is true or false, not merely whether it survives;
- keep a board of 3–7 live concepts and compare every new idea with each one;
- generate views from `objects x representations x invariants x operations x
  symmetries/quotients x scales`;
- specify a claimed connection's source, target, map, preserved predicate,
  destroyed information, needed sidecar, and cheapest decisive test;
- give every small mathematical compulsion a cheap hostile probe, then pursue a
  positive signal until it produces structure or a recorded stopping reason;
- revisit the concept board after every meaningful pull or computation and ask
  how each new object changes every other live lane;
- promote successful research moves to `00-navigation/META-PATTERNS.md` only
  with triggers, counterindications, and evidence from distinct threads.

For a true result, record the mechanism, quantifiers, equality/failure boundary,
dependencies, and possible generalization. For a false result, record the
minimal witness, first failed implication, strongest survivor, repaired form,
missing coordinate, and new question.

## Validity gate

Before canonizing work, audit types; necessary/sufficient/iff directions;
symmetries and orbit representatives; quotient losses; canonical adversarial
examples; and the actual consequence rather than an intermediate statistic.
Computations need an explicit universe, inherited filters, positive and hostile
controls, a reproduction command, and preferably an independent path. Lean
claims need a build, axiom audit, satisfiable hypotheses, and root-import reach.

Tournament Analysis is encouraged only when an intrinsic binary relation exists.
Declare vertices, pairwise observable, orientation gauge, ties, preserved target,
lost data, and sidecar. Challenge runner/arc vertices: gaps, sections, boundaries,
events, residues, Fourier modes, circuits, and proof obligations may be better.
Never force ties into a cosmetic tournament or treat a sufficient tournament
certificate as an equivalence.

## Concurrent Git protocol

Treat `origin/main` as the live shared surface. Pull/rebase only from a clean
tree. Read incoming commits as mathematical signal and integrate real
connections. Reserve scarce IDs and filenames with honest stubs only after
checking filename, YAML ID, indexes, and remote history. An empty theorem stub
must say `RESERVED / UNPROVED EMPTY STUB` in its status and body, have no proved
dependencies, and be routed separately from proved canon.

Never prepend session prose above a maintained warning, current digest, or
historical-boundary banner. Insert immediately below that banner or in the
historical suffix; `agents/check_docs.py` treats byte-zero displacement as a
truth-surface failure.

Push small coherent checkpoints after reservations, meaningful computations,
proof/refutation milestones, and before risky rebases. In an isolated worktree:

```bash
python3 agents/checkpoint_session.py --message "[instance-id]: checkpoint - [state]"
```

In a shared dirty checkout, stage only your explicit paths; never sweep another
session's work with `git add -A`. See `00-navigation/CONCURRENT-SESSIONS.md` and
`agents/README.md`.

## Mandatory close-out

Before the final response, inspect `git status --short --branch`. If the session
made intentional changes, stage only those paths, commit, and push to the
upstream (or `origin/main` for a session worktree). In every session, fetch and
verify `HEAD` is not ahead of its push target; a no-change session needs no
empty commit. Use `agents/finish_session.py` only in an isolated worktree. If a
required push fails, the session is not complete: leave a local commit and
report the exact blocker.
