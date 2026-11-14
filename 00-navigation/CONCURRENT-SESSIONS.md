# Concurrent Session Playbook

**Current process — refreshed 2026-07-22.** `origin/main` is the shared live
research surface. Publish small coherent states often, but never absorb another
session's uncommitted files merely to checkpoint yours.

## Default: work on main in the checkout you were given

Owner directive 2026-11-13: no worktrees and no extra clones. Every session
works on `main` in the directory it was started in:

```bash
git pull --rebase --autostash origin main
python3 agents/start_session.py --topic "<target statement or object>"
```

The checkout is shared with other live sessions. Stage only the paths you
wrote (`git add <paths>`, then `git commit -m ... -- <paths>`); never
`git add -A`, `git stash`, `git reset --hard`, or `git checkout -- <path>` on
another session's files. Commit and push small pieces often so the shared tree
stays nearly clean.

## Startup and sync

Sync with autostash (other sessions' uncommitted edits are set aside and
restored):

```bash
git pull --rebase --autostash origin main
```

Then read the bounded startup packet and the exact canon it routes to. Do not
resolve a rebase by discarding unknown changes. If the tree is dirty, determine
ownership first and leave other agents' files alone; commit your own paths
and let `--autostash` carry the rest.

When a fetch/rebase lands new work, ask:

- Does it change this statement, invariant, runner family, scope, or namespace?
- Does it correct, subsume, independently verify, or refute the current route?
- Does it supply a new object, operation, sidecar, extremal, or decisive test?
- Should the connection enter the concept board, theorem, hypothesis, result,
  guardrail, frontier, or reflection?

This is the point of live mainline research: incoming work is potential
mathematics, not just Git noise.

## Claim scarce names honestly

Before reserving `THM-*`, `HYP-*`, `MISTAKE-*`, tangent IDs, script/output
names, or session territory:

1. search the filename namespace;
2. search frontmatter IDs and all indexes;
3. fetch and inspect current `origin/main`; and
4. cite an existing item by ID **plus slug/path** because legacy collisions
   remain.

An empty stub states what is reserved and why, has no dependencies, and says
`RESERVED / UNPROVED EMPTY STUB`. A provisional proof candidate may reserve a
theorem slug only when its status, missing audit, and candidate prerequisites
are explicit; it remains outside the proof graph. Broader speculation belongs
in hypotheses/tangents.

## Checkpoint rhythm

Checkpoint after a reservation, a meaningful computation, a proof/refutation,
a corrected hypothesis, a reusable connection, and before a risky integration.
During long computations, checkpoint the harness/input specification first;
fetch and inspect new work at natural yields without changing the inputs of an
already-running reproducibility job.

```bash
python3 agents/checkpoint_session.py \
  --message "[instance-id]: checkpoint - [brief mathematical state]" \
  --paths <files this session wrote>
```

A coherent checkpoint may be negative: a counterexample, failed implication,
smallest witness, repaired statement, or stopping certificate is useful state.

## Append-heavy files

`.gitattributes` uses the union merge driver only for genuinely append-only
archives. It is deliberately disabled for bounded `SESSION-LOG.md`, the current
hypothesis digest, and `MISTAKES.md`: those truth surfaces require semantic
reconciliation. After any automatic union merge elsewhere:

- check that no heading or ID was duplicated;
- ensure correction/retraction order is intelligible;
- do not delete another agent's content merely for cosmetic uniformity; and
- repair current truth in the rolling frontier/guardrails rather than rewriting
  historical provenance.

Recent correction-ID collisions are the canonical warning: a clean rebase can
still create a semantically colliding ledger. Canonical 230/231 are the inverse-
class and observable-entropy corrections; 237 is the HYP-8905 JC bridge; 238
is the empty-full-row dyadic transport error; and 239 is reversed guard
containment. Always fetch and search filename, heading, and remote history.

## Computation and formalization

Reserve script and output names before long runs. Store the exact universe,
filters, positive/hostile controls, command, source/output pair, and hashes for
load-bearing results. A later pull may inspire a new run, but never silently
changes the interpretation of a frozen run.

For Lean, checkpoint both the module and the root-import wiring. A green
standalone file is not yet a project-level formalization claim.

## Close-out

Claude-style sessions may use:

```bash
python3 agents/finish_session.py \
  --to all \
  --subject "[instance-id]: [one-line result]" \
  --body "Status, mechanism, evidence, corrections, and exact next obligations." \
  --commit-msg "[instance-id]: [one-line git summary]" \
  --paths <files this session wrote>
```

Otherwise: inspect status, stage explicit paths, commit, rebase/fetch, push the
intended target, and verify `HEAD` is not ahead of it. If a push fails, keep the
local commit and report the exact blocker. Never use a destructive reset to
make a concurrent session look clean.

See [`../agents/README.md`](../agents/README.md) for tool details and failure
recovery. There is deliberately no global stateful Stop hook.
