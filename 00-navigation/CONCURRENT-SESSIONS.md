# Concurrent Session Playbook

**Current process — refreshed 2026-07-21.** `origin/main` is the shared live
research surface. Publish small coherent states often, but never absorb another
session's uncommitted files merely to checkpoint yours.

## Default: one worktree per session

From a clean main checkout:

```bash
cd "$(bash agents/new_session_worktree.sh short-topic)"
git status --short --branch
python3 agents/start_session.py --topic "<target statement or object>"
```

The helper fetches `origin/main`, creates a unique `codex/session-*` branch, and
sets `origin/main` as its upstream. The private working directory prevents one
session's stage/commit from sweeping another session's files. Platforms that
already supply an isolated clone/worktree need not create another.

If you must use a shared checkout, inspect every changed path before staging.
Use explicit `git add <paths>`; `git add -A` is allowed only after verifying
that every change belongs to the same session.

## Startup and sync

Sync only from a clean tree:

```bash
git fetch origin
git rebase origin/main
```

Then read the bounded startup packet and the exact canon it routes to. Do not
resolve a rebase by discarding unknown changes. If the tree is dirty, determine
ownership first; use a separate worktree rather than stashing another agent's
state.

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

In an isolated session worktree:

```bash
python3 agents/checkpoint_session.py \
  --message "[instance-id]: checkpoint - [brief mathematical state]"
```

A coherent checkpoint may be negative: a counterexample, failed implication,
smallest witness, repaired statement, or stopping certificate is useful state.

## Append-heavy files

`.gitattributes` uses the union merge driver for several logs and indexes. It
preserves both concurrent hunks but does not adjudicate truth, order entries,
or remove duplicate lines. After any automatic union merge:

- check that no heading or ID was duplicated;
- ensure correction/retraction order is intelligible;
- do not delete another agent's content merely for cosmetic uniformity; and
- repair current truth in the rolling frontier/guardrails rather than rewriting
  historical provenance.

The recent MISTAKE-226/227 collisions are the canonical warning: a clean rebase
can still create a semantically colliding ledger. Canonical 226 is the
augmentation/energy error, 227 the AP-chain saturation error, 233 the modular-
cusp bridge, and 235 the S102 weighted-fiber/AP-reduction correction.

## Computation and formalization

Reserve script and output names before long runs. Store the exact universe,
filters, positive/hostile controls, command, source/output pair, and hashes for
load-bearing results. A later pull may inspire a new run, but never silently
changes the interpretation of a frozen run.

For Lean, checkpoint both the module and the root-import wiring. A green
standalone file is not yet a project-level formalization claim.

## Close-out

In an isolated worktree, Claude-style sessions may use:

```bash
python3 agents/finish_session.py \
  --to all \
  --subject "[instance-id]: [one-line result]" \
  --body "Status, mechanism, evidence, corrections, and exact next obligations." \
  --commit-msg "[instance-id]: [one-line git summary]"
```

Otherwise: inspect status, stage explicit paths, commit, rebase/fetch, push the
intended target, and verify `HEAD` is not ahead of it. If a push fails, keep the
local commit and report the exact blocker. Never use a destructive reset to
make a concurrent session look clean.

See [`../agents/README.md`](../agents/README.md) for tool details and failure
recovery. There is deliberately no global stateful Stop hook.
