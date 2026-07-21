# Agent Session Tools

This directory contains the repository's explicit session lifecycle tools.
There is deliberately no automatic Stop hook: close-out is an agent action,
not a background mutation of the checkout.

## Start in an isolated worktree

Concurrent sessions should not share uncommitted files. From the main checkout,
create a session worktree and enter it:

```bash
cd "$(bash agents/new_session_worktree.sh short-topic)"
git status --short --branch
```

The script fetches `origin`, creates a unique local branch named
`codex/session-short-topic` (adding a suffix on collision), and makes
`origin/main` its upstream. It also copies the untracked `.machine-id` when one
is present. A failed fetch or a non-worktree path collision stops with an error
instead of silently starting from stale or unrelated state.

If a platform already gives the session an isolated clone or worktree, use that
checkout directly. In a shared checkout, inspect `git status` before staging:
uncommitted files may belong to another live session.

## Read current work

After fetching, inspect the newest changes that touch the current object rather
than treating concurrent work as noise:

```bash
git log --oneline -10 --name-only
python3 agents/processor.py --check
```

The message processor uses the machine-local `.machine-id` and untracked
`agents/.read-log.json`. Do not commit either runtime file.

## Checkpoint useful partial results

Checkpoint after reserving scarce identifiers, completing a meaningful
computation, proving or refuting a claim, and before risky integration:

```bash
python3 agents/checkpoint_session.py \
  --message "[instance-id]: checkpoint - [brief state]"
```

The checkpoint helper stages the whole current worktree. This is safe in the
isolated session worktree; in a shared checkout, first verify every changed file
is intentional.

## Finish and push

Claude-style sessions can send a handoff, commit, rebase on concurrent mainline
work, push, and verify in one command:

```bash
python3 agents/finish_session.py \
  --to all \
  --subject "[instance-id]: [one-line summary]" \
  --body "Findings, evidence, remaining gaps, and next actions." \
  --commit-msg "[instance-id]: [one-line git summary]"
```

Because branches made by `new_session_worktree.sh` track `origin/main`, the
closer pushes `HEAD` to main and can verify that the session is no longer ahead
of the shared upstream. If a push or rebase fails, keep the local commit and
report the exact blocker; do not use destructive reset commands.

For a session that made no repository changes, still run and report:

```bash
git status --short --branch
git fetch origin
git rev-list --left-right --count origin/main...HEAD
```

Remove clean temporary worktrees later with:

```bash
bash agents/cleanup_session_worktrees.sh
```

## Why Stop hooks are absent

The former Codex and Claude Stop hooks invoked
`python3 agents/check_session_closed.py` through a relative path and failed in
sandboxed or non-repository working directories. Their shared tracked
`agents/.session-state.json` also mixed independent sessions, dirtied clean
trees, and could mistake another agent's message for the current handoff.

Close-out therefore remains mandatory policy, but enforcement is explicit and
stateless: inspect status, run the checkpoint or finish helper, push, and verify
the upstream. Runtime session state must remain local and untracked.
