# Agent Session Tools

This directory contains the repository's explicit session lifecycle tools.
There is deliberately no automatic Stop hook: close-out is an agent action,
not a background mutation of the checkout.

## Work on main in the checkout you were given (no worktrees)

Owner directive 2026-11-13: every session works directly on `main` in the
directory it was started in. Do not create worktrees or clones;
`agents/new_session_worktree.sh` now refuses. (Thirty full worktree checkouts
filled the volume on 2026-11-04; the 2026-11-13 cleanup removed ten more.)

The checkout is shared with other live sessions, so:

- stage and commit only the paths you wrote (`--paths` on the helpers below);
- never `git add -A`, `git stash`, `git reset --hard`, or `git checkout -- <path>`
  on files you did not write;
- sync with `git pull --rebase --autostash origin main` (autostash sets other
  sessions' uncommitted edits aside and restores them).

## Read current work

After fetching, generate a bounded packet and inspect the newest changes that
touch the current object rather than treating concurrent work as noise. Use
`--strict` when beginning from a supposedly clean/current worktree; it fails on
cached behindness, dirt, or a missing canonical route:

```bash
git fetch origin
python3 agents/start_session.py --strict --topic "<target statement or object>"
python3 agents/processor.py --status --peek --limit 8
```

The message processor uses the machine-local `.machine-id` and untracked
`agents/.read-log.json`. Do not commit either runtime file. `--check` and
`--status` show only the newest 12 unread messages per mailbox and 40 lines per
message by default; `--peek` leaves their read state unchanged. Use `--all` or
`--body-lines 0` only for an intentional archive read.

Validate the maintained startup surface after changing its routes or policies:

```bash
python3 agents/check_docs.py
```

Do not prepend generated session rows above `CURRENT DIGEST` or
`CURRENT-TRUTH WARNING` banners. Those byte-zero sentinels protect current truth
from chronological insertion. Put new rows below the maintained digest or in
the explicitly historical suffix. The startup router lists every `RESERVED`
file separately, whether empty or carrying a provisional candidate; content
and canon location do not upgrade its status.

Bounded current digests and `MISTAKES.md` intentionally do not use Git's union
merge driver. If they conflict, reconcile status and IDs semantically; never
paste an historical ledger back into a bounded startup file.

## Checkpoint useful partial results

Checkpoint after reserving scarce identifiers, completing a meaningful
computation, proving or refuting a claim, and before risky integration:

```bash
python3 agents/checkpoint_session.py \
  --message "[instance-id]: checkpoint - [brief state]" \
  --paths <files this session wrote>
```

The helper stages exactly `--paths`, commits only those paths (a pathspec
commit, so nobody else's staged entries are swept in), rebases with
`--autostash`, pushes to `origin/main`, and verifies. Without `--paths` it
refuses in the main checkout.

## Finish and push

Claude-style sessions can send a handoff, commit, rebase on concurrent mainline
work, push, and verify in one command:

```bash
python3 agents/finish_session.py \
  --to all \
  --subject "[instance-id]: [one-line summary]" \
  --body "Findings, evidence, remaining gaps, and next actions." \
  --commit-msg "[instance-id]: [one-line git summary]" \
  --paths <files this session wrote>
```

The closer adds the letter files it writes to `--paths`, pushes `HEAD` to
`origin/main`, and verifies that the session is no longer ahead of it. If a push or rebase fails, keep the local commit and
report the exact blocker; do not use destructive reset commands.

For a session that made no repository changes, still run and report:

```bash
git status --short --branch
git fetch origin
git rev-list --left-right --count origin/main...HEAD
```

Leftover worktrees from before 2026-11-13: `python3 agents/worktree_hygiene.py`
reports them; `--remove` deletes the clean merged ones.

## Why Stop hooks are absent

The former Stop hooks ran `check_session_closed.py` by a relative path (failing
outside the repository) and shared one tracked state file across sessions.
Close-out stays mandatory but explicit and stateless: inspect status, run the
checkpoint or finish helper, push, verify the upstream; session state stays
local and untracked.
