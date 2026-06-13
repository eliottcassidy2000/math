# Concurrent Session Playbook

This repository is a shared research instrument, not a queue. Multiple agents
may be working from `main` at the same time, and the correct response is to
publish small, coherent states often enough that everyone can see the living
frontier.

## Operating Principle

Push to `origin/main` whenever a useful unit of state exists. A useful unit can
be a proof attempt, a refutation, a hypothesis-number reservation, a result
file, a script checkpoint, a newly found connection, or a session-log update.
Do not wait for a perfect ending if the work has become informative.

Frequent pushes do three jobs:

1. They reserve scarce names such as `HYP-*`, `THM-*`, `T*`, result filenames,
   and session IDs before another agent independently uses them.
2. They turn long investigations into visible increments, reducing painful
   late-session rebases.
3. They let concurrent discoveries cross-pollinate while the ideas are still
   live.

## Checkpoint Rhythm

Make an early checkpoint after startup if you reserve any IDs or create any
stub files. During long sessions, checkpoint every 30-60 minutes, before a long
computation, after a meaningful computational result, after confirming or
refuting a hypothesis, and before any risky rebase or large edit.

Preferred command:

```bash
python3 agents/checkpoint_session.py \
  --message "[instance-id]: checkpoint - [brief state]"
```

The checkpoint helper stages intentional repo changes, commits if needed,
pushes the current branch, and verifies the branch is not ahead of upstream.
Use a direct `git add -A && git commit && git push` sequence if the helper is
not available.

## Claiming Numbers and Names

When you know a session will need a scarce namespace, claim it early with a
small commit. Examples:

- Add a one-paragraph placeholder hypothesis detail file plus an index row.
- Add a result filename to `05-knowledge/results/INDEX.md` before launching a
  long run.
- Add a theorem stub only when the certainty threshold is already met; otherwise
  use a hypothesis, tangent, or backlog lead.
- Add a session-log stub if a long investigation is underway and other agents
  need to know what territory is occupied.

The placeholder must state what is claimed, why, and what evidence is still
missing. Do not use a placeholder to smuggle speculation into canon.

## Rebase Noise Is Signal

When a rebase or fetch brings in new work, treat it as possible evidence, not
mere conflict noise.

Before continuing, quickly ask:

- Does the new work touch the same invariant, runner family, namespace, theorem,
  script, or proof route?
- Does it refute, rename, sharpen, or duplicate my current hypothesis?
- Does it create a new pairwise observable, gauge, tie path, result file, or
  application bridge that my session should compare against?
- Is there a tangent or backlog lead that should be updated because the two
  threads unexpectedly meet?

If the answer is yes, integrate the connection explicitly in the notes, script
output, hypothesis file, reflection, or session log. If the answer is no, keep
the new work intact and continue. Never delete another agent's fresh work just
to make your rebase feel cleaner.

## Mathematical Standard

Novelty pressure does not lower the burden of proof. Keep the certainty scale
intact:

- Canon only gets proved or proof-sketch-level results.
- Hypotheses record both positive and negative evidence.
- Computations save scripts and outputs together.
- Tournament Analysis remains the default computational lens when it fits:
  pairwise observable, switch/gauge, tie Hamiltonian path, and fingerprints.

The project wants large-scale connection search and rigorous mathematics at the
same time. A good session should leave behind either a stronger theorem, a
clearer failed route, a useful computational artifact, or a newly connected
piece of the map.

## Close-Out

End-of-session close-out is still mandatory. A session can checkpoint many
times and must still finish with the normal letter, commit, push, and
verification that the branch is no longer ahead of upstream.

## Infrastructure: working-tree isolation + conflict-free logs (added 2026-06-01)

Two mechanisms now absorb most concurrency pain so close-out no longer stalls on
manual conflicts:

### 1. Union-merge on shared logs (automatic — nothing to do)

`.gitattributes` marks the append-heavy coordination files (`SESSION-LOG.md`, the
`INDEX.md` files, `MISTAKES.md`, `TANGENTS.md`, `INVESTIGATION-BACKLOG.md`,
`OPEN-QUESTIONS.md`, this file) with the built-in `merge=union` driver. When two
sessions both append, a rebase keeps **both** sides' lines automatically — no
conflict markers, no manual fixup. `finish_session.py` now loops the
fetch/rebase/push retry and, on any conflict the union driver cannot resolve,
**aborts the rebase to leave a clean tree** (your commit stays intact) instead of
stranding the repo mid-rebase. You normally never see a conflict again.

Caveat: union merge keeps both hunks but not a guaranteed order, and can leave a
duplicate line if two sessions add the *same* line. Glance at the top of
`SESSION-LOG.md` after close-out; cosmetic dedupe is fine, never delete another
agent's entry.

### 2. Per-session worktrees (opt-in — recommended for remote-control sessions)

Nomad-launched sessions already isolate (they clone into `/tmp`). Remote-control
sessions share the one `~/math` checkout, so one session's `git add -A` can sweep
another's **uncommitted** files. To get the same isolation, start your session in
its own git worktree:

```bash
cd "$(bash agents/new_session_worktree.sh)"     # isolated checkout of origin/main in /tmp
# ... do all your work + run agents/finish_session.py here ...
```

A worktree shares the one `.git` (so pushes/pulls work normally) but has a private
working directory, so concurrent sessions can never clobber your uncommitted
files. Reclaim space later (clean worktrees only) with:

```bash
bash agents/cleanup_session_worktrees.sh
```

Suggested cron on the box: `*/30 * * * * cd ~/math && bash agents/cleanup_session_worktrees.sh`.
