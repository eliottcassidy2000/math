# Repository Instructions for Codex

This repository is a persistent research workspace. Do not end a session from
this checkout with local work stranded on the machine.

## Concurrent Mainline Operating Mode

Treat `origin/main` as the shared live research surface. Long investigations
should push small coherent commits repeatedly, not only at final close-out. Push
after claiming scarce IDs or filenames, after meaningful computations, after
confirming/refuting a hypothesis, and before risky rebases. Prefer:

```bash
python3 agents/checkpoint_session.py \
  --message "[instance-id]: checkpoint - [brief state]"
```

Claim hypothesis numbers, theorem numbers, tangent IDs, result filenames, and
session-log territory early with honest stubs when a session will need them.
Make the stub explicit about what is known, what is still missing, and why the
namespace is being reserved.

When `git fetch`, `git pull`, or `git rebase` reveals new work from another
agent, treat that work as signal. Before continuing, check whether it connects
to the current invariant, proof route, script, runner family, or application
thread. If it does, integrate the connection in the relevant hypothesis,
backlog, reflection, or session log. If it does not, leave it intact and
continue.

Read `00-navigation/CONCURRENT-SESSIONS.md` for the full playbook.

## Research Default: Tournament Analysis

When adding or updating computational research scripts, aim to include
Tournament Analysis whenever it is meaningfully available. Declare the pairwise
observable, the switch/gauge that turns pair data into a binary relation, and
the tie Hamiltonian path. Report tournament fingerprints such as score
histograms, directed cycles, SCCs, edge flips, and Hamiltonian-path counts when
feasible. If a script cannot use Tournament Analysis cleanly, state why in its
methodology or notes.

## Assumption-Challenge Default

For LRC or Tournament Analysis sessions, do not assume tournament vertices must
be runners or arcs. Before settling on a mapping, explicitly consider alternate
vertex sets such as runners, gaps, fixed circle sections, section boundaries,
wall-crossing events, residues, cover arcs, Fourier modes, matroid circuits,
and proof obligations. Record which LRC predicate the quotient preserves, what
information it destroys, and at least one challenged assumption in the relevant
reflection, hypothesis, script note, or session log when the session is
exploratory.

## Mandatory GitHub Close-Out

Before the final response in every future Codex session from this repository:

1. Run `git status --short --branch`.
2. Stage all intentional repo changes with `git add -A`.
3. Commit if there is anything to commit.
4. Push the current branch to its upstream with `git push`. If no upstream is
   set, use `git push -u origin $(git branch --show-current)`.
5. Verify the branch is no longer ahead of its upstream.

If the push fails, do not treat the session as complete. Report the exact
blocker and leave the work committed locally so the next session can push or
repair it.

For Claude-style sessions, prefer the repo closer:

```bash
python3 agents/finish_session.py \
  --to all \
  --subject "[instance-id]: [one-line summary]" \
  --body "Detailed findings and handoff." \
  --commit-msg "[instance-id]: [one-line git summary]"
```

That script sends the session letter, commits, and pushes the current branch.
