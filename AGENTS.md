# Repository Instructions for Codex

This repository is a persistent research workspace. Do not end a session from
this checkout with local work stranded on the machine.

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
