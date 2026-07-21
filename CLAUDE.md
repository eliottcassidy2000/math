# Claude Code adapter

Claude Code reads this file automatically. The shared research and Git contract
is [`AGENTS.md`](AGENTS.md); read and follow it first. Do not reconstruct project
status from this adapter.

## Identity and bounded startup

If `.machine-id` exists, use its contents in session IDs and provenance. If it
does not, use a short temporary identity and do not invent a registered agent.

```bash
python3 agents/start_session.py --topic "<your task>"
python3 agents/processor.py --status
```

Read targeted messages relevant to the task. The inbox and broadcast archives
are very large; do not run an unbounded historical message dump as warm-up.
Process `inbox/` only when the owner supplied new intake files or the current
task asks for it.

## Claude-specific handoff

Use an isolated worktree whenever possible. At a useful checkpoint, follow
`AGENTS.md`; at final close-out, a Claude session may send a concise letter and
publish with:

```bash
python3 agents/finish_session.py \
  --to all \
  --subject "[instance-id]: [one-line result]" \
  --body "What changed; proof status; decisive evidence; next obligations." \
  --commit-msg "[instance-id]: [one-line git summary]"
```

The closer is safe only when the worktree contains this session's changes.
There is deliberately no blocking Stop hook: successful close-out is verified
by Git state, not by global timestamp files or another agent's messages.

For setup, messaging, recovery, and worktree details, read
[`agents/README.md`](agents/README.md).
