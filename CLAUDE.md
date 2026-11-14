# Claude Code adapter

Claude Code reads this file automatically. The shared research and Git contract
is [`AGENTS.md`](AGENTS.md); read and follow it first. Do not reconstruct project
status from this adapter.

## Identity and bounded startup

If `.machine-id` exists, use its contents in session IDs and provenance. If it
does not, use a short temporary identity and do not invent a registered agent.

```bash
python3 agents/start_session.py --topic "<your task>"
python3 agents/processor.py --status --peek --limit 8
```

Read targeted messages relevant to the task. The inbox and broadcast archives
are very large; do not run an unbounded historical message dump as warm-up.
Process `inbox/` only when the owner supplied new intake files or the current
task asks for it.

## Owner directives (2026-11-13)

No worktrees: work on `main` in the given directory (no `EnterWorktree`, no
`use_worktree: true`); stage only files you wrote. Ideas, not audits: no
auditor subagents or re-verification passes; push to close proofs.

## Claude-specific handoff

At a useful checkpoint, follow `AGENTS.md`; at final close-out, a Claude
session may send a concise letter and publish with:

```bash
python3 agents/finish_session.py \
  --to all \
  --subject "[instance-id]: [one-line result]" \
  --body "What changed; proof status; decisive evidence; next obligations." \
  --commit-msg "[instance-id]: [one-line git summary]" \
  --paths <every file this session wrote>
```

`--paths` is required in the shared checkout (the closer commits only those
paths plus its letter files). There is deliberately no blocking Stop hook.


For setup, messaging, and recovery details, read
[`agents/README.md`](agents/README.md).
