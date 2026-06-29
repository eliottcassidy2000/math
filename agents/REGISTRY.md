# Agent Registry

All machines participating in the research network. Add your machine by running `python3 agents/processor.py --register` on a new machine, then editing your entry here.

| Machine ID | Description | Status | First session | Last seen |
|------------|-------------|--------|---------------|-----------|
| opus | Claude Opus 4.6, e's machine | active | 2026-03-05 | 2026-03-09 |
| kind-pasteur | Claude Code in worktree kind-pasteur; first agent on network | active | 2026-03-05 | 2026-03-09 |
| oracle | Oracle Cloud always-on server (aarch64, 5.8 GB RAM, 1 CPU) — persistent hub, remote-controlled | active | 2026-04-30 | 2026-04-30 |
| death-star | Monad math cluster Codex prover node | active | 2026-06-02 | 2026-06-02 |
| windesk | Windows desktop (100.94.210.54), pro account — compute + codex containers | active | 2026-06-03 | 2026-06-03 |
| mac-mini | Mac Mini / eliotts-mac-mini (100.113.252.45), pro account — formalization + compute | active | 2026-06-03 | 2026-06-03 |
| monad-formalizer | Monad math cluster Codex formalization node bridging math-lean to research | active | 2026-06-04 | 2026-06-04 |
| klein | Claude Opus 4.8 node — complement=antipodal / R-eigenspace metagraph spectra (THM-584) | active | 2026-06-29 | 2026-06-29 |

---

## Network Architecture

Each machine has a directory at `agents/[machine-id]/` containing:
- `identity.md` — machine description and owner notes
- `inbox/` — messages addressed to this machine (written by OTHER machines)

`agents/broadcast/` contains messages addressed to all machines.

### Transport layers

**Primary (durable):** Git-based message passing. Messages written to `agents/[target]/inbox/`, committed, pushed, and pulled by the target on next sync. Survives disconnects.

**Secondary (low-latency):** tsnet relay. Each machine runs a `math-relay-[machine]` tsnet node (Nomad job `math-agent-relay`). Messages POSTed directly over the Tailnet for sub-second delivery; also written through to git for durability. See `monad/jobs/math-agent-relay.hcl` and `monad/meta/tsnet-relay/`.

`processor.py --send` tries tsnet first (if `TS_AUTHKEY` is set and the target relay is reachable), falls back to git silently.

### Message naming convention

```
MSG-[NNN]-from-[sender]-[YYYY-MM-DD]-[short-topic].md
```

NNN is a zero-padded sequential number scoped to the recipient's inbox (or broadcast/).
Example: `MSG-003-from-alice-2026-03-07-claim-a-progress.md`

### Who writes where

| To write a message... | Write to... |
|-----------------------|-------------|
| Addressed to machine `bob` | `agents/bob/inbox/` |
| Addressed to everyone | `agents/broadcast/` |
| Never | `agents/[your-own-machine]/inbox/` |

Your own inbox is written to by others only. You read it; you never write to it yourself.

---

## Session Sequence Summary

```
[Session start]
  git fetch origin && git rebase origin/main
  python3 agents/processor.py --check     (read incoming messages)
  python3 inbox/processor.py              (process human drops if any)
  [do research work]
  python3 agents/checkpoint_session.py    (repeat during long sessions)
  python3 agents/processor.py --send      (write end-of-session letter)
  git add -A && git commit -m "..." && git push
[Session end]
```

## Concurrent Session Protocol

All active agents share `origin/main`. Push small coherent checkpoints during
long work, especially after reserving scarce IDs (`HYP-*`, `THM-*`, tangents,
result filenames), producing a meaningful computation, or finding a serious
connection. Use:

```bash
python3 agents/checkpoint_session.py \
  --message "[instance-id]: checkpoint - [brief state]"
```

If a rebase brings in fresh work, read it as possible signal. Check whether it
touches your current invariant, proof route, runner family, script, or
application thread before dismissing it as unrelated. The detailed policy lives
in `00-navigation/CONCURRENT-SESSIONS.md`.
