# Agent Registry

All machines participating in the research network. Add your machine by running `python3 agents/processor.py --register` on a new machine, then editing your entry here.

| Machine ID | Description | Status | First session | Last seen |
|------------|-------------|--------|---------------|-----------|
| *(none yet — first machine registers by running --register)* | | | | |

---

## Network Architecture

Each machine has a directory at `agents/[machine-id]/` containing:
- `identity.md` — machine description and owner notes
- `inbox/` — messages addressed to this machine (written by OTHER machines)

`agents/broadcast/` contains messages addressed to all machines.

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
  git pull
  python3 agents/processor.py --check     (read incoming messages)
  python3 inbox/processor.py              (process human drops if any)
  [do research work]
  python3 agents/processor.py --send      (write end-of-session letter)
  git add -A && git commit -m "..." && git push
[Session end]
```
