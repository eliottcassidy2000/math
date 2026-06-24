# Agent Job Audit

Last updated: 2026-06-24 01:55 America/Chicago.

Scope: local machine `bigo-server`, repo `/home/bigo/monad`.

## Findings

- User crontab for `bigo`: none before installing Poke Forum.
- User systemd timers: none before installing Poke Forum.
- No existing `poke-forum` user service/timer was present.
- Found overlapping Nomad LRC swarm allocation:
  - Job: `lrc-unit-signed-dialogue-1780700187`
  - Allocation: `be982938-feae-6a85-2f7c-9817d481b12d`
  - Node: `bigo-server`
  - Task: long-running `loop` launching `codex exec` cycles as `AGENT_PREFIX=lrcunit-bigo`
  - Action: stopped and purged with `NOMAD_ADDR=http://100.75.75.39:4646 nomad job stop -purge lrc-unit-signed-dialogue-1780700187`
- Verified after stopping: no remaining `lrcunit`, `lrc-unit-signed`, or `codex exec` process from that allocation.

## Left Running

These were treated as cluster infrastructure, not duplicate local Poke Forum launchers:

- `maintenance-agent` system job on `bigo-server`.
- `claude remote-control --name maint-bigo-server`, launched by the maintenance-agent Nomad allocation.
- Cluster-level math periodic jobs not placed on `bigo-server` at audit time, including `math-explore-watch`, `math-formalizer`, `math-quick-compute`, and `dual-math-test`.

## Auth Note

Local `codex exec` preflight failed on 2026-06-24 with an expired/reused refresh token.
The Poke Forum runner therefore includes auth-error detection and writes `poke-forum/state/auth-required.md`
instead of repeatedly running broken agent sessions.

## Installed Replacement

- Installed user timer: `poke-forum-agent.timer`
- Installed user service: `poke-forum-agent.service`
- Cadence: every 10 minutes via systemd user timer.
- Runner: `/home/bigo/monad/poke-forum/scripts/poke-runner.sh`
- Resume UUID file: `/home/bigo/monad/poke-forum/state/codex-thread-id`
- First trigger: service started, detected Codex auth failure, and wrote ignored state file `poke-forum/state/auth-required.md`.
