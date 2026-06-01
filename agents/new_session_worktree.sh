#!/usr/bin/env bash
# new_session_worktree.sh — give a math session its OWN isolated git worktree.
#
# Why: multiple remote-control sessions sharing the one ~/math checkout clobber
# each other's UNCOMMITTED files (one session's `git add -A` sweeps another's
# work; see CONCURRENT-SESSIONS.md). A worktree is a separate working directory
# on the same .git, so uncommitted files never collide. (Nomad-launched sessions
# already isolate via /tmp clones; this is the equivalent for remote-control.)
#
# Usage:
#   cd "$(bash agents/new_session_worktree.sh [tag])"     # then work + finish_session
#   <tag>  short unique label (default: pid-timestamp)
#
# Work in the worktree, then run agents/finish_session.py as normal — it pushes
# HEAD:main with the union-merge rebase-retry loop. Reclaim space afterwards with
# agents/cleanup_session_worktrees.sh.
set -euo pipefail
REPO="$(git -C "$(dirname "$(readlink -f "$0")")/.." rev-parse --show-toplevel)"
TAG="${1:-$$-$(date +%s)}"
WT="/tmp/math-wt-${TAG}"
if [ -e "$WT" ]; then echo "$WT"; exit 0; fi          # idempotent
git -C "$REPO" fetch origin -q || true
git -C "$REPO" worktree add --detach "$WT" origin/main >/dev/null 2>&1
# carry the machine identity so the session keeps its name even if .machine-id is untracked
[ -f "$REPO/.machine-id" ] && cp "$REPO/.machine-id" "$WT/.machine-id" 2>/dev/null || true
echo "$WT"
