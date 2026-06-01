#!/usr/bin/env bash
# cleanup_session_worktrees.sh — prune finished/stale per-session worktrees made
# by new_session_worktree.sh. Safe: only removes worktrees that are CLEAN (no
# uncommitted changes) so unfinished work is never discarded. Run from cron or by
# hand. Also prunes git's stale worktree registrations.
set -euo pipefail
REPO="$(git -C "$(dirname "$(readlink -f "$0")")/.." rev-parse --show-toplevel)"
git -C "$REPO" worktree prune
removed=0
for wt in /tmp/math-wt-*; do
  [ -d "$wt" ] || continue
  if [ -z "$(git -C "$wt" status --porcelain 2>/dev/null)" ]; then
    if git -C "$REPO" worktree remove --force "$wt" 2>/dev/null; then
      echo "removed clean worktree $wt"; removed=$((removed+1))
    fi
  else
    echo "KEEP (dirty/unfinished): $wt"
  fi
done
echo "pruned; removed $removed clean worktree(s)."
