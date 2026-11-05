#!/usr/bin/env bash
# cleanup_session_worktrees.sh — compatibility wrapper around agents/worktree_hygiene.py.
#
#   bash agents/cleanup_session_worktrees.sh               # report every worktree (read-only)
#   bash agents/cleanup_session_worktrees.sh <path>        # remove that worktree if clean and merged
#   bash agents/cleanup_session_worktrees.sh --all-clean   # remove every clean, merged worktree
#
# The Python script checks that a worktree has no uncommitted or untracked files
# (a whitespace-only line-ending phantom is tolerated) and that its HEAD is
# contained in origin/main before removing it; dirty or unmerged worktrees are
# reported and kept. Worktrees are no longer restricted to /tmp (they must not
# live there: see new_session_worktree.sh).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
if [ "$#" -eq 0 ]; then
  exec python3 "$SCRIPT_DIR/worktree_hygiene.py"
elif [ "$1" = "--all-clean" ]; then
  exec python3 "$SCRIPT_DIR/worktree_hygiene.py" --remove
else
  status=0
  for wt in "$@"; do
    python3 "$SCRIPT_DIR/worktree_hygiene.py" --remove --path "$wt" || status=1
  done
  exit $status
fi
