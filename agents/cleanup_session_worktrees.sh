#!/usr/bin/env bash
# cleanup_session_worktrees.sh — inspect or remove temporary session worktrees.
#
# With no arguments this is a dry-run inventory. Pass exact /tmp/math-wt-* paths
# to remove those clean worktrees, or --all-clean only after checking that no
# still-starting session owns a clean directory. Dirty worktrees are never
# removed. Stale Git registrations are pruned in every mode.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(git -C "$SCRIPT_DIR/.." rev-parse --show-toplevel)"
git -C "$REPO" worktree prune

REMOVE=0
if [ "${1:-}" = "--all-clean" ]; then
  REMOVE=1
  shift
  set -- /tmp/math-wt-*
elif [ "$#" -gt 0 ]; then
  REMOVE=1
fi

removed=0
found=0
for wt in "${@:-}"; do
  [ -d "$wt" ] || continue
  case "$wt" in
    /tmp/math-wt-*) ;;
    *) echo "REFUSE (not a /tmp/math-wt-* path): $wt" >&2; continue ;;
  esac
  found=$((found+1))
  if [ -z "$(git -C "$wt" status --porcelain 2>/dev/null)" ]; then
    if [ "$REMOVE" -eq 0 ]; then
      echo "CLEAN (not removed): $wt"
    elif git -C "$REPO" worktree remove --force "$wt" 2>/dev/null; then
      echo "removed clean worktree $wt"
      removed=$((removed+1))
    fi
  else
    echo "KEEP (dirty/unfinished): $wt"
  fi
done

if [ "$REMOVE" -eq 0 ]; then
  # The quoted default above is empty, so inventory the glob explicitly.
  for wt in /tmp/math-wt-*; do
    [ -d "$wt" ] || continue
    found=$((found+1))
    if [ -z "$(git -C "$wt" status --porcelain 2>/dev/null)" ]; then
      echo "CLEAN (not removed): $wt"
    else
      echo "KEEP (dirty/unfinished): $wt"
    fi
  done
  [ "$found" -gt 0 ] || echo "no /tmp/math-wt-* session worktrees found"
  echo "dry run only; pass exact paths or --all-clean to remove clean entries"
else
  echo "pruned; removed $removed requested clean worktree(s)."
fi
