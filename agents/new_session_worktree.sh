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
#   <tag>  short label (default: pid-timestamp); unsafe characters are replaced
#          and a numeric suffix is added if the local branch name already exists
#
# The worktree gets a unique local branch named codex/session-<tag> whose
# upstream is origin/main. Work in it, then run agents/finish_session.py as
# normal: the closer pushes HEAD to main and can verify against that upstream.
# Reclaim space afterwards with agents/cleanup_session_worktrees.sh.
set -euo pipefail
REPO="$(git -C "$(dirname "$(readlink -f "$0")")/.." rev-parse --show-toplevel)"
RAW_TAG="${1:-$$-$(date +%s)}"
TAG="$(printf '%s' "$RAW_TAG" | tr -c '[:alnum:]_-' '-')"
[ -n "$TAG" ] || TAG="$$-$(date +%s)"
WT="/tmp/math-wt-${TAG}"
if [ -e "$WT" ]; then
  if git -C "$WT" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "$WT"                                            # idempotent
    exit 0
  fi
  echo "ERROR: $WT exists but is not a git worktree" >&2
  exit 1
fi

BASE_BRANCH="codex/session-${TAG}"
BRANCH="$BASE_BRANCH"
SUFFIX=2
while git -C "$REPO" show-ref --verify --quiet "refs/heads/$BRANCH"; do
  BRANCH="${BASE_BRANCH}-${SUFFIX}"
  SUFFIX=$((SUFFIX + 1))
done

git -C "$REPO" fetch origin -q
git -C "$REPO" worktree add --quiet -b "$BRANCH" "$WT" origin/main
git -C "$WT" branch --set-upstream-to=origin/main "$BRANCH" >/dev/null
# carry the machine identity so the session keeps its name even if .machine-id is untracked
[ -f "$REPO/.machine-id" ] && cp "$REPO/.machine-id" "$WT/.machine-id" 2>/dev/null || true
echo "$WT"
