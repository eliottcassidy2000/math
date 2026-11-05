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
#
# Location: a SIBLING of the repository (<parent>/math-wt-<tag>), never /tmp:
# /tmp worktrees are deleted by the periodic cleanup and by reboots, which lost
# untracked lane files, and the ones that survived (thirty full 3 GB checkouts)
# filled the whole data volume on 2026-11-04 and stopped every session. Override
# the parent with MATH_WT_PARENT=<dir> if you must.
#
# Every worktree is a full checkout. REMOVE YOURS AT CLOSE-OUT, after
# finish_session.py has pushed and from outside the worktree:
#     python3 agents/worktree_hygiene.py --remove --path <worktree>
# (it refuses dirty or unmerged worktrees; the plain call reports them all).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd -P)"
REPO="$(git -C "$SCRIPT_DIR/.." rev-parse --show-toplevel)"
RAW_TAG="${1:-$$-$(date +%s)}"
TAG="$(printf '%s' "$RAW_TAG" | tr -c '[:alnum:]_-' '-')"
[ -n "$TAG" ] || TAG="$$-$(date +%s)"
PARENT="${MATH_WT_PARENT:-$(dirname "$REPO")}"
WT="${PARENT}/math-wt-${TAG}"
case "$WT" in /tmp/*|/private/tmp/*) echo "ERROR: refusing to create a session worktree under /tmp ($WT); set MATH_WT_PARENT to a persistent directory" >&2; exit 1;; esac
# hygiene warnings: many worktrees or little free space
NWT=$(git -C "$REPO" worktree list --porcelain | grep -c '^worktree ')
FREE_GB=$(df -k "$PARENT" | awk 'NR==2 {printf "%d", $4/1048576}')
if [ "$NWT" -ge 6 ]; then echo "WARNING: $NWT worktrees are registered; remove merged ones with: python3 agents/worktree_hygiene.py --remove" >&2; fi
if [ "${FREE_GB:-0}" -lt 20 ]; then echo "WARNING: only ${FREE_GB} GB free on the volume of $PARENT; run python3 agents/worktree_hygiene.py --remove before working" >&2; fi
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
echo "created $WT; remove it at close-out (after the push, from outside it): python3 agents/worktree_hygiene.py --remove --path $WT" >&2
echo "$WT"
