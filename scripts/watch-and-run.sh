#!/usr/bin/env bash
#
# Local session watcher: polls for new commits and launches Claude Code sessions.
# Run in background: nohup ./scripts/watch-and-run.sh &
#
# Requires: claude CLI installed, ANTHROPIC_API_KEY set, .machine-id configured
#
# Configuration via environment:
#   POLL_INTERVAL  - seconds between git fetch checks (default: 60)
#   MAX_TURNS      - max tool calls per session (default: 50)
#   MODEL          - claude model to use (default: claude-sonnet-4-6)
#   ONCE           - set to "1" to run one session and exit (no polling)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

# Ensure remote uses SSH (never prompts for credentials)
CURRENT_URL=$(git remote get-url origin)
if [[ "$CURRENT_URL" == https://* ]]; then
  SSH_URL=$(echo "$CURRENT_URL" | sed 's|https://github.com/|git@github.com:|')
  git remote set-url origin "$SSH_URL"
  echo "Switched remote to SSH: $SSH_URL"
fi

POLL_INTERVAL="${POLL_INTERVAL:-60}"
MAX_TURNS="${MAX_TURNS:-50}"
MODEL="${MODEL:-claude-sonnet-4-6}"
LOCKFILE="$REPO_ROOT/.session.lock"
LOGDIR="$REPO_ROOT/scripts/logs"

mkdir -p "$LOGDIR"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/watcher.log"
}

cleanup() {
  rm -f "$LOCKFILE"
  log "Watcher stopped."
}
trap cleanup EXIT

# Ensure .machine-id exists
if [ ! -f .machine-id ]; then
  echo "ERROR: .machine-id not found. Run setup first." >&2
  exit 1
fi
MACHINE_ID=$(cat .machine-id)

run_session() {
  local trigger_msg="$1"

  # Lockfile prevents concurrent sessions
  if [ -f "$LOCKFILE" ]; then
    log "Session already running (lockfile exists). Skipping."
    return
  fi
  echo $$ > "$LOCKFILE"

  local session_log="$LOGDIR/session-$(date +%Y%m%d-%H%M%S).log"
  log "Starting session: $trigger_msg"
  log "Session log: $session_log"

  # Pull latest
  git pull --rebase 2>&1 | tee -a "$session_log"

  # Run Claude Code
  claude -p "Follow the CLAUDE.md startup sequence exactly. ${trigger_msg}. When finished, write a session letter, update SESSION-LOG.md, commit, and push." \
    --model "$MODEL" \
    --max-turns "$MAX_TURNS" \
    --allowedTools "Bash(git *),Bash(python3 *),Read,Write,Edit,Glob,Grep" \
    2>&1 | tee -a "$session_log"

  local exit_code=${PIPESTATUS[0]}

  # Push any uncommitted changes
  if [ -n "$(git status --porcelain)" ]; then
    git add -A
    git diff --cached --quiet || git commit -m "${MACHINE_ID}-$(date +%Y-%m-%d)-auto: automated session

Co-Authored-By: Claude <noreply@anthropic.com>"
    git push 2>&1 | tee -a "$session_log" || log "Push failed — will retry next cycle"
  fi

  rm -f "$LOCKFILE"
  log "Session completed (exit $exit_code)"
}

# Single-run mode
if [ "${ONCE:-0}" = "1" ]; then
  run_session "Work on the highest-priority open question or respond to inbox messages"
  exit 0
fi

# Polling loop
log "Watcher started for machine '$MACHINE_ID' (poll every ${POLL_INTERVAL}s)"
LAST_REMOTE_SHA=$(git rev-parse origin/main 2>/dev/null || echo "none")

while true; do
  git fetch origin main --quiet 2>/dev/null || true
  CURRENT_REMOTE_SHA=$(git rev-parse origin/main 2>/dev/null || echo "none")

  if [ "$CURRENT_REMOTE_SHA" != "$LAST_REMOTE_SHA" ]; then
    # Get info about the new commit(s)
    NEW_COMMITS=$(git log --oneline "$LAST_REMOTE_SHA".."$CURRENT_REMOTE_SHA" 2>/dev/null || echo "new commits")
    LAST_AUTHOR=$(git log -1 --format='%an' origin/main 2>/dev/null || echo "unknown")

    # Don't trigger on our own commits
    if [ "$LAST_AUTHOR" = "$(git config user.name)" ]; then
      log "Skipping: commit by self ($LAST_AUTHOR)"
      LAST_REMOTE_SHA="$CURRENT_REMOTE_SHA"
    else
      log "New commits detected from '$LAST_AUTHOR':"
      log "$NEW_COMMITS"
      LAST_REMOTE_SHA="$CURRENT_REMOTE_SHA"
      run_session "New commits detected. Check messages and continue research"
    fi
  fi

  sleep "$POLL_INTERVAL"
done
