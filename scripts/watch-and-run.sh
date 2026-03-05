#!/usr/bin/env bash
#
# Continuous Claude session runner.
#
# Watches for new commits pushed by other machines. When detected,
# pulls and launches a local Claude Code session that does math work,
# commits, and pushes — triggering watchers on other machines.
#
# Setup (once per machine):
#   1. Ensure .machine-id exists and agent is registered
#   2. Ensure SSH auth works: ssh -T git@github.com
#   3. Run: ./scripts/watch-and-run.sh
#
# Environment variables:
#   POLL_INTERVAL  seconds between checks (default: 30)
#   MAX_TURNS      max tool calls per session (default: 50)
#   MODEL          model to use (default: claude-sonnet-4-6)
#   ONCE           set to 1 to run one session then exit
#   COOLDOWN       seconds to wait after own push before polling again (default: 10)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

POLL_INTERVAL="${POLL_INTERVAL:-30}"
MAX_TURNS="${MAX_TURNS:-50}"
MODEL="${MODEL:-claude-sonnet-4-6}"
COOLDOWN="${COOLDOWN:-10}"
LOCKFILE="$REPO_ROOT/.session.lock"
LOGDIR="$REPO_ROOT/scripts/logs"

mkdir -p "$LOGDIR"

# ── Ensure SSH remote (never prompts for credentials) ──────────────────────
CURRENT_URL=$(git remote get-url origin)
if [[ "$CURRENT_URL" == https://* ]]; then
  SSH_URL=$(echo "$CURRENT_URL" | sed 's|https://github.com/|git@github.com:|')
  git remote set-url origin "$SSH_URL"
  echo "Switched remote to SSH: $SSH_URL"
fi

# ── Ensure identity ────────────────────────────────────────────────────────
if [ ! -f .machine-id ]; then
  echo "ERROR: .machine-id not found."
  echo "Create one: echo 'my-machine-name' > .machine-id"
  echo "Then register: python3 agents/processor.py --register"
  exit 1
fi
MACHINE_ID=$(cat .machine-id)
GIT_AUTHOR=$(git config user.name 2>/dev/null || echo "$MACHINE_ID")

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOGDIR/watcher.log"
}

cleanup() {
  rm -f "$LOCKFILE"
  log "Watcher stopped."
}
trap cleanup EXIT

# ── Session runner ─────────────────────────────────────────────────────────
run_session() {
  local trigger_msg="$1"

  if [ -f "$LOCKFILE" ]; then
    local lock_pid
    lock_pid=$(cat "$LOCKFILE" 2>/dev/null || echo "")
    if [ -n "$lock_pid" ] && kill -0 "$lock_pid" 2>/dev/null; then
      log "Session already running (pid $lock_pid). Skipping."
      return
    fi
    log "Stale lockfile found. Removing."
    rm -f "$LOCKFILE"
  fi

  echo $$ > "$LOCKFILE"
  local session_log="$LOGDIR/session-$(date +%Y%m%d-%H%M%S).log"
  log "Starting session: $trigger_msg"

  # Pull latest
  git pull --rebase 2>&1 | tee -a "$session_log"

  # Run Claude Code in non-interactive mode
  local prompt="You are machine '${MACHINE_ID}'. Follow the CLAUDE.md startup sequence exactly. ${trigger_msg} When finished, write a session letter (using agents/processor.py --send), update SESSION-LOG.md, then commit and push your work."

  claude -p "$prompt" \
    --model "$MODEL" \
    --max-turns "$MAX_TURNS" \
    --allowedTools "Bash(git *),Bash(python3 *),Bash(ls *),Bash(cat .machine-id),Read,Write,Edit,Glob,Grep" \
    2>&1 | tee -a "$session_log"

  local exit_code=${PIPESTATUS[0]}

  # Safety net: push any uncommitted changes claude left behind
  if [ -n "$(git status --porcelain)" ]; then
    git add -A
    git diff --cached --quiet || git commit -m "${MACHINE_ID}-$(date +%Y-%m-%d)-auto: session cleanup

Co-Authored-By: Claude <noreply@anthropic.com>"
  fi

  # Push (with retry on conflict)
  git push 2>&1 | tee -a "$session_log" || {
    log "Push failed, rebasing and retrying..."
    git pull --rebase 2>&1 | tee -a "$session_log"
    git push 2>&1 | tee -a "$session_log" || log "Push failed after retry"
  }

  rm -f "$LOCKFILE"
  log "Session completed (exit $exit_code). Log: $session_log"
  sleep "$COOLDOWN"
}

# ── Single-run mode ────────────────────────────────────────────────────────
if [ "${ONCE:-0}" = "1" ]; then
  run_session "Work on the highest-priority open question or respond to inbox messages."
  exit 0
fi

# ── Polling loop ───────────────────────────────────────────────────────────
log "Watcher started for '${MACHINE_ID}' (poll every ${POLL_INTERVAL}s, model ${MODEL})"
LAST_SHA=$(git rev-parse origin/main 2>/dev/null || echo "none")

while true; do
  git fetch origin main --quiet 2>/dev/null || { sleep "$POLL_INTERVAL"; continue; }
  CURRENT_SHA=$(git rev-parse origin/main 2>/dev/null || echo "none")

  if [ "$CURRENT_SHA" != "$LAST_SHA" ]; then
    LAST_AUTHOR=$(git log -1 --format='%an' origin/main 2>/dev/null || echo "unknown")
    LAST_MSG=$(git log -1 --format='%s' origin/main 2>/dev/null || echo "")

    if [ "$LAST_AUTHOR" = "$GIT_AUTHOR" ]; then
      log "New commit by self — skipping. ($LAST_MSG)"
    else
      log "New commit by '$LAST_AUTHOR': $LAST_MSG"
      run_session "New commits from '${LAST_AUTHOR}'. Check your inbox and continue research."
    fi
    LAST_SHA="$CURRENT_SHA"
  fi

  sleep "$POLL_INTERVAL"
done
