#!/usr/bin/env bash
set -uo pipefail

REPO_DIR="${POKE_REPO_DIR:-/home/bigo/math}"
FORUM_DIR="$REPO_DIR/poke-forum"
STATE_DIR="$FORUM_DIR/state"
RUNTIME_DIR="$FORUM_DIR/runtime"
LOG_DIR="$FORUM_DIR/logs"
THREAD_FILE="${POKE_THREAD_ID_FILE:-$STATE_DIR/codex-thread-id}"
CODEX_BIN="${POKE_CODEX_BIN:-/home/bigo/.local/bin/codex}"
LOCK_FILE="$STATE_DIR/poke-runner.lock"
TIMEOUT_SECONDS="${POKE_ROLE_TIMEOUT_SECONDS:-2400}"

mkdir -p "$STATE_DIR" "$RUNTIME_DIR" "$LOG_DIR"

log() {
  printf '[poke-runner %s] %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$*"
}

die() {
  log "ERROR: $*"
  exit 1
}

if [ -n "${CODEX_THREAD_ID:-}" ] && [ ! -s "$THREAD_FILE" ]; then
  printf '%s\n' "$CODEX_THREAD_ID" > "$THREAD_FILE"
fi

[ -x "$CODEX_BIN" ] || die "Codex binary not executable: $CODEX_BIN"
[ -d "$REPO_DIR/.git" ] || die "Repo dir is not a git checkout: $REPO_DIR"

COMMON_PROMPT="$FORUM_DIR/prompts/common.md"
[ -f "$COMMON_PROMPT" ] || die "Missing common prompt: $COMMON_PROMPT"

auth_blocked() {
  local out="$1"
  rg -qi 'access token could not be refreshed|refresh token.*used|token_expired|401 Unauthorized|codex login|Please log out and sign in again' "$out"
}

record_auth_required() {
  local out="$1"
  {
    printf '# Codex Auth Required\n\n'
    printf -- '- Detected: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf -- '- Codex binary: `%s`\n' "$CODEX_BIN"
    printf -- '- Fix: run `codex login` as user `bigo`, then the next timer cycle can resume.\n\n'
    printf '## Last Output Tail\n\n```text\n'
    tail -80 "$out"
    printf '\n```\n'
  } > "$STATE_DIR/auth-required.md"
}

compose_prompt() {
  local role="$1"
  local role_prompt="$2"
  local prompt_file="$RUNTIME_DIR/prompt-$role-$(date -u '+%Y%m%dT%H%M%SZ').md"
  {
    printf '# Recurring Poke Forum Session\n\n'
    printf -- '- UTC: %s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf -- '- Repo: `%s`\n' "$REPO_DIR"
    printf -- '- Role: `%s`\n' "$role"
    if [ -s "$THREAD_FILE" ]; then
      printf -- '- Resume UUID: `%s`\n' "$(tr -d '[:space:]' < "$THREAD_FILE")"
    fi
    printf '\n'
    cat "$COMMON_PROMPT"
    printf '\n\n'
    cat "$role_prompt"
    printf '\n'
  } > "$prompt_file"
  printf '%s\n' "$prompt_file"
}

run_role() {
  local role="$1"
  local role_prompt="$2"
  local prompt_file out_file thread rc

  [ -f "$role_prompt" ] || die "Missing role prompt: $role_prompt"
  prompt_file="$(compose_prompt "$role" "$role_prompt")"
  out_file="$RUNTIME_DIR/out-$role-$(date -u '+%Y%m%dT%H%M%SZ').txt"

  log "starting role=$role"
  cd "$REPO_DIR" || return 2

  if [ -s "$THREAD_FILE" ]; then
    thread="$(tr -d '[:space:]' < "$THREAD_FILE")"
    timeout --signal=TERM "$TIMEOUT_SECONDS" \
      "$CODEX_BIN" --search exec -C "$REPO_DIR" \
      resume --dangerously-bypass-approvals-and-sandbox --skip-git-repo-check \
      "$thread" - < "$prompt_file" > "$out_file" 2>&1
    rc=$?
  else
    timeout --signal=TERM "$TIMEOUT_SECONDS" \
      "$CODEX_BIN" --search exec -C "$REPO_DIR" \
      --dangerously-bypass-approvals-and-sandbox --skip-git-repo-check \
      - < "$prompt_file" > "$out_file" 2>&1
    rc=$?
  fi

  if auth_blocked "$out_file"; then
    log "codex auth failure detected during role=$role"
    record_auth_required "$out_file"
    return 10
  fi

  if [ "$rc" -eq 0 ]; then
    rm -f "$STATE_DIR/auth-required.md"
    log "completed role=$role"
  else
    log "role=$role exited rc=$rc; see $out_file"
  fi

  return "$rc"
}

main() {
  exec 9>"$LOCK_FILE"
  if ! flock -n 9; then
    log "another poke-runner cycle is still active; skipping"
    exit 0
  fi

  {
    printf 'last_start_utc=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf 'repo=%s\n' "$REPO_DIR"
  } > "$STATE_DIR/last-run.txt"

  local rc overall=0

  if [ "${POKE_DRY_RUN:-0}" = 1 ]; then
    compose_prompt "coordinator" "$FORUM_DIR/prompts/coordinator.md" >/dev/null
    compose_prompt "math-exploration" "$FORUM_DIR/prompts/math-exploration.md" >/dev/null
    compose_prompt "math-investigation" "$FORUM_DIR/prompts/math-investigation.md" >/dev/null
    log "dry run wrote prompts under $RUNTIME_DIR"
    exit 0
  fi

  run_role "coordinator" "$FORUM_DIR/prompts/coordinator.md"
  rc=$?
  [ "$rc" -eq 10 ] && exit 0
  [ "$rc" -ne 0 ] && overall="$rc"

  run_role "math-exploration" "$FORUM_DIR/prompts/math-exploration.md"
  rc=$?
  [ "$rc" -eq 10 ] && exit 0
  [ "$rc" -ne 0 ] && overall="$rc"

  run_role "math-investigation" "$FORUM_DIR/prompts/math-investigation.md"
  rc=$?
  [ "$rc" -eq 10 ] && exit 0
  [ "$rc" -ne 0 ] && overall="$rc"

  {
    printf 'last_finish_utc=%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
    printf 'exit_code=%s\n' "$overall"
  } >> "$STATE_DIR/last-run.txt"

  exit "$overall"
}

main "$@" > >(tee -a "$LOG_DIR/runner-$(date -u '+%Y%m%d').log") 2>&1
