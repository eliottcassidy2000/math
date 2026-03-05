#!/usr/bin/env python3
"""
check_session_closed.py — Stop hook for Claude Code.

Runs when Claude ends its response. Checks whether the session was properly
closed (letter sent + pushed). If not, prints a warning to stderr so it
appears in the Claude Code UI.

Detection heuristic:
- Look at recent git log for a push since session start (tracked via .session-state.json)
- Check if any outbox message was written in this session
- If neither happened, warn Claude.

This script is NON-BLOCKING: it exits 0 regardless, so it never prevents
Claude from ending. It only prints a reminder.
"""

import json
import subprocess
from datetime import datetime
from pathlib import Path

REPO_ROOT  = Path(__file__).resolve().parent.parent
STATE_FILE = REPO_ROOT / "agents" / ".session-state.json"
AGENTS_DIR = REPO_ROOT / "agents"


def get_state():
    if STATE_FILE.exists():
        try:
            return json.loads(STATE_FILE.read_text())
        except Exception:
            pass
    return {}


def save_state(state):
    STATE_FILE.write_text(json.dumps(state, indent=2))


def get_git_head():
    r = subprocess.run(
        "git rev-parse HEAD", shell=True, cwd=REPO_ROOT,
        capture_output=True, text=True
    )
    return r.stdout.strip() if r.returncode == 0 else None


def get_remote_head():
    r = subprocess.run(
        "git rev-parse @{u}", shell=True, cwd=REPO_ROOT,
        capture_output=True, text=True
    )
    return r.stdout.strip() if r.returncode == 0 else None


def count_outbox_messages(since_time: str) -> int:
    """Count messages written to any inbox since session started."""
    count = 0
    for inbox in AGENTS_DIR.glob("*/inbox/MSG-*.md"):
        try:
            mtime = inbox.stat().st_mtime
            mdt = datetime.fromtimestamp(mtime).isoformat()
            if mdt > since_time:
                count += 1
        except Exception:
            pass
    # Also check broadcast
    broadcast = AGENTS_DIR / "broadcast"
    for msg in broadcast.glob("MSG-*.md"):
        try:
            mtime = msg.stat().st_mtime
            mdt = datetime.fromtimestamp(mtime).isoformat()
            if mdt > since_time:
                count += 1
        except Exception:
            pass
    return count


def main():
    state = get_state()
    now = datetime.now().isoformat()

    # First call in this session: record start state
    if "session_start" not in state:
        state["session_start"] = now
        state["start_commit"] = get_git_head()
        state["pushed"] = False
        state["message_sent"] = False
        save_state(state)
        # Don't warn on first call — session is just starting
        return

    # Check if pushed
    local  = get_git_head()
    remote = get_remote_head()
    if local and remote and local == remote:
        state["pushed"] = True

    # Check if a message was sent during this session
    since = state.get("session_start", "2000-01-01")
    msgs  = count_outbox_messages(since)
    if msgs > 0:
        state["message_sent"] = True

    save_state(state)

    pushed       = state.get("pushed", False)
    message_sent = state.get("message_sent", False)

    if pushed and message_sent:
        # Clean state for next session
        save_state({})
        return

    # Print warning to stderr — appears in Claude Code as a hook message
    missing = []
    if not message_sent:
        missing.append("send a session letter (python3 agents/finish_session.py)")
    if not pushed:
        missing.append("push to remote (git push, or use finish_session.py)")

    if missing:
        print(
            "\n\u26a0\ufe0f  SESSION NOT CLOSED PROPERLY\n"
            "Before ending this session, you must:\n" +
            "\n".join(f"  - {m}" for m in missing) +
            "\n\nRun: python3 agents/finish_session.py --to all "
            '--subject "..." --body "..." --commit-msg "..."'
            "\n",
            flush=True
        )


if __name__ == "__main__":
    main()
