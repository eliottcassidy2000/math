#!/usr/bin/env python3
"""
finish_session.py — End-of-session closer for Claude agents.

Performs all three mandatory end-of-session steps in one command:
  1. Sends a summary message to another agent (or broadcast)
  2. git add -A && git commit
  3. git push to the current branch's upstream (with rebase-retry on conflict)

Usage (interactive):
    python3 agents/finish_session.py

Usage (non-interactive, for Claude):
    python3 agents/finish_session.py \\
        --to all \\
        --subject "Session summary: <one-line>" \\
        --body-file /tmp/letter.md \\
        --commit-msg "kind-pasteur-2026-03-05-S2: <one-line summary>"

If --body-file is omitted you can use --body for short inline messages.

Exit codes:
    0  all steps succeeded
    1  argument error
    2  message delivery failed
    3  git commit failed
    4  git push failed (even after rebase retry)
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
PYTHON    = sys.executable  # use same Python that's running this script


def run(cmd_list, cwd=REPO_ROOT, capture=False, extra_env=None):
    """Run a command given as a list (no shell injection)."""
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    if extra_env:
        env.update(extra_env)
    result = subprocess.run(
        cmd_list, cwd=cwd,
        capture_output=capture, text=True, env=env
    )
    return result


def run_shell(cmd_str, cwd=REPO_ROOT, capture=False):
    """Run a shell string command (git, etc.)."""
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    result = subprocess.run(
        cmd_str, cwd=cwd, shell=True,
        capture_output=capture, text=True, env=env
    )
    return result


def send_message(to, subject, body_source, machine_id):
    """Deliver the session letter via processor.py."""
    processor = str(REPO_ROOT / "agents" / "processor.py")

    cmd = [PYTHON, processor, "--send", "--to", to, "--subject", subject]
    if body_source.startswith("file:"):
        cmd += ["--body-file", body_source[5:]]
    else:
        cmd += ["--body", body_source]

    result = run(cmd, capture=True)
    if result.returncode != 0:
        print(f"ERROR: Message delivery failed:\n{result.stderr}", file=sys.stderr)
        return False
    print(result.stdout, end="")
    return True


def current_branch():
    """Return the checked-out branch name, or None for detached HEAD."""
    r = run(["git", "branch", "--show-current"], capture=True)
    branch = (r.stdout or "").strip()
    return branch or None


def upstream_ref():
    """Return the current branch upstream ref, e.g. origin/main, or None."""
    r = run(
        ["git", "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"],
        capture=True,
    )
    if r.returncode != 0:
        return None
    upstream = (r.stdout or "").strip()
    return upstream or None


def push_current_branch(max_attempts=6):
    """Push HEAD to the upstream, with a fetch/rebase retry LOOP for concurrent-
    push races. Append-heavy coordination files auto-resolve via the `merge=union`
    driver in .gitattributes; on any conflict the union driver cannot resolve, the
    rebase is ABORTED so the working tree is left clean (never mid-rebase) and the
    local commit stays intact for a manual retry."""
    branch = current_branch()
    upstream = upstream_ref()
    if upstream:
        remote, remote_branch = upstream.split("/", 1)
    elif branch:
        remote, remote_branch = "origin", branch
    else:  # detached HEAD (e.g. after a prior aborted rebase): default to origin/main
        remote, remote_branch = "origin", "main"

    for attempt in range(1, max_attempts + 1):
        print(f"\n-- git push (to {remote}/{remote_branch}, attempt {attempt}/{max_attempts}) --")
        if branch and upstream_ref() is None:
            r = run(["git", "push", "-u", remote, branch], capture=True)
        else:
            r = run(["git", "push", remote, f"HEAD:{remote_branch}"], capture=True)
        if r.returncode == 0:
            print(r.stdout or "  Pushed successfully.")
            return True

        print(f"  Push rejected (remote moved) — fetch + rebase onto {remote}/{remote_branch}...")
        rf = run(["git", "fetch", remote], capture=True)
        if rf.returncode != 0:
            print(f"ERROR: git fetch failed.\n{rf.stdout}\n{rf.stderr}", file=sys.stderr)
            return False

        rr = run(["git", "rebase", f"{remote}/{remote_branch}"], capture=True)
        if rr.returncode != 0:
            print("  Rebase conflict the union driver could not auto-resolve —"
                  " aborting rebase to keep a clean tree.")
            run(["git", "rebase", "--abort"], capture=True)
            print(f"ERROR: rebase conflict on close.\n{rr.stdout}\n{rr.stderr}", file=sys.stderr)
            print("  Your commit is intact on a clean branch; resolve the file(s) and "
                  "re-run finish_session, or push manually.", file=sys.stderr)
            return False
        # rebase succeeded (possibly via union merge) — loop and try the push again.

    print(f"ERROR: push still failing after {max_attempts} attempts (heavy contention).",
          file=sys.stderr)
    return False


def verify_not_ahead():
    """Verify that HEAD is not ahead of its configured upstream."""
    upstream = upstream_ref()
    if upstream is None:
        print("ERROR: Cannot verify ahead/behind state without an upstream.", file=sys.stderr)
        return False

    remote = upstream.split("/", 1)[0]
    rf = run(["git", "fetch", remote], capture=True)
    if rf.returncode != 0:
        print(f"ERROR: git fetch for verification failed.\n{rf.stdout}\n{rf.stderr}", file=sys.stderr)
        return False

    r = run(["git", "rev-list", "--left-right", "--count", f"{upstream}...HEAD"], capture=True)
    if r.returncode != 0:
        print(f"ERROR: git ahead/behind verification failed.\n{r.stdout}\n{r.stderr}", file=sys.stderr)
        return False

    parts = (r.stdout or "").strip().split()
    if len(parts) != 2:
        print(f"ERROR: Unexpected ahead/behind output: {r.stdout!r}", file=sys.stderr)
        return False

    behind, ahead = (int(parts[0]), int(parts[1]))
    if ahead != 0:
        print(f"ERROR: Branch is still ahead of {upstream} by {ahead} commit(s).", file=sys.stderr)
        return False

    suffix = f" and behind by {behind}" if behind else ""
    print(f"\n-- git verification -----------------------------")
    print(f"  Verified: HEAD is not ahead of {upstream}{suffix}.")
    return True


def git_commit_push(commit_msg):
    """Stage all, commit, push with rebase-retry."""
    print("\n-- git add -A ------------------------------------")
    r = run_shell("git add -A")
    if r.returncode != 0:
        print("ERROR: git add failed.", file=sys.stderr)
        return False

    print("\n-- git commit ------------------------------------")
    r = run(["git", "commit", "-m", commit_msg], capture=True)
    if r.returncode != 0:
        combined = (r.stdout or "") + (r.stderr or "")
        if "nothing to commit" in combined:
            print("  (Nothing new to commit — skipping.)")
        else:
            print(f"ERROR: git commit failed.\n{r.stdout}\n{r.stderr}", file=sys.stderr)
            return False
    else:
        print(r.stdout, end="")

    if not push_current_branch():
        return False

    return verify_not_ahead()



def interactive_mode(machine_id):
    """Guide Claude through the finish-session steps interactively."""
    print(f"\n{'='*60}")
    print(f"  Finish Session — {machine_id}")
    print(f"{'='*60}")
    print(
        "\nThis script will:\n"
        "  1. Send your session letter to another agent\n"
        "  2. Commit all staged/unstaged changes\n"
        "  3. Push to origin\n"
    )

    to = input("Send letter to (agent name or 'all'): ").strip()
    subject = input("Subject: ").strip()
    print("Letter body (enter lines; blank line twice to finish):")
    lines = []
    while True:
        line = input()
        if line == "" and lines and lines[-1] == "":
            lines.pop()
            break
        lines.append(line)
    body = "\n".join(lines).strip()
    commit_msg = input("\nCommit message (one line): ").strip()

    return to, subject, body, commit_msg


def main():
    # Get machine identity for defaults
    machine_id_file = REPO_ROOT / ".machine-id"
    machine_id = machine_id_file.read_text().strip() if machine_id_file.exists() else "unknown"

    parser = argparse.ArgumentParser(
        description="End-of-session closer: send letter + commit + push"
    )
    parser.add_argument("--to",         help="Recipient (agent name or 'all')")
    parser.add_argument("--subject",    help="Letter subject line")
    parser.add_argument("--body",       help="Letter body (inline string)")
    parser.add_argument("--body-file",  help="Path to a file containing the letter body")
    parser.add_argument("--commit-msg", help="Git commit message")
    args = parser.parse_args()

    non_interactive = args.to and args.subject and (args.body or args.body_file) and args.commit_msg

    if non_interactive:
        to = args.to
        subject = args.subject
        commit_msg = args.commit_msg
        if args.body_file:
            body_source = f"file:{args.body_file}"
        else:
            body_source = args.body
    else:
        to, subject, body, commit_msg = interactive_mode(machine_id)
        body_source = body

    # Step 1: send message
    print(f"\n-- Sending letter to '{to}' ----------------------")
    if not send_message(to, subject, body_source, machine_id):
        sys.exit(2)

    # Step 2+3: commit and push
    if not git_commit_push(commit_msg):
        sys.exit(3)

    print(f"\n{'='*60}")
    print(f"  Session closed successfully.")
    print(f"  Letter delivered to: {to}")
    print(f"  Changes pushed to origin.")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
