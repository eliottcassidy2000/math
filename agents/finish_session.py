#!/usr/bin/env python3
"""
finish_session.py — End-of-session closer for Claude agents.

Performs all three mandatory end-of-session steps in one command:
  1. Sends a summary message to another agent (or broadcast)
  2. git add -- <your paths> (+ the letter files it just wrote) && git commit
  3. git push to the current branch's upstream (with rebase-retry on conflict)

Owner directive 2026-11-13: sessions work directly on `main` in the checkout
they are given, never in a new worktree. That checkout is SHARED with other
live sessions, so pass --paths with exactly the files this session wrote;
without --paths the script refuses to stage anything in the main checkout
(it falls back to `git add -A` only inside a linked worktree).

Usage (interactive):
    python3 agents/finish_session.py

Usage (non-interactive, for Claude):
    python3 agents/finish_session.py \\
        --to all \\
        --subject "Session summary: <one-line>" \\
        --body-file /tmp/letter.md \\
        --commit-msg "kind-pasteur-2026-03-05-S2: <one-line summary>" \\
        --paths 05-knowledge/results/my_note.md 01-canon/theorems/THM-XXXX-slug.md

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


def in_linked_worktree():
    """True inside a linked worktree, False in the main checkout."""
    a = run(["git", "rev-parse", "--git-dir"], capture=True)
    b = run(["git", "rev-parse", "--git-common-dir"], capture=True)
    if a.returncode != 0 or b.returncode != 0:
        return False
    gd = (REPO_ROOT / a.stdout.strip()).resolve()
    cd = (REPO_ROOT / b.stdout.strip()).resolve()
    return gd != cd


def dirty_paths():
    """Set of untracked (not ignored) and modified tracked paths right now."""
    out = set()
    r = run(["git", "ls-files", "--others", "--exclude-standard"], capture=True)
    out.update(x for x in (r.stdout or "").splitlines() if x)
    r = run(["git", "diff", "--name-only"], capture=True)
    out.update(x for x in (r.stdout or "").splitlines() if x)
    return out


def sparse_flag():
    """`git add --sparse` lets a sparse checkout stage paths outside its cone."""
    r = run(["git", "config", "--get", "core.sparseCheckout"], capture=True)
    return ["--sparse"] if (r.stdout or "").strip() == "true" else []


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

        # --autostash: the shared main checkout normally holds other live
        # sessions' uncommitted edits; they are set aside and restored.
        rr = run(["git", "rebase", "--autostash", f"{remote}/{remote_branch}"], capture=True)
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


def git_commit_push(commit_msg, paths=None):
    """Stage this session's paths, commit only them, push with rebase-retry.

    paths=None keeps the old `git add -A` behaviour, but only inside a linked
    worktree; in the shared main checkout it refuses rather than sweep other
    sessions' files into this commit."""
    if paths is None:
        if not in_linked_worktree():
            print("ERROR: this is the shared main checkout; pass --paths with only "
                  "the files this session wrote (no `git add -A` here).", file=sys.stderr)
            return False
        print("\n-- git add -A (linked worktree) -------------------")
        r = run_shell("git add -A")
        if r.returncode != 0:
            print("ERROR: git add failed.", file=sys.stderr)
            return False
        commit_cmd = ["git", "commit", "-m", commit_msg]
    else:
        paths = sorted(set(paths))
        if not paths:
            print("  (No paths given — nothing to commit.)")
            return push_current_branch() and verify_not_ahead()
        print("\n-- git add -- <paths> -----------------------------")
        for x in paths:
            print(f"  {x}")
        r = run(["git", "add", *sparse_flag(), "--", *paths], capture=True)
        if r.returncode != 0:
            print(f"ERROR: git add failed.\n{r.stdout}\n{r.stderr}", file=sys.stderr)
            return False
        # Pathspec commit: other sessions' staged entries are never swept in.
        commit_cmd = ["git", "commit", "-m", commit_msg, "--", *paths]

    print("\n-- git commit ------------------------------------")
    r = run(commit_cmd, capture=True)
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
    parser.add_argument("--paths", nargs="*", default=None,
                        help="Exactly the files this session wrote (required in the "
                             "shared main checkout). The letter files are added automatically.")
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

    if args.paths is None and not in_linked_worktree():
        print("ERROR: shared main checkout: pass --paths <files this session wrote>.",
              file=sys.stderr)
        sys.exit(1)

    # Step 1: send message (remember which files the letter creates)
    before = dirty_paths()
    print(f"\n-- Sending letter to '{to}' ----------------------")
    if not send_message(to, subject, body_source, machine_id):
        sys.exit(2)
    letter_files = sorted(dirty_paths() - before)

    # Step 2+3: commit and push
    paths = None if args.paths is None else list(args.paths) + letter_files
    if not git_commit_push(commit_msg, paths):
        sys.exit(3)

    print(f"\n{'='*60}")
    print(f"  Session closed successfully.")
    print(f"  Letter delivered to: {to}")
    print(f"  Changes pushed to origin.")
    try:
        import subprocess as _sp
        _top = _sp.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True).stdout.strip()
        _main = _sp.run(["git", "worktree", "list", "--porcelain"], capture_output=True, text=True).stdout.split("\n")[0].split(" ", 1)[1]
        if _top and _main and os.path.realpath(_top) != os.path.realpath(_main):
            print(f"  This is a session worktree ({_top}). Every worktree is a full checkout:")
            print(f"  remove it now, from the main checkout, after leaving it:")
            print(f"      cd {_main} && python3 agents/worktree_hygiene.py --remove --path {_top}")
    except Exception:
        pass
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
