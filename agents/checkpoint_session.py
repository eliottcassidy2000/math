#!/usr/bin/env python3
"""
checkpoint_session.py - lightweight mid-session checkpoint for concurrent agents.

Stages intentional repo changes, commits if there is anything to commit, pushes
the current branch to its upstream, and verifies the branch is no longer ahead
of upstream. Unlike finish_session.py, this does not send a session letter.
"""

import argparse
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from finish_session import git_commit_push, run  # noqa: E402


def print_status():
    result = run(["git", "status", "--short", "--branch"], capture=True)
    if result.returncode == 0:
        print(result.stdout, end="")
    else:
        print(result.stderr, file=sys.stderr, end="")


def main():
    parser = argparse.ArgumentParser(
        description="Commit and push a lightweight concurrent-session checkpoint."
    )
    parser.add_argument(
        "-m",
        "--message",
        required=True,
        help="Checkpoint commit message, usually prefixed with [instance-id].",
    )
    parser.add_argument(
        "--no-status",
        action="store_true",
        help="Skip printing git status before the checkpoint.",
    )
    args = parser.parse_args()

    if not args.no_status:
        print("-- git status -----------------------------------")
        print_status()

    if not git_commit_push(args.message):
        sys.exit(4)

    print("\nCheckpoint pushed and verified.")


if __name__ == "__main__":
    main()
