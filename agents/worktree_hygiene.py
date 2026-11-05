#!/usr/bin/env python3
"""worktree_hygiene.py -- list, check and remove git worktrees of this repository.

Every worktree is a full 3.3 GB checkout; thirty of them filled the 228 GB data volume on
2026-11-04 and stopped every session on the machine. Rule (AGENTS.md, "Worktree hygiene"):
a session removes its own worktree at close-out once its branch is contained in origin/main,
and never creates one under /tmp or /private/tmp (periodic cleanup and reboots delete them,
losing untracked lane files).

    python3 agents/worktree_hygiene.py            # report every worktree (safe, read-only)
    python3 agents/worktree_hygiene.py --remove   # remove the clean ones whose HEAD is in origin/main
    python3 agents/worktree_hygiene.py --remove --path <worktree>   # just that one (same safety checks)

A worktree is removed only if (a) `git status --porcelain --untracked-files=all` is empty, or the
only modified file is a line-ending phantom (its diff is empty under --ignore-all-space), and
(b) its HEAD is an ancestor of origin/main (after `git fetch origin`). Dirty or unmerged worktrees
are reported with the reason and left alone: untracked files that were never pushed are somebody's
work. Exit status 1 if any worktree was kept for a reason that needs a human decision.
"""
import argparse, os, shutil, subprocess, sys

def git(args, cwd=None, check=True):
    p = subprocess.run(["git"] + args, cwd=cwd, capture_output=True, text=True)
    if check and p.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} failed in {cwd or os.getcwd()}: {p.stderr.strip()}")
    return p

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--remove", action="store_true", help="remove the clean, merged worktrees (default: report only)")
    ap.add_argument("--path", help="restrict to one worktree path")
    args = ap.parse_args()
    root = git(["rev-parse", "--show-toplevel"]).stdout.strip()
    main_wt = git(["worktree", "list", "--porcelain"], cwd=root).stdout.split("\n")[0].split(" ", 1)[1]
    git(["fetch", "-q", "origin"], cwd=root, check=False)
    git(["worktree", "prune"], cwd=root, check=False)
    entries = []
    cur = {}
    for line in git(["worktree", "list", "--porcelain"], cwd=root).stdout.splitlines() + [""]:
        if not line:
            if cur: entries.append(cur); cur = {}
            continue
        k, _, v = line.partition(" ")
        cur[k] = v
    kept_for_human = 0
    total_kb = 0
    for e in entries:
        wt = e.get("worktree")
        if wt == main_wt or (args.path and os.path.realpath(wt) != os.path.realpath(args.path)):
            continue
        if not os.path.isdir(wt):
            print(f"MISSING  {wt}  (pruned)"); continue
        head = e.get("HEAD", "?")
        branch = e.get("branch", "detached")
        du = subprocess.run(["du", "-sk", wt], capture_output=True, text=True).stdout.split()
        kb = int(du[0]) if du else 0
        total_kb += kb
        status = git(["status", "--porcelain", "--untracked-files=all"], cwd=wt, check=False).stdout.splitlines()
        phantom = False
        if status and all(s.startswith(" M") for s in status):
            phantom = git(["diff", "--ignore-all-space", "--quiet"], cwd=wt, check=False).returncode == 0
        merged = git(["merge-base", "--is-ancestor", head, "origin/main"], cwd=root, check=False).returncode == 0
        in_tmp = wt.startswith("/tmp/") or wt.startswith("/private/tmp/")
        reason = []
        if status and not phantom: reason.append(f"dirty({len(status)} files)")
        if not merged: reason.append("HEAD not in origin/main")
        tag = "OK" if not reason else "KEEP"
        print(f"{tag:5} {kb/1e6:5.1f} GB  {wt}  [{branch.replace('refs/heads/', '')}]{'  (under /tmp: at risk)' if in_tmp else ''}  {' '.join(reason)}{'  (phantom line-ending diff only)' if phantom else ''}")
        if reason:
            kept_for_human += 1
            continue
        if args.remove:
            p = git(["worktree", "remove"] + (["--force"] if phantom else []) + [wt], cwd=root, check=False)
            print("      removed" if p.returncode == 0 else f"      FAILED: {p.stderr.strip()}")
    git(["worktree", "prune"], cwd=root, check=False)
    free = shutil.disk_usage(root).free / 1e9
    print(f"worktrees total {total_kb/1e6:.1f} GB; free space on the volume {free:.1f} GB")
    if free < 20:
        print("WARNING: less than 20 GB free; remove merged worktrees before starting a session")
    sys.exit(1 if kept_for_human else 0)

if __name__ == "__main__":
    main()
