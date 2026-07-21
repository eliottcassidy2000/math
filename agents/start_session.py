#!/usr/bin/env python3
"""Print a bounded, read-only research startup packet.

This deliberately does not mutate message read-state, pull, rebase, reserve an
identifier, or rewrite a session file. It reports Git freshness, canonical
routes, recent commits, current correction headings, and topic-matched files so
an agent can decide what to read without scanning the entire corpus.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parent.parent

ROUTES = (
    "00-navigation/START-HERE.md",
    "00-navigation/CURRENT-FRONTIER.md",
    "01-canon/ACTIVE-GUARDRAILS.md",
    "00-navigation/RESEARCH-PROTOCOL.md",
    "00-navigation/META-PATTERNS.md",
    "05-knowledge/reference/CORE-PAPERS.md",
    "00-navigation/PROBLEM-LEDGER.md",
)

SEARCH_ROOTS = (
    "00-navigation/CURRENT-FRONTIER.md",
    "01-canon/ACTIVE-GUARDRAILS.md",
    "05-knowledge/reference/CORE-PAPERS.md",
    "00-navigation/PROBLEM-LEDGER.md",
    "00-navigation/LRC14-PROOF-MAP.md",
    "01-canon/theorems",
    "05-knowledge/hypotheses",
    "07-reflections",
)


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )


def one_line(*args: str) -> str:
    result = run(*args)
    return result.stdout.strip() if result.returncode == 0 else ""


def succeeds(*args: str) -> bool:
    return run(*args).returncode == 0


def section(title: str) -> None:
    print(f"\n== {title} ==")


def git_state(recent: int) -> None:
    section("Git state")
    branch = one_line("git", "branch", "--show-current") or "DETACHED"
    head = one_line("git", "rev-parse", "--short=12", "HEAD") or "unknown"
    upstream = one_line(
        "git", "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"
    )
    target = upstream or (
        "origin/main"
        if succeeds("git", "show-ref", "--verify", "refs/remotes/origin/main")
        else ""
    )
    print(f"branch: {branch}")
    print(f"HEAD:   {head}")
    print(f"target: {target or 'none configured'}")

    if target:
        counts = one_line(
            "git", "rev-list", "--left-right", "--count", f"{target}...HEAD"
        ).split()
        if len(counts) == 2:
            behind, ahead = counts
            print(f"relative to target: behind {behind}, ahead {ahead}")
            if behind != "0":
                print("ACTION: sync from a clean tree before relying on frontier state.")

    status = one_line("git", "status", "--short")
    if status:
        rows = status.splitlines()
        print(f"working tree: DIRTY ({len(rows)} paths; first 12 shown)")
        for row in rows[:12]:
            print(f"  {row}")
        if len(rows) > 12:
            print(f"  ... {len(rows) - 12} more")
        print("ACTION: identify ownership; do not sweep unrelated files into a commit.")
    else:
        print("working tree: clean")

    section(f"Recent commits ({recent})")
    log = one_line(
        "git", "log", f"-{recent}", "--oneline", "--decorate", "--no-merges"
    )
    print(log or "(no Git history available)")


def route_packet() -> None:
    section("Canonical startup route")
    for path in ROUTES:
        marker = "OK" if (REPO / path).exists() else "MISSING"
        print(f"[{marker}] {path}")
    print("Historical logs/reflections are searched for provenance, not read as truth.")


def topic_pattern(topic: str) -> str:
    tokens = re.findall(r"[A-Za-z0-9_.()+/-]+", topic)
    useful = [token for token in tokens if len(token) >= 3]
    terms = [topic.strip(), *useful]
    seen: set[str] = set()
    escaped: list[str] = []
    for term in terms:
        key = term.casefold()
        if term and key not in seen:
            seen.add(key)
            escaped.append(re.escape(term))
    return "|".join(escaped) or re.escape(topic)


def topic_hits(topic: str, max_matches: int) -> None:
    section(f"Topic packet: {topic}")
    pattern = topic_pattern(topic)

    current = run(
        "rg",
        "--line-number",
        "--ignore-case",
        "--max-count",
        "3",
        "--",
        pattern,
        *SEARCH_ROOTS[:5],
    )
    current_rows = current.stdout.splitlines() if current.returncode in (0, 1) else []
    print("Current-route matches:")
    if current_rows:
        for row in current_rows[:max_matches]:
            print(f"  {row}")
        if len(current_rows) > max_matches:
            print(f"  ... {len(current_rows) - max_matches} more route matches")
    else:
        print("  (none; try exact constants, theorem IDs, quantifiers, or synonyms)")

    files = run(
        "rg",
        "--files-with-matches",
        "--ignore-case",
        "--",
        pattern,
        *SEARCH_ROOTS[5:],
    )
    file_rows = files.stdout.splitlines() if files.returncode in (0, 1) else []
    priorities = {"01-canon/theorems": 0, "05-knowledge/hypotheses": 1, "07-reflections": 2}
    file_rows.sort(key=lambda p: (next((v for k, v in priorities.items() if p.startswith(k)), 9), p))
    print("Detailed files (canon before hypotheses before reflections):")
    if file_rows:
        for path in file_rows[:max_matches]:
            print(f"  {path}")
        if len(file_rows) > max_matches:
            print(f"  ... {len(file_rows) - max_matches} more; refine the query")
    else:
        print("  (none)")


def recent_guardrails(limit: int) -> None:
    section("Front-of-ledger correction headings")
    path = REPO / "01-canon/MISTAKES.md"
    if not path.exists():
        print("(missing mistakes ledger)")
        return
    headings: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("## MISTAKE-"):
            headings.append(line[3:])
            if len(headings) == limit:
                break
    for heading in headings:
        print(f"- {heading}")
    print("Read ACTIVE-GUARDRAILS first; open the full entry only when relevant.")


def identity_note() -> None:
    section("Coordination")
    machine = REPO / ".machine-id"
    if machine.exists():
        print(f"machine identity: {machine.read_text().strip()}")
        print("Use `python3 agents/processor.py --status --peek --limit 8`;")
        print("remove --peek only when intentionally consuming the displayed messages.")
    else:
        print("No .machine-id: use a temporary session name; do not invent registration.")
    print("Legacy IDs collide: cite every theorem as ID + slug/path.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topic", required=True, help="Target statement, object, or theorem ID")
    parser.add_argument("--recent", type=int, default=8, help="Recent commits to print (1-20)")
    parser.add_argument(
        "--max-matches", type=int, default=12, help="Maximum rows/files per topic group (1-30)"
    )
    args = parser.parse_args()
    if not 1 <= args.recent <= 20 or not 1 <= args.max_matches <= 30:
        parser.error("--recent must be 1..20 and --max-matches must be 1..30")

    print("MATH RESEARCH STARTUP PACKET (read-only, bounded)")
    git_state(args.recent)
    route_packet()
    recent_guardrails(8)
    topic_hits(args.topic, args.max_matches)
    identity_note()
    return 0


if __name__ == "__main__":
    sys.exit(main())
