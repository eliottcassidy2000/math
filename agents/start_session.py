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


STOPWORDS = {
    "agent", "agents", "audit", "current", "documentation", "exact",
    "explore", "make", "math", "mathematical", "problem", "proof",
    "repo", "repository", "research", "result", "results", "session",
    "theorem", "update", "work", "working",
}


def topic_terms(topic: str) -> list[tuple[str, int]]:
    """Return specificity-weighted literal terms, with identifiers first."""
    identifiers = [
        f"{match.group(1).upper()}-{match.group(2)}"
        for match in re.finditer(
            r"\b(THM|HYP|MISTAKE|LEM|OPEN-Q)-?(\d+)\b",
            topic,
            flags=re.IGNORECASE,
        )
    ]
    tokens = re.findall(r"[A-Za-z0-9][A-Za-z0-9_.()+/-]*", topic)
    useful = [
        token
        for token in tokens
        if len(token) >= 3 and token.casefold() not in STOPWORDS
    ]

    weighted: list[tuple[str, int]] = []
    if len(useful) >= 2:
        weighted.append((" ".join(useful), 30))
    for identifier in identifiers:
        weighted.append((identifier, 100))
    for token in useful:
        structured = bool(re.search(r"[A-Za-z].*\d|\d.*[A-Za-z]", token))
        weighted.append((token, 50 if structured else min(14, len(token))))

    seen: set[str] = set()
    deduplicated: list[tuple[str, int]] = []
    for term, weight in weighted:
        key = term.casefold()
        if key and key not in seen:
            seen.add(key)
            deduplicated.append((term, weight))
    return deduplicated or [(topic.strip(), 1)]


def match_score(text: str, terms: list[tuple[str, int]]) -> tuple[int, int]:
    folded = text.casefold()
    matches = [(term, weight) for term, weight in terms if term.casefold() in folded]
    return sum(weight for _, weight in matches), len(matches)


def topic_hits(topic: str, max_matches: int) -> None:
    section(f"Topic packet: {topic}")
    terms = topic_terms(topic)
    has_identifier = any(weight >= 100 for _, weight in terms)
    required_matches = 1 if has_identifier or len(terms) == 1 else 2
    print("rank terms: " + ", ".join(term for term, _ in terms[:8]))

    current_rows: list[tuple[int, int, str]] = []
    for route_rank, relative in enumerate(SEARCH_ROOTS[:5]):
        path = REPO / relative
        if not path.is_file():
            continue
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8", errors="replace").splitlines(), 1
        ):
            score, count = match_score(line, terms)
            identifier_hit = any(
                weight >= 100 and term.casefold() in line.casefold()
                for term, weight in terms
            )
            if count >= required_matches and (not has_identifier or identifier_hit):
                current_rows.append((score, route_rank, f"{relative}:{line_number}:{line}"))
    current_rows.sort(key=lambda row: (-row[0], row[1], row[2]))

    print("Current-route matches (ranked):")
    if current_rows:
        for _, _, row in current_rows[:max_matches]:
            print(f"  {row}")
        if len(current_rows) > max_matches:
            print(f"  ... {len(current_rows) - max_matches} more qualifying route matches")
    else:
        print("  (none; try exact constants, theorem IDs, quantifiers, or synonyms)")

    file_scores: dict[str, tuple[int, int]] = {}
    identifier_paths: set[str] = set()
    for term, weight in terms:
        files = run(
            "rg", "--files-with-matches", "--fixed-strings", "--ignore-case",
            "--", term, *SEARCH_ROOTS[5:]
        )
        if files.returncode not in (0, 1):
            continue
        for path in files.stdout.splitlines():
            score, count = file_scores.get(path, (0, 0))
            filename_bonus = weight if term.casefold() in Path(path).name.casefold() else 0
            file_scores[path] = (score + weight + filename_bonus, count + 1)
            if weight >= 100:
                identifier_paths.add(path)

    priorities = {"01-canon/theorems": 0, "05-knowledge/hypotheses": 1, "07-reflections": 2}
    file_rows = [
        (score, path)
        for path, (score, count) in file_scores.items()
        if count >= required_matches and (not has_identifier or path in identifier_paths)
    ]
    file_rows.sort(
        key=lambda row: (
            -row[0],
            next((v for k, v in priorities.items() if row[1].startswith(k)), 9),
            row[1],
        )
    )
    print("Detailed files (ranked by specificity; canon wins ties):")
    if file_rows:
        for _, path in file_rows[:max_matches]:
            print(f"  {path}")
        if len(file_rows) > max_matches:
            print(f"  ... {len(file_rows) - max_matches} more qualifying files; refine the query")
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
