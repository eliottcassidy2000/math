#!/usr/bin/env python3
"""Print a bounded, read-only research startup packet.

This deliberately does not mutate message read-state, pull, rebase, reserve an
identifier, or rewrite a session file. It reports Git freshness, canonical
routes, recent commits, current correction headings, and topic-matched files so
an agent can decide what to read without scanning the entire corpus.
"""

from __future__ import annotations

import argparse
import heapq
import re
import subprocess
import sys
from pathlib import Path

from doc_surface import CANONICAL_ROUTES, TOPIC_SEARCH_ROUTES, TOPIC_SEARCH_TREES


REPO = Path(__file__).resolve().parent.parent
MAX_TOPIC_BYTES = 500
MAX_TOPIC_TERMS = 12


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


def git_state(recent: int) -> bool:
    section("Git state")
    branch = one_line("git", "branch", "--show-current") or "DETACHED"
    head = one_line("git", "rev-parse", "--short=12", "HEAD") or "unknown"
    upstream = one_line(
        "git", "rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}"
    )
    live_main = (
        "origin/main"
        if succeeds("git", "show-ref", "--verify", "refs/remotes/origin/main")
        else ""
    )
    target = upstream or live_main
    print(f"branch: {branch}")
    print(f"HEAD:   {head}")
    print(f"target (cached refs): {target or 'none configured'}")
    healthy = bool(target)

    if target:
        counts = one_line(
            "git", "rev-list", "--left-right", "--count", f"{target}...HEAD"
        ).split()
        if len(counts) == 2:
            behind, ahead = counts
            print(f"relative to target: behind {behind}, ahead {ahead}")
            if behind != "0":
                print("ACTION: sync from a clean tree before relying on frontier state.")
                healthy = False

    if live_main and target != live_main:
        counts = one_line(
            "git", "rev-list", "--left-right", "--count", f"{live_main}...HEAD"
        ).split()
        if len(counts) == 2:
            behind, ahead = counts
            print(f"relative to live origin/main: behind {behind}, ahead {ahead}")
            if behind != "0":
                print("ACTION: incorporate the live shared surface before research.")
                healthy = False

    status = one_line("git", "status", "--short")
    if status:
        rows = status.splitlines()
        print(f"working tree: DIRTY ({len(rows)} paths; first 12 shown)")
        for row in rows[:12]:
            print(f"  {row}")
        if len(rows) > 12:
            print(f"  ... {len(rows) - 12} more")
        print("ACTION: identify ownership; do not sweep unrelated files into a commit.")
        healthy = False
    else:
        print("working tree: clean")

    section(f"Recent commits ({recent})")
    log = one_line(
        "git", "log", f"-{recent}", "--oneline", "--decorate", "--no-merges"
    )
    print(log or "(no Git history available)")
    print("Freshness is cached-only: run `git fetch origin` before strict reliance.")
    return healthy


def route_packet() -> bool:
    section("Canonical startup route")
    healthy = True
    for path in CANONICAL_ROUTES:
        marker = "OK" if (REPO / path).exists() else "MISSING"
        print(f"[{marker}] {path}")
        healthy = healthy and marker == "OK"
    print("Historical logs/reflections are searched for provenance, not read as truth.")
    return healthy


STOPWORDS = {
    "agent", "agents", "audit", "current", "documentation", "exact",
    "explore", "make", "math", "mathematical", "problem", "proof",
    "repo", "repository", "research", "result", "results", "session",
    "theorem", "update", "work", "working",
}


def topic_terms(topic: str) -> list[tuple[str, int]]:
    """Return specificity-weighted literal terms, with identifiers first."""
    if not topic.strip():
        raise ValueError("topic must contain a specific mathematical object or statement")
    if len(topic.encode("utf-8")) > MAX_TOPIC_BYTES:
        raise ValueError(f"topic exceeds the {MAX_TOPIC_BYTES}-byte safety limit")
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
    if not deduplicated:
        raise ValueError("topic contains only generic/stop words; add an object, constant, or ID")
    return deduplicated[:MAX_TOPIC_TERMS]


def match_score(text: str, terms: list[tuple[str, int]]) -> tuple[int, int]:
    folded = text.casefold()
    matches = [(term, weight) for term, weight in terms if term.casefold() in folded]
    return sum(weight for _, weight in matches), len(matches)


def topic_hits(
    topic: str, max_matches: int, terms: list[tuple[str, int]]
) -> None:
    section(f"Topic packet: {topic}")
    has_identifier = any(weight >= 100 for _, weight in terms)
    required_matches = 1 if has_identifier or len(terms) == 1 else 2
    print("rank terms: " + ", ".join(term for term, _ in terms[:8]))

    current_heap: list[tuple[int, int, str]] = []
    current_count = 0
    for route_rank, relative in enumerate(TOPIC_SEARCH_ROUTES):
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
                current_count += 1
                row = f"{relative}:{line_number}:{line}"
                item = (score, -route_rank, row)
                if len(current_heap) < max_matches:
                    heapq.heappush(current_heap, item)
                elif item > current_heap[0]:
                    heapq.heapreplace(current_heap, item)
    current_rows = [(score, -neg_rank, row) for score, neg_rank, row in current_heap]
    current_rows.sort(key=lambda row: (-row[0], row[1], row[2]))

    print("Current-route matches (ranked):")
    if current_rows:
        for _, _, row in current_rows[:max_matches]:
            print(f"  {row}")
        if current_count > len(current_rows):
            print(f"  ... {current_count - len(current_rows)} more qualifying route matches")
    else:
        print("  (none; try exact constants, theorem IDs, quantifiers, or synonyms)")

    file_scores: dict[str, tuple[int, int]] = {}
    identifier_paths: set[str] = set()
    for term, weight in terms:
        files = run(
            "rg", "--files-with-matches", "--fixed-strings", "--ignore-case",
            "--", term, *TOPIC_SEARCH_TREES
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
            row[1].endswith("/INDEX.md"),
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
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("## MISTAKE-"):
                headings.append(line[3:].rstrip())
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
    parser.add_argument(
        "--strict", action="store_true",
        help="Exit nonzero when cached Git state is dirty/behind or a canonical route is missing",
    )
    args = parser.parse_args()
    if not 1 <= args.recent <= 20 or not 1 <= args.max_matches <= 30:
        parser.error("--recent must be 1..20 and --max-matches must be 1..30")

    try:
        terms = topic_terms(args.topic)
    except ValueError as exc:
        parser.error(str(exc))

    print("MATH RESEARCH STARTUP PACKET (read-only, bounded)")
    git_ok = git_state(args.recent)
    routes_ok = route_packet()
    recent_guardrails(8)
    topic_hits(args.topic, args.max_matches, terms)
    identity_note()
    if args.strict and not (git_ok and routes_ok):
        print("\nSTRICT CHECK FAILED: resolve the actions above and rerun.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
