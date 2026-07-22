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

from doc_surface import (
    ALWAYS_READ_ROUTES,
    MAX_EMITTED_MATCH_BYTES,
    ON_DEMAND_ROUTES,
    SEARCHABLE_PREFIX_END,
    TOPIC_ARTIFACT_TREES,
    TOPIC_SEARCH_ROUTES,
    TOPIC_SEARCH_TREES,
    WORKFLOW_ROUTES,
)


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
    tiers = (
        ("always read", ALWAYS_READ_ROUTES),
        ("mathematical workflow", WORKFLOW_ROUTES),
        ("on demand", ON_DEMAND_ROUTES),
    )
    for label, routes in tiers:
        print(f"{label}:")
        for path in routes:
            marker = "OK" if (REPO / path).exists() else "MISSING"
            print(f"  [{marker}] {path}")
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
    if any(ord(character) < 32 or ord(character) == 127 for character in topic):
        raise ValueError("topic must not contain control characters")
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


def truncate_utf8(value: str, limit: int = MAX_EMITTED_MATCH_BYTES) -> str:
    """Collapse whitespace and truncate without splitting a UTF-8 codepoint."""
    compact = " ".join(value.split())
    encoded = compact.encode("utf-8")
    if len(encoded) <= limit:
        return compact
    clipped = encoded[: max(0, limit - 3)]
    while True:
        try:
            return clipped.decode("utf-8") + "..."
        except UnicodeDecodeError:
            clipped = clipped[:-1]


def maintained_route_lines(relative: str) -> list[str]:
    text = (REPO / relative).read_text(encoding="utf-8", errors="replace")
    marker = SEARCHABLE_PREFIX_END.get(relative)
    if marker and marker in text:
        text = text.split(marker, 1)[0]
    return text.splitlines()


STATUS_WORDS = (
    "REFUTED", "SUPERSEDED", "PROVED", "FINITE-EXACT", "CITED",
    "CONDITIONAL", "VERIFIED", "PARTIAL", "CLAIMED", "OPEN", "MIXED",
)


def file_status(relative: str) -> str:
    if relative.startswith("07-reflections/"):
        return "HISTORY"
    path = REPO / relative
    try:
        prefix = path.read_text(encoding="utf-8", errors="replace")[:5000]
    except OSError:
        return "UNKNOWN"
    match = re.search(r"^status:\s*(?:>\s*)?(.*)$", prefix, flags=re.MULTILINE)
    probe = match.group(1) if match else prefix[:400]
    if match and not probe.strip():
        probe = prefix[match.end():match.end() + 400]
    upper = probe.upper()
    return next((word for word in STATUS_WORDS if word in upper), "UNLABELLED")


def scored_files(
    terms: list[tuple[str, int]], trees: tuple[str, ...]
) -> tuple[dict[str, tuple[int, int]], set[str]]:
    scores: dict[str, tuple[int, int]] = {}
    identifier_paths: set[str] = set()
    for term, weight in terms:
        files = run(
            "rg", "--files-with-matches", "--fixed-strings", "--ignore-case",
            "--", term, *trees
        )
        if files.returncode not in (0, 1):
            continue
        for relative in files.stdout.splitlines():
            score, count = scores.get(relative, (0, 0))
            filename_bonus = weight if term.casefold() in Path(relative).name.casefold() else 0
            scores[relative] = (score + weight + filename_bonus, count + 1)
            if weight >= 100:
                identifier_paths.add(relative)
    return scores, identifier_paths


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
        for line_number, line in enumerate(maintained_route_lines(relative), 1):
            score, count = match_score(line, terms)
            identifier_hit = any(
                weight >= 100 and term.casefold() in line.casefold()
                for term, weight in terms
            )
            if count >= required_matches and (not has_identifier or identifier_hit):
                current_count += 1
                row = truncate_utf8(f"{relative}:{line_number}: {line}")
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

    file_scores, identifier_paths = scored_files(terms, TOPIC_SEARCH_TREES)
    file_rows = [
        (score, path)
        for path, (score, count) in file_scores.items()
        if count >= required_matches and (not has_identifier or path in identifier_paths)
    ]
    groups = (
        ("Canon", "01-canon/theorems/"),
        ("Hypotheses", "05-knowledge/hypotheses/"),
        ("Historical reflections", "07-reflections/"),
    )
    per_group = max(1, max_matches // len(groups))
    print("Detailed files (truth precedence first; specificity within group):")
    for label, prefix in groups:
        rows = sorted(
            (row for row in file_rows if row[1].startswith(prefix)),
            key=lambda row: (row[1].endswith("/INDEX.md"), -row[0], row[1]),
        )
        print(f"  {label}:")
        if not rows:
            print("    (none)")
            continue
        for _, relative in rows[:per_group]:
            print(f"    [{file_status(relative)}] {relative}")
        if len(rows) > per_group:
            print(f"    ... {len(rows) - per_group} more; refine the query")

    artifact_scores, artifact_ids = scored_files(terms, TOPIC_ARTIFACT_TREES)
    artifact_rows = [
        (score, relative)
        for relative, (score, count) in artifact_scores.items()
        if count >= required_matches and (not has_identifier or relative in artifact_ids)
    ]
    artifact_rows.sort(key=lambda row: (-row[0], row[1]))
    artifact_limit = max(2, min(6, max_matches // 2))
    print("Niche/artifact hits (evidence or inspiration; not truth authority):")
    if artifact_rows:
        for _, relative in artifact_rows[:artifact_limit]:
            print(f"  {relative}")
        if len(artifact_rows) > artifact_limit:
            print(f"  ... {len(artifact_rows) - artifact_limit} more; refine the query")
    else:
        print("  (none)")


def recent_guardrails(limit: int, terms: list[tuple[str, int]]) -> None:
    section("Front-of-ledger correction headings")
    path = REPO / "01-canon/MISTAKES.md"
    if not path.exists():
        print("(missing mistakes ledger)")
        return
    text = path.read_text(encoding="utf-8", errors="replace")
    entries = list(
        re.finditer(r"^## (MISTAKE-\d+[^\n]*)\n(.*?)(?=^## MISTAKE-|\Z)", text, re.M | re.S)
    )
    headings = [match.group(1) for match in entries[:limit]]
    for heading in headings:
        print(f"- {heading}")
    matched: list[tuple[int, str]] = []
    required_matches = 1 if any(weight >= 100 for _, weight in terms) or len(terms) == 1 else 2
    for match in entries:
        score, count = match_score(match.group(0), terms)
        if count >= required_matches:
            matched.append((score, match.group(1)))
    matched.sort(key=lambda row: (-row[0], row[1]))
    newest = set(headings)
    relevant = [heading for _, heading in matched if heading not in newest][:4]
    print("Topic-matched older corrections:")
    if relevant:
        for heading in relevant:
            print(f"- {heading}")
    else:
        print("- (none beyond the recent headings)")
    print("Read ACTIVE-GUARDRAILS first; open the full entry only when relevant.")


def session_posture(topic: str) -> None:
    section("Session posture")
    print(f"Anchor: {truncate_utf8(topic, 220)}")
    print("Choose one underexplored Niche and one freely generated Wildcard.")
    print("Keep a 3–7 concept board; compare each new result against every item.")
    print("Explain the mechanism or failure anatomy, not only the verdict.")
    print("Type connections as map / preserved predicate / loss / sidecar / test.")


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
    session_posture(args.topic)
    recent_guardrails(8, terms)
    topic_hits(args.topic, args.max_matches, terms)
    identity_note()
    if args.strict and not (git_ok and routes_ok):
        print("\nSTRICT CHECK FAILED: resolve the actions above and rerun.")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
