#!/usr/bin/env python3
"""Validate the small maintained surface that future research agents read.

The repository intentionally keeps thousands of historical Markdown files.
This checker does not pretend they are all current. It protects the bounded
startup route, its local links, runtime-state policy, and a few headline truth
sentinels that would be costly for a future agent to get wrong.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from urllib.parse import unquote


REPO = Path(__file__).resolve().parent.parent

STARTUP_DOCS = (
    "AGENTS.md",
    "CLAUDE.md",
    "README.md",
    "00-navigation/START-HERE.md",
    "00-navigation/CURRENT-FRONTIER.md",
    "00-navigation/RESEARCH-PROTOCOL.md",
    "00-navigation/META-PATTERNS.md",
    "00-navigation/CONCURRENT-SESSIONS.md",
    "01-canon/ACTIVE-GUARDRAILS.md",
    "01-canon/README.md",
    "05-knowledge/README.md",
    "05-knowledge/reference/CORE-PAPERS.md",
    "07-reflections/README.md",
    "agents/README.md",
)

LINE_BUDGETS = {
    "AGENTS.md": 140,
    "CLAUDE.md": 70,
    "README.md": 120,
    "00-navigation/START-HERE.md": 150,
    "00-navigation/CURRENT-FRONTIER.md": 450,
    "01-canon/ACTIVE-GUARDRAILS.md": 180,
    "agents/README.md": 140,
}

HEADLINE_SENTINELS = {
    "00-navigation/START-HERE.md": (
        "14 total runners",
        "is **OPEN**",
        "uniform good period `q <= 25` is false",
        "uniform emptiness of the twelve-speed sporadic tight branch remains open",
        "unrestricted GMC(2)",
        "**PROVED in repo",
    ),
    "00-navigation/CURRENT-FRONTIER.md": (
        "## LRC(14)",
        "**OPEN.**",
        "first good period `q=27`",
        "disc_v>=6|G'_{~v}|^2",
        "proves NC2 and hence unrestricted GMC(2)",
        "GMC is false for every dimension at least 3",
    ),
    "01-canon/ACTIVE-GUARDRAILS.md": (
        "No uniform `q <= 25` good-period theorem",
        "Uniform twelve-speed sporadic emptiness is OPEN",
        "HYP-8815 is a heuristic, not a disproof characterization",
        "NC2/GMC(2) is proved, not fully formalized",
    ),
}

LINK_RE = re.compile(r"!?(?:\[[^\]]*\])\(([^)\n]+)\)")


def git(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ("git", *args), cwd=REPO, text=True, capture_output=True, check=False
    )


def local_link_target(raw: str) -> str | None:
    target = raw.strip()
    if target.startswith("<") and target.endswith(">"):
        target = target[1:-1]
    if target.startswith(("http://", "https://", "mailto:", "#")):
        return None
    # Markdown titles are not used on the maintained surface. Flag rather than
    # guessing if one ever appears and breaks this simple parser.
    target = unquote(target.split("#", 1)[0])
    return target or None


def main() -> int:
    errors: list[str] = []

    for relative in STARTUP_DOCS:
        path = REPO / relative
        if not path.is_file():
            errors.append(f"missing maintained startup document: {relative}")
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        budget = LINE_BUDGETS.get(relative)
        if budget is not None and len(lines) > budget:
            errors.append(f"{relative}: {len(lines)} lines exceeds startup budget {budget}")

        for raw in LINK_RE.findall(text):
            target = local_link_target(raw)
            if target is None:
                continue
            resolved = Path(target) if Path(target).is_absolute() else path.parent / target
            if not resolved.exists():
                errors.append(f"{relative}: broken local link {raw!r}")

    for relative, sentinels in HEADLINE_SENTINELS.items():
        path = REPO / relative
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for sentinel in sentinels:
            if sentinel not in text:
                errors.append(f"{relative}: missing headline truth sentinel {sentinel!r}")

    start = (REPO / "00-navigation/START-HERE.md").read_text(
        encoding="utf-8", errors="replace"
    )
    if re.search(r"refreshed[^\n]*against\s+`[0-9a-f]{7,40}`", start, re.I):
        errors.append("START-HERE.md: static commit pin will become stale; startup script prints HEAD")

    for relative in (".codex/hooks.json", ".claude/settings.local.json"):
        path = REPO / relative
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"{relative}: invalid or unreadable JSON: {exc}")
            continue
        serialized = json.dumps(payload)
        if "check_session_closed.py" in serialized or '"Stop"' in serialized:
            errors.append(f"{relative}: stateful Stop hook has returned")

    tracked_state = git("ls-files", "--", "agents/.session-state.json").stdout.strip()
    if tracked_state:
        errors.append("agents/.session-state.json is tracked; runtime state must stay local")

    for relative in (
        "agents/start_session.py",
        "agents/check_docs.py",
        "agents/new_session_worktree.sh",
    ):
        if not os.access(REPO / relative, os.X_OK):
            errors.append(f"{relative}: expected executable bit is missing")

    mistakes = (REPO / "01-canon/MISTAKES.md").read_text(
        encoding="utf-8", errors="replace"
    )
    recent_ids = [
        value
        for value in re.findall(r"^## MISTAKE-(\d+)\b", mistakes, re.MULTILINE)
        if int(value) >= 209
    ]
    for value, count in sorted(Counter(recent_ids).items(), key=lambda item: int(item[0])):
        if count > 1:
            errors.append(f"MISTAKES.md: current-range MISTAKE-{value} occurs {count} times")

    if errors:
        print("Agent-facing documentation check FAILED:")
        for error in errors:
            print(f"- {error}")
        return 1

    print(
        f"Agent-facing documentation check passed: {len(STARTUP_DOCS)} maintained "
        "documents, local links, headline sentinels, configs, runtime policy, and IDs."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
