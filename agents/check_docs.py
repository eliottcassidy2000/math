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

from doc_surface import (
    LINE_BUDGETS,
    MAX_STARTUP_LINE_BYTES,
    STARTUP_BYTE_BUDGET,
    STARTUP_DOCS,
)


REPO = Path(__file__).resolve().parent.parent

HEADLINE_SENTINELS = {
    "README.md": (
        "LRC(14) is OPEN",
        "THM-2051 now closes the relation-dissociated",
        "the final `nc2 : NC2` theorem remain",
    ),
    "00-navigation/START-HERE.md": (
        "14 total runners",
        "is **OPEN**",
        "uniform good period `q <= 25` is false",
        "uniform emptiness of the twelve-speed sporadic tight branch remains open",
        "unrestricted GMC(2)",
        "**PROVED in repo",
        "two-pair Poisson conjecture false",
        "THM-2047](../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md)",
        "Every counterexample lies on a bounded rational plane",
    ),
    "00-navigation/CURRENT-FRONTIER.md": (
        "## LRC(14)",
        "**OPEN.**",
        "first good period `q=27`",
        "disc_v>=6|G'_{~v}|^2",
        "proves NC2 and hence unrestricted GMC(2)",
        "GMC is false for every dimension at least 3",
        "THM-2044",
        "proves the corresponding labelled phase-height carrier",
        "PROVED finite-circuit alternative",
        "THM-2052",
        "gmc2_of_nc2` is conditional",
    ),
    "01-canon/ACTIVE-GUARDRAILS.md": (
        "No uniform `q <= 25` good-period theorem",
        "Uniform twelve-speed sporadic emptiness is OPEN",
        "HYP-8815 is a heuristic, not a disproof characterization",
        "A shared Pascal array is not a geometric bridge",
        "Braid localization does not factor every wall object",
        "Poisson rank two is not DC(2) or planar JC",
        "A thickened safe set is not an ordinary toric complement",
        "Antisymmetry is not the whole tournament-game or torus theorem",
        "Diagonal additive energy is not the LRC relation lattice",
        "NC2/GMC(2) is proved, not fully formalized",
    ),
    "00-navigation/LRC14-PROOF-MAP.md": (
        "## 2026-07-21 current control panel",
        "rank-11 bounded support-at-most-3 code W",
        "### Mandatory hostile controls",
    ),
    "00-navigation/LRC-TECHNIQUE-INDEX.md": (
        "SEARCHABLE ATLAS, NOT STARTUP TRUTH",
        "## 2026-07-21 current-use overlay",
        "## LTI-532 - Rank-eleven relation-code dispatcher",
    ),
}

# Exact retired summaries that are especially likely to send a new session
# down a closed or misattributed route. Historical ledgers may preserve them;
# the bounded maintained surface may not.
FORBIDDEN_STARTUP_TEXT = {
    "README.md": (
        "Every counterexample is small-relation structured",
    ),
    "00-navigation/START-HERE.md": (
        "every counterexample has a 2--5-term integer relation of height at most `2^20`",
    ),
    "01-canon/theorems/THM-2010-new-tournament-invariant-sequences.md": (
        "Breen--Stover--Yates",
        "Breen–Stover–Yates",
    ),
}

RUNTIME_STATE_PATHS = (
    ".machine-id",
    "agents/.read-log.json",
    "agents/.session-state.json",
)

DOCUMENTED_PYTHON_TOOLS = (
    "agents/start_session.py",
    "agents/processor.py",
    "agents/checkpoint_session.py",
    "agents/finish_session.py",
)

DOCUMENTED_SHELL_TOOLS = (
    "agents/new_session_worktree.sh",
    "agents/cleanup_session_worktrees.sh",
)

LINK_RE = re.compile(r"!?(?:\[[^\]]*\])\(([^)\n]+)\)")


def git(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ("git", *args), cwd=REPO, text=True, capture_output=True, check=False
    )


def local_link_parts(raw: str) -> tuple[str, str] | None:
    target = raw.strip()
    if target.startswith("<") and target.endswith(">"):
        target = target[1:-1]
    if target.startswith(("http://", "https://", "mailto:")):
        return None
    # Markdown titles are not used on the maintained surface. Flag rather than
    # guessing if one ever appears and breaks this simple parser.
    path, separator, fragment = target.partition("#")
    return unquote(path), unquote(fragment) if separator else ""


def markdown_anchors(path: Path) -> set[str]:
    anchors: set[str] = set()
    counts: Counter[str] = Counter()
    text = path.read_text(encoding="utf-8", errors="replace")
    for match in re.finditer(r"^#{1,6}\s+(.+?)\s*#*\s*$", text, re.MULTILINE):
        heading = match.group(1)
        heading = re.sub(r"\[([^]]+)\]\([^)]+\)", r"\1", heading)
        heading = re.sub(r"<[^>]+>", "", heading)
        heading = heading.replace("`", "").replace("*", "").replace("_", "_")
        slug = re.sub(r"[^\w\- ]", "", heading.casefold())
        slug = re.sub(r"\s", "-", slug).strip("-")
        if not slug:
            continue
        suffix = counts[slug]
        counts[slug] += 1
        anchors.add(slug if suffix == 0 else f"{slug}-{suffix}")
    return anchors


def read_required(relative: str, errors: list[str]) -> str:
    path = REPO / relative
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        errors.append(f"{relative}: required truth source is unreadable: {exc}")
        return ""


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
        byte_count = len(text.encode("utf-8"))
        if byte_count > STARTUP_BYTE_BUDGET:
            errors.append(
                f"{relative}: {byte_count} bytes exceeds startup budget "
                f"{STARTUP_BYTE_BUDGET}"
            )
        longest_line = max((len(line.encode("utf-8")) for line in lines), default=0)
        if longest_line > MAX_STARTUP_LINE_BYTES:
            errors.append(
                f"{relative}: longest line is {longest_line} bytes; startup maximum "
                f"is {MAX_STARTUP_LINE_BYTES}"
            )

        for raw in LINK_RE.findall(text):
            parts = local_link_parts(raw)
            if parts is None:
                continue
            target, fragment = parts
            resolved = path if not target else (
                Path(target) if Path(target).is_absolute() else path.parent / target
            )
            if not resolved.exists():
                errors.append(f"{relative}: broken local link {raw!r}")
                continue
            if (
                fragment
                and resolved.is_file()
                and resolved.suffix.casefold() == ".md"
                and not re.fullmatch(r"L\d+", fragment, flags=re.IGNORECASE)
                and fragment.casefold() not in markdown_anchors(resolved)
            ):
                errors.append(f"{relative}: broken Markdown fragment {raw!r}")

    for relative, sentinels in HEADLINE_SENTINELS.items():
        path = REPO / relative
        if not path.is_file():
            errors.append(f"{relative}: headline truth source is missing")
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for sentinel in sentinels:
            if sentinel not in text:
                errors.append(f"{relative}: missing headline truth sentinel {sentinel!r}")

    for relative, forbidden in FORBIDDEN_STARTUP_TEXT.items():
        text = read_required(relative, errors)
        for phrase in forbidden:
            if phrase in text:
                errors.append(f"{relative}: retired startup claim has returned: {phrase!r}")

    start = read_required("00-navigation/START-HERE.md", errors)
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

    for relative in RUNTIME_STATE_PATHS:
        tracked_state = git("ls-files", "--", relative).stdout.strip()
        if tracked_state:
            errors.append(f"{relative} is tracked; runtime state must stay local")

    for relative in (*DOCUMENTED_PYTHON_TOOLS, *DOCUMENTED_SHELL_TOOLS, "agents/check_docs.py"):
        if not os.access(REPO / relative, os.X_OK):
            errors.append(f"{relative}: expected executable bit is missing")

    for relative in DOCUMENTED_PYTHON_TOOLS:
        smoke = subprocess.run(
            (sys.executable, str(REPO / relative), "--help"),
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        if smoke.returncode != 0:
            errors.append(f"{relative}: documented --help smoke test failed")
    for relative in DOCUMENTED_SHELL_TOOLS:
        smoke = subprocess.run(
            ("bash", "-n", str(REPO / relative)),
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        if smoke.returncode != 0:
            errors.append(f"{relative}: shell syntax check failed")

    mistakes = read_required("01-canon/MISTAKES.md", errors)
    recent_ids = [
        value
        for value in re.findall(r"^## MISTAKE-(\d+)\b", mistakes, re.MULTILINE)
        if int(value) >= 209
    ]
    for value, count in sorted(Counter(recent_ids).items(), key=lambda item: int(item[0])):
        if count > 1:
            errors.append(f"MISTAKES.md: current-range MISTAKE-{value} occurs {count} times")

    hypotheses = read_required("05-knowledge/hypotheses/INDEX.md", errors)
    digest_marker = "# Hypothesis Log — Index"
    if not hypotheses.startswith("> **CURRENT DIGEST"):
        errors.append("hypotheses/INDEX.md: current digest must remain first")
    if digest_marker not in hypotheses:
        errors.append("hypotheses/INDEX.md: historical boundary is missing")
    current_digest = hypotheses.split(digest_marker, 1)[0]
    current_hypothesis_ids = [
        value
        for value in re.findall(r"^- \*\*HYP-(\d+)\b", current_digest, re.MULTILINE)
        if int(value) >= 8800
    ]
    for value, count in sorted(
        Counter(current_hypothesis_ids).items(), key=lambda item: int(item[0])
    ):
        if count > 1:
            errors.append(
                f"hypotheses/INDEX.md: current-range HYP-{value} occurs {count} times"
            )

    session_log = read_required("00-navigation/SESSION-LOG.md", errors)
    if not session_log.startswith("> **CURRENT-TRUTH WARNING"):
        errors.append("SESSION-LOG.md: current-truth warning must remain first")

    empty_topic = subprocess.run(
        (sys.executable, str(REPO / "agents/start_session.py"), "--topic", ""),
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    if empty_topic.returncode == 0:
        errors.append("start_session.py: an empty topic must be rejected")
    generic_topic = subprocess.run(
        (
            sys.executable,
            str(REPO / "agents/start_session.py"),
            "--topic",
            "math proof research",
        ),
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    if generic_topic.returncode == 0:
        errors.append("start_session.py: an all-generic topic must be rejected")

    if errors:
        print("Agent-facing documentation check FAILED:")
        for error in errors:
            print(f"- {error}")
        return 1

    print(
        f"Agent-facing documentation check passed: {len(STARTUP_DOCS)} maintained "
        "documents, bounded size, local links, truth sentinels, configs, runtime "
        "policy, tool/startup smoke tests, and current IDs."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
