#!/usr/bin/env python3
"""Validate the bounded agent-facing truth surface.

The checker deliberately derives theorem status from frontmatter.  It keeps a
single small table only for recent frontier promotions whose accidental
rollback would misroute a new session; prose checks are mechanism-level, not
word-for-word snapshots of one editorial rendering.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

from doc_surface import (
    LINE_BUDGETS,
    MAX_STARTUP_LINE_BYTES,
    PREFIX_LINE_BUDGETS,
    SEARCHABLE_PREFIX_END,
    STARTUP_BYTE_BUDGET,
    STARTUP_DOCS,
)
from start_session import UNPROVED_CANDIDATE_STATUSES, file_status


ROOT = Path(__file__).resolve().parent.parent
errors: list[str] = []


def fail(message: str) -> None:
    errors.append(message)


def read(relative: str) -> str:
    path = ROOT / relative
    try:
        return path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        fail(f"{relative}: unreadable: {exc}")
        return ""


def require(relative: str, *needles: str) -> None:
    body = read(relative)
    for needle in needles:
        if needle not in body:
            fail(f"{relative}: missing current-truth token {needle!r}")


def forbid(relative: str, *needles: str) -> None:
    body = read(relative)
    for needle in needles:
        if needle in body:
            fail(f"{relative}: stale/forbidden token {needle!r}")


def maintained_lines(relative: str) -> list[str]:
    body = read(relative)
    marker = SEARCHABLE_PREFIX_END.get(relative)
    if marker:
        if marker not in body:
            fail(f"{relative}: missing historical boundary {marker!r}")
        else:
            body = body.split(marker, 1)[0]
    return body.splitlines()


def theorem_path(number: int) -> str | None:
    hits = sorted((ROOT / "01-canon/theorems").glob(f"THM-{number}-*.md"))
    if len(hits) != 1:
        fail(f"THM-{number}: expected one theorem file, found {len(hits)}")
        return None
    return hits[0].relative_to(ROOT).as_posix()


def run(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=ROOT,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )


# Bounded startup surface.
for relative in STARTUP_DOCS:
    if not (ROOT / relative).is_file():
        fail(f"{relative}: missing startup document")

for relative, budget in LINE_BUDGETS.items():
    count = len(maintained_lines(relative))
    if count > budget:
        fail(f"{relative}: {count} maintained lines exceeds budget {budget}")

for relative, budget in PREFIX_LINE_BUDGETS.items():
    count = len(maintained_lines(relative))
    if count > budget:
        fail(f"{relative}: {count} maintained-prefix lines exceeds budget {budget}")

startup_bytes = 0
for relative in STARTUP_DOCS:
    raw = (ROOT / relative).read_bytes() if (ROOT / relative).is_file() else b""
    # The repository's evidence convention is LF-normalized bytes.  Count the
    # same canonical surface on LF and CRLF checkouts instead of charging one
    # extra byte per line on Windows.
    normalized = raw.replace(b"\r\n", b"\n")
    startup_bytes += len(normalized)
    for number, line in enumerate(normalized.splitlines(), 1):
        if len(line) > MAX_STARTUP_LINE_BYTES:
            fail(f"{relative}:{number}: {len(line)} bytes exceeds line budget")
if startup_bytes > STARTUP_BYTE_BUDGET:
    fail(f"startup surface: {startup_bytes} bytes exceeds {STARTUP_BYTE_BUDGET}")

conflicts = run("git", "grep", "-n", "-E", r"^(<<<<<<<|>>>>>>>)")
if conflicts.returncode == 0:
    fail("tracked files contain merge-conflict markers")
elif conflicts.returncode not in (0, 1):
    fail(f"git grep conflict scan failed: {conflicts.stderr.strip()}")

# Byte-zero current-truth sentinels.
for relative, prefix in {
    "00-navigation/SESSION-LOG.md": "> **CURRENT-TRUTH WARNING",
    "05-knowledge/hypotheses/INDEX.md": "> **CURRENT DIGEST",
    "07-reflections/bypassing-the-gmc2-dvdk-dependency-for-the-unique-channel-class-boxeph-S230.md": "> **CORRECTION",
    "07-reflections/the-dvdk-residual-is-one-unramified-hensel-small-root-product-a-formalization-map-for-thm-2067-deathstar-S106.md": "> **STATUS CORRECTION",
    "07-reflections/eliminating-dvdk-for-the-residual-12-percent-via-monomial-certificates-boxeph-S231.md": "> **SCOPE CORRECTION",
}.items():
    if not read(relative).startswith(prefix):
        fail(f"{relative}: byte-zero sentinel displaced")

# Stable headline truth.  Synonyms may change; these load-bearing tokens may not.
require(
    "00-navigation/START-HERE.md",
    "LRC(14)", "OPEN", "q <= 25", "THM-2084", "THM-2081--2087",
    "THM-2088--2093", "GMC2Main.gmc2", "singlePolyCrux_holds",
    "Anchor / Niche / Wildcard", "MISTAKE-260", "MISTAKE-261", "28,393", "THM-2363", "THM-2368", "THM-2371", "THM-2376",
)
require(
    "00-navigation/CURRENT-FRONTIER.md",
    "LRC(14)", "THM-2081", "THM-2082", "THM-2083", "THM-2084",
    "THM-2085", "THM-2086", "THM-2085/2087", "THM-2088", "THM-2089",
    "THM-2090", "THM-2092", "MISTAKE-240", "DvdK1", "28,393", "THM-2363", "THM-2368", "THM-2371", "THM-2376",
)
require(
    "01-canon/ACTIVE-GUARDRAILS.md",
    "No uniform `q<=25` theorem", "THM-2095", "MISTAKE-240",
    "Support is not an indexed multiset", "CLAIMED",
    "MISTAKE-260", "MISTAKE-261", "28,393", "THM-2363",
)
require(
    "00-navigation/SESSION-LOG.md",
    "LRC(14)", "OPEN", "THM-2092", "MISTAKE-240", "MISTAKE-241",
    "DvdK1", "GMC2Main.gmc2", "singlePolyCrux_holds",
    "Anchor / Niche / Wildcard",
)
require(
    "05-knowledge/hypotheses/INDEX.md",
    "THM-2091", "THM-2094", "THM-2095", "THM-2096", "CLAIMED STUB", "THM-2092", "HYP-8931", "HYP-8932",
    "HYP-8935", "MISTAKE-240", "MISTAKE-241", "DvdK1",
    "GMC2Main.gmc2", "singlePolyCrux_holds",
)
require(
    "01-canon/MISTAKES.md",
    "MISTAKE-240", "lambda=0", "F=empty", "MISTAKE-241",
    "floating root asymptotics", "GMC2OrbitProduct.lean",
)
require(
    "05-knowledge/reference/CORE-PAPERS.md",
    "Vaaler", "Theorem 19", "THM-2085", "Does not prove",
)

forbid(
    "00-navigation/START-HERE.md",
    "THM-2084/2085/2086 are RESERVED",
)
forbid(
    "01-canon/ACTIVE-GUARDRAILS.md",
    "THM-2084, THM-2085, and THM-2086 are RESERVED",
)
forbid(
    "00-navigation/CURRENT-FRONTIER.md",
    "A two-target graph lands",
    "complete one-sparse no-landing locus",
)
forbid(
    "00-navigation/LRC14-PROOF-MAP.md",
    "A two-target graph lands",
    "complete one-sparse no-landing locus",
)
forbid(
    "05-knowledge/hypotheses/INDEX.md",
    "THM-2091 / THM-2092 (CLAIMED STUBS",
)

# One centralized recent-frontier status gate.
expected_recent = {
    2084: "PROVED", 2085: "PROVED", 2086: "PROVED", 2087: "PROVED",
    2088: "PROVED", 2089: "PROVED", 2090: "PROVED", 2091: "PROVED",
    2092: "PROVED", 2093: "PROVED", 2094: "PROVED", 2095: "PROVED",
    2096: "PROVED", 2097: "PROVED",
    2118: "PROVED", 2119: "PROVED", 2120: "PROVED", 2121: "PROVED",
    2122: "PROVED", 2123: "PROVED", 2124: "PROVED", 2125: "PROVED",
    2126: "PROVED", 2127: "PROVED", 2129: "PROVED", 2130: "PROVED",
    2350: "PROVED", 2354: "PROVED", 2356: "PROVED", 2360: "PROVED", 2362: "PROVED",
    2363: "PROVED", 2364: "PROVED",
}
for number, expected in expected_recent.items():
    relative = theorem_path(number)
    if relative is None:
        continue
    actual = file_status(relative)
    if actual != expected:
        fail(f"{relative}: expected status {expected}, found {actual}")

# No proved theorem may import an explicitly unproved candidate.
theorem_files = sorted((ROOT / "01-canon/theorems").glob("THM-*.md"))
status_by_id: dict[str, str] = {}
body_by_id: dict[str, str] = {}
path_by_id: dict[str, str] = {}
# The historical corpus contains documented legacy collisions.  The current
# allocation epoch begins above all of them; enforce uniqueness there so a
# raced reservation cannot silently overwrite the proof-graph dictionary.
UNIQUE_THEOREM_ID_FLOOR = 3300
for path in theorem_files:
    body = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"^id:\s*(THM-\d+)\s*$", body, re.MULTILINE)
    if match:
        relative = path.relative_to(ROOT).as_posix()
        theorem_id = match.group(1)
        if theorem_id in path_by_id:
            if int(theorem_id.removeprefix("THM-")) >= UNIQUE_THEOREM_ID_FLOOR:
                fail(
                    f"{theorem_id}: duplicate theorem ID in "
                    f"{path_by_id[theorem_id]} and {relative}"
                )
            continue
        path_by_id[theorem_id] = relative
        status_by_id[theorem_id] = file_status(relative)
        body_by_id[theorem_id] = body
for theorem_id, status in status_by_id.items():
    if status != "PROVED":
        continue
    header = body_by_id[theorem_id].split("---", 2)[1] if "---" in body_by_id[theorem_id] else ""
    dependency_block = re.search(
        r"^depends_on:\s*\n"
        r"(?P<rows>(?:[ \t]+-[^\r\n]*(?:\r?\n|$))*)",
        header,
        re.MULTILINE,
    )
    dependencies = re.findall(
        r"^\s*-\s*(THM-\d+)(?:-[^\s#]+)?(?:\s+#.*)?\s*$",
        dependency_block.group("rows") if dependency_block else "",
        re.MULTILINE,
    )
    for dependency in dependencies:
        if status_by_id.get(dependency) in UNPROVED_CANDIDATE_STATUSES:
            fail(f"{theorem_id}: proved dependency graph imports {dependency}")

# Router smoke test: an actual current unproved theorem must be visibly
# quarantined.  Select from frontmatter instead of pinning a theorem that can
# later be promoted (THM-741 was the original fixture and is now proved).
smoke_candidates = sorted(
    (
        (status != "CLAIMED", theorem_id, status)
        for theorem_id, status in status_by_id.items()
        if status in UNPROVED_CANDIDATE_STATUSES
    )
)
if not smoke_candidates:
    fail("start_session.py: no RESERVED/CLAIMED theorem available for router smoke")
else:
    _, smoke_id, smoke_status = smoke_candidates[0]
    smoke = run(
        sys.executable,
        "agents/start_session.py",
        "--topic",
        smoke_id,
        "--max-matches",
        "8",
    )
    if smoke.returncode != 0:
        fail(f"start_session.py smoke failed: {smoke.stderr.strip()}")
    else:
        for token in (
            "Unproved candidates (not results)",
            f"[{smoke_status}]",
            smoke_id,
        ):
            if token not in smoke.stdout:
                fail(f"start_session.py: unproved-route smoke lacks {token!r}")

posture = run(sys.executable, "agents/start_session.py", "--topic", "LRC(14)", "--max-matches", "4")
if posture.returncode != 0:
    fail(f"start_session.py LRC smoke failed: {posture.stderr.strip()}")
else:
    for token in ("Inheritance:", "Niche", "Wildcard", "META-PATTERNS"):
        if token not in posture.stdout:
            fail(f"start_session.py: method posture lacks {token!r}")

if errors:
    print("Agent-facing documentation check FAILED:")
    for error in errors:
        print(f"- {error}")
    raise SystemExit(1)

print(
    "Agent-facing documentation check PASSED: bounded startup surface, "
    "frontmatter-derived frontier status, correction sentinels, and router smoke."
)
