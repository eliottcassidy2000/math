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
    CANONICAL_ROUTES,
    FRONTIER_EPOCH,
    LINE_BUDGETS,
    MAX_STARTUP_PACKET_BYTES,
    MAX_STARTUP_LINE_BYTES,
    PREFIX_LINE_BUDGETS,
    SEARCHABLE_PREFIX_END,
    STARTUP_BYTE_BUDGET,
    STARTUP_DOCS,
)


REPO = Path(__file__).resolve().parent.parent

HEADLINE_SENTINELS = {
    "README.md": (
        "LRC(14) is OPEN",
        "THM-2051 now closes the relation-dissociated",
        "exact transverse deck `D_N(m)`",
        "THM-2055/2056",
        "THM-2057 closes two scaled AP one-tail planes",
        "THM-2059 gives an exact",
        "`GMC2NC2Capstone`",
        "unique minimum-mass balanced channel",
        "MISTAKE-230--233",
        "HYP-8885",
    ),
    "00-navigation/START-HERE.md": (
        "14 total runners",
        "is **OPEN**",
        "uniform good period `q <= 25` is false",
        "uniform emptiness of the twelve-speed sporadic tight branch remains open",
        "unrestricted GMC(2)",
        "**PROVED in repo",
        "two-pair Poisson conjecture false",
        "A unique minimal balanced channel needs no DvdK seed",
        "THM-2047](../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md)",
        "Every counterexample is in a finite labelled code/deck/fan atlas",
        "THM-2057 closes two scaled AP one-tail families",
        "THM-2059 exactly joins arbitrary-clock",
        "THM-2058 is an unproved empty stub",
        "MISTAKE-230--233",
        "HYP-8885",
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
        "transverse deck D_N(m)",
        "Kelvin-Farey certificate",
        "THM-2057",
        "THM-2059",
        "RESERVED / UNPROVED EMPTY STUB",
        "TournamentH7.GMC2Formalization",
        "TournamentH7.GMC2NC2Capstone",
        "unique minimal balanced channel",
        "HYP-8890",
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
        "HYP-8878 removes that citation only",
        "THM-2059 only joins packets",
        "Equal ranks do not identify lattices or tournaments",
        "Paley spectra do not assign LRC roles to small primes",
        "THM-2053 has no Heegner discriminant `-7`",
        "THM-2058 is a RESERVED / UNPROVED EMPTY STUB",
        "Fiber cardinality is observable-relative Hartley ambiguity",
        "Scaled clocks are not automatically modular cusps or Frobenius",
        "`f14`, genus, and discriminant `-7` are not the LRC obstruction",
        "HYP-8885 may use “cusp” only as a",
    ),
    "00-navigation/LRC14-PROOF-MAP.md": (
        "## 2026-07-21 current control panel",
        "bounded support-at-most-3 code of rank >=11",
        "transverse deck D_N(m)>=1/14",
        "finite Farey address; listed rays unresolved",
        "scaled AP one-tail clock/binding leaves",
        "### Mandatory hostile controls",
    ),
    "00-navigation/LRC-TECHNIQUE-INDEX.md": (
        "SEARCHABLE ATLAS, NOT STARTUP TRUTH",
        "## 2026-07-21 current-use overlay",
        "## LTI-532 - Rank-eleven relation-code dispatcher",
        "Transverse deck",
        "Missing-clock binding",
        "CRT packet join",
    ),
}

POLICY_SENTINELS = {
    "AGENTS.md": (
        "Anchor / Niche / Wildcard",
        "determine **why** a claim is true or false",
        "keep a board of 3–7 live concepts",
        "source, target, map, preserved predicate",
        "give every small mathematical compulsion a cheap hostile probe",
        "00-navigation/META-PATTERNS.md",
        "RESERVED / UNPROVED EMPTY STUB",
        "Never prepend session prose above a maintained warning",
    ),
    "00-navigation/RESEARCH-PROTOCOL.md": (
        "The session portfolio: Anchor / Niche / Wildcard",
        "Keep an active concept board",
        "Generate perspectives procedurally",
        "The connection contract",
        "Explain why, not only whether",
        "Make the process cumulative",
    ),
}

# Exact retired summaries that are especially likely to send a new session
# down a closed or misattributed route. Historical ledgers may preserve them;
# the bounded maintained surface may not.
FORBIDDEN_STARTUP_TEXT = {
    "README.md": (
        "Every counterexample is small-relation structured",
        "height at most `2^21`",
        "pointed rank-or-Euler transport is the live prize",
        "root wiring, and the final `nc2 : NC2` theorem remain",
    ),
    "00-navigation/START-HERE.md": (
        "every counterexample has a 2--5-term integer relation of height at most `2^20`",
        "height at most `2^21`",
        "The remaining jump is pointed Euler survival",
        "They are not root-imported",
    ),
    "00-navigation/CURRENT-FRONTIER.md": (
        "but are not root-imported",
        "6. **Pointed plane transport.**",
        "The incoming THM-2054 relative-Fejer program",
        "THM-2058 supplies the primitive phase-packet",
    ),
    "01-canon/ACTIVE-GUARDRAILS.md": (
        "height-`2^21` relation",
        "MISTAKE-225 corrects HYP-8865",
    ),
    "00-navigation/LRC14-PROOF-MAP.md": (
        "height <=2^21",
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


def maintained_text(relative: str, errors: list[str]) -> str:
    text = read_required(relative, errors)
    marker = SEARCHABLE_PREFIX_END.get(relative)
    if not marker:
        return text
    if marker not in text:
        errors.append(f"{relative}: maintained-prefix boundary {marker!r} is missing")
        return text
    return text.split(marker, 1)[0]


def validate_surface(
    relative: str,
    text: str,
    errors: list[str],
    *,
    line_budget: int | None,
) -> None:
    path = REPO / relative
    lines = text.splitlines()
    if line_budget is not None and len(lines) > line_budget:
        errors.append(
            f"{relative}: {len(lines)} maintained lines exceeds budget {line_budget}"
        )
    byte_count = len(text.encode("utf-8"))
    if byte_count > STARTUP_BYTE_BUDGET:
        errors.append(
            f"{relative}: {byte_count} maintained bytes exceeds budget "
            f"{STARTUP_BYTE_BUDGET}"
        )
    longest_line = max((len(line.encode("utf-8")) for line in lines), default=0)
    if longest_line > MAX_STARTUP_LINE_BYTES:
        errors.append(
            f"{relative}: longest maintained line is {longest_line} bytes; maximum "
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


def main() -> int:
    errors: list[str] = []

    for relative in STARTUP_DOCS:
        path = REPO / relative
        if not path.is_file():
            errors.append(f"missing maintained startup document: {relative}")
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        validate_surface(
            relative, text, errors, line_budget=LINE_BUDGETS.get(relative)
        )

    for relative in dict.fromkeys((*CANONICAL_ROUTES, *SEARCHABLE_PREFIX_END)):
        if relative in STARTUP_DOCS:
            continue
        if not (REPO / relative).is_file():
            errors.append(f"missing maintained route: {relative}")
            continue
        validate_surface(
            relative,
            maintained_text(relative, errors),
            errors,
            line_budget=PREFIX_LINE_BUDGETS.get(relative),
        )

    for relative, sentinels in HEADLINE_SENTINELS.items():
        path = REPO / relative
        if not path.is_file():
            errors.append(f"{relative}: headline truth source is missing")
            continue
        text = maintained_text(relative, errors)
        for sentinel in sentinels:
            if sentinel not in text:
                errors.append(f"{relative}: missing headline truth sentinel {sentinel!r}")

    for relative, forbidden in FORBIDDEN_STARTUP_TEXT.items():
        text = maintained_text(relative, errors)
        normalized = " ".join(text.split())
        for phrase in forbidden:
            if phrase in text or " ".join(phrase.split()) in normalized:
                errors.append(f"{relative}: retired startup claim has returned: {phrase!r}")

    for relative, sentinels in POLICY_SENTINELS.items():
        text = read_required(relative, errors)
        for sentinel in sentinels:
            if sentinel not in text:
                errors.append(f"{relative}: missing research-policy sentinel {sentinel!r}")

    for relative in (
        "README.md",
        "00-navigation/START-HERE.md",
        "00-navigation/CURRENT-FRONTIER.md",
        "01-canon/ACTIVE-GUARDRAILS.md",
    ):
        if FRONTIER_EPOCH not in read_required(relative, errors):
            errors.append(f"{relative}: frontier epoch is not {FRONTIER_EPOCH}")

    meta = read_required("00-navigation/META-PATTERNS.md", errors)
    meta_cards = list(
        re.finditer(r"^## ([^\n]+)\n(.*?)(?=^## |\Z)", meta, flags=re.M | re.S)
    )
    required_card_fields = (
        "**Trigger:**", "**Action:**", "**Mechanism:**",
        "**Counterindication:**", "**Evidence:**",
    )
    for card in meta_cards:
        missing = [field for field in required_card_fields if field not in card.group(2)]
        if missing:
            errors.append(
                f"META-PATTERNS.md: card {card.group(1)!r} lacks {', '.join(missing)}"
            )

    papers = read_required("05-knowledge/reference/CORE-PAPERS.md", errors)
    paper_cards = list(
        re.finditer(r"^### ([^\n]+)\n(.*?)(?=^### |^## |\Z)", papers, flags=re.M | re.S)
    )
    for card in paper_cards:
        body = card.group(2)
        if "**Primary / freshness:**" not in body:
            continue
        required = (
            "**Imported role:**",
            "**Does not prove:**",
        )
        for field in required:
            if field not in body:
                errors.append(f"CORE-PAPERS.md: card {card.group(1)!r} lacks {field}")
        if not re.search(r"\*\*Repo (?:consumer|consumers|landing point|landing points):\*\*", body):
            errors.append(f"CORE-PAPERS.md: card {card.group(1)!r} lacks a repo consumer")

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

    for value in (230, 231, 232, 233, 234):
        if f"## MISTAKE-{value}" not in mistakes:
            errors.append(f"MISTAKES.md: current correction MISTAKE-{value} is missing")

    reserved = read_required(
        "01-canon/theorems/THM-2058-primitive-phase-packets-and-deck-fan-intervals.md",
        errors,
    )
    for sentinel in (
        "RESERVED / UNPROVED EMPTY STUB",
        "No theorem statement, proof, script, or result is present",
        "depends_on: []",
    ):
        if sentinel not in reserved:
            errors.append(f"THM-2058 reservation lacks {sentinel!r}")
    if re.search(r"^status:\s*>?\s*PROVED\b", reserved, re.MULTILINE):
        errors.append("THM-2058 reservation is mislabeled PROVED")

    crt_packet = read_required(
        "01-canon/theorems/THM-2059-crt-fiber-product-phase-packet.md", errors
    )
    for sentinel in (
        "status: >\n  PROVED.",
        "P_N(C;a,w)=sum_(j mod d) alpha_j beta_j",
        "zero overlap rejects only that clock grid",
    ):
        if sentinel not in crt_packet:
            errors.append(f"THM-2059 carrier lacks {sentinel!r}")

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
    for value in (8878, 8885, 8890, 8895, 8900):
        if f"HYP-{value}" not in current_digest:
            errors.append(
                f"hypotheses/INDEX.md: current digest lacks incoming HYP-{value} routing"
            )
    if "THM-2059 / CRT PHASE-PACKET CARRIER" not in current_digest:
        errors.append("hypotheses/INDEX.md: current digest lacks THM-2059 routing")

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

    control_topic = subprocess.run(
        (
            sys.executable,
            str(REPO / "agents/start_session.py"),
            "--topic",
            "LRC14\nforged-heading",
        ),
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    if control_topic.returncode == 0:
        errors.append("start_session.py: a topic with control characters must be rejected")

    routing_smoke = subprocess.run(
        (
            sys.executable,
            str(REPO / "agents/start_session.py"),
            "--topic",
            "Jacobian conjecture dimension three",
            "--recent",
            "1",
            "--max-matches",
            "8",
        ),
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    if routing_smoke.returncode != 0:
        errors.append("start_session.py: bounded routing smoke test failed")
    packet = routing_smoke.stdout
    if len(packet.encode("utf-8")) > MAX_STARTUP_PACKET_BYTES:
        errors.append(
            f"start_session.py: packet exceeds {MAX_STARTUP_PACKET_BYTES} bytes"
        )
    longest_packet_line = max(
        (len(line.encode("utf-8")) for line in packet.splitlines()), default=0
    )
    if longest_packet_line > MAX_STARTUP_LINE_BYTES:
        errors.append(
            f"start_session.py: emitted line is {longest_packet_line} bytes; maximum "
            f"is {MAX_STARTUP_LINE_BYTES}"
        )
    if "JC_n) — DISPROVED" in packet or "PROBLEM-LEDGER.md:51:" in packet:
        errors.append("start_session.py: legacy problem-ledger suffix leaked as current")
    for sentinel in ("== Session posture ==", "truth precedence first"):
        if sentinel not in packet:
            errors.append(f"start_session.py: routing smoke lacks {sentinel!r}")

    correction_smoke = subprocess.run(
        (
            sys.executable,
            str(REPO / "agents/start_session.py"),
            "--topic",
            "Heegner determinant gate",
            "--recent",
            "1",
            "--max-matches",
            "4",
        ),
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
    )
    if "MISTAKE-229" not in correction_smoke.stdout:
        errors.append("start_session.py: topic-matched MISTAKE-229 was not surfaced")

    entropy_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "arithmetic entropy tournament score fiber",
            "--recent", "1", "--max-matches", "6",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    if "MISTAKE-231" not in entropy_smoke.stdout:
        errors.append("start_session.py: entropy topic did not surface MISTAKE-231")

    modular_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "f14 modular cusp LRC obstruction",
            "--recent", "1", "--max-matches", "6",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    if "MISTAKE-233" not in modular_smoke.stdout:
        errors.append("start_session.py: modular-cusp topic did not surface MISTAKE-233")

    semigroup_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "coprime interval return semigroup DvdK",
            "--recent", "1", "--max-matches", "6",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    if "MISTAKE-234" not in semigroup_smoke.stdout:
        errors.append("start_session.py: return-semigroup topic did not surface MISTAKE-234")

    reservation_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "THM-2058 primitive phase packets",
            "--recent", "1", "--max-matches", "12",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    reservation_packet = reservation_smoke.stdout
    if "Reservations (not results):" not in reservation_packet or "[RESERVED]" not in reservation_packet:
        errors.append("start_session.py: THM-2058 is not routed as a reservation")
    canon_block = reservation_packet.partition("  Canon:")[2].partition(
        "  Reservations (not results):"
    )[0]
    if "THM-2058" in canon_block:
        errors.append("start_session.py: THM-2058 leaked into the proved canon group")
    for sentinel in (
        "[OPEN] 05-knowledge/hypotheses/HYP-8846-lrc14-pointed-plane-transport.md",
        "[OPEN] 05-knowledge/hypotheses/HYP-8871-lrc14-owner-sector-klein-sail-automaton.md",
        "[LEDGER] 05-knowledge/hypotheses/INDEX.md",
    ):
        if sentinel not in reservation_packet:
            errors.append(
                "start_session.py: folded-status routing lacks " + repr(sentinel)
            )

    crt_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "THM-2059 CRT phase packet",
            "--recent", "1", "--max-matches", "8",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    if (
        "[PROVED] 01-canon/theorems/THM-2059-crt-fiber-product-phase-packet.md"
        not in crt_smoke.stdout
    ):
        errors.append("start_session.py: THM-2059 is not routed as proved canon")

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
