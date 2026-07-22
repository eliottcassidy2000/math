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
from start_session import file_status


REPO = Path(__file__).resolve().parent.parent

HEADLINE_SENTINELS = {
    "README.md": (
        "LRC(14) is OPEN",
        "THM-2051 now closes the relation-dissociated",
        "exact transverse deck `D_N(m)`",
        "THM-2055/2056",
        "THM-2057 closes two scaled AP one-tail planes",
        "THM-2058/2059 give phase packets and",
        "THM-2060/2064 isolate a dyadic seam",
        "THM-2065 makes circuit-free templates finite",
        "THM-2066",
        "THM-2068 compresses that bounded",
        "THM-2069",
        "THM-2072",
        "THM-2073/2076",
        "depth-at-most-five",
        "THM-2075",
        "THM-2078",
        "THM-2079",
        "hereditarily primitive terminal",
        "`GMC2NC2Capstone`",
        "THM-2067 supplies its bare one-variable constant-term",
        "`DvdK1` and `HeightWitnessSupplier` remain explicit interfaces",
        "THM-2070 explains why a cofinite",
        "JC is false from dimension three; JC(2) and DC(2) remain open",
        "THM-2071 closes every quadratic-fiber",
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
        "THM-2067",
        "`DvdK1` and `HeightWitnessSupplier` remain explicit",
        "Support return is not coefficient noncancellation",
        "JC is false from dimension three; JC(2) and DC(2) remain open",
        "THM-2071 closes every quadratic-fiber",
        "two-pair Poisson conjecture false",
        "A unique minimal balanced channel needs no DvdK seed",
        "THM-2047](../01-canon/theorems/THM-2047-phase-height-toric-arrangement-for-lrc.md)",
        "Every counterexample is in a finite labelled code/deck/fan atlas",
        "THM-2057 closes two AP one-tail families",
        "THM-2065 makes circuit-free templates finite",
        "THM-2066 closes divisor-complete cores through `24`",
        "seven-clock minimum only inside `15..34`",
        "THM-2073/2076",
        "depth-at-most-five safe-child tower",
        "THM-2075",
        "THM-2078",
        "THM-2079",
        "unbounded hereditarily primitive terminal remains open",
        "THM-2069 proves the code/cogirth wheel",
        "keep THM-2074 reserved",
        "MISTAKE-230--233",
        "HYP-8885",
        "MISTAKE-235",
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
        "THM-2058's coprime intervals",
        "THM-2059",
        "THM-2060",
        "THM-2061",
        "THM-2062",
        "THM-2064",
        "THM-2065",
        "THM-2066's complementary owner-word bank",
        "THM-2068 proves that `{25,26,27,28,32,33,34}` is a minimum sub-bank inside",
        "THM-2072 proves",
        "THM-2073/2076 give a depth-at-most-five safe-child",
        "THM-2075 preserves nonempty-core component/Euler/owner data",
        "THM-2078",
        "THM-2079 pairs addresses",
        "hereditarily primitive quotient core",
        "THM-2069 proves the code/cogirth wheel",
        "THM-2074 (density one) remains reserved",
        "TournamentH7.GMC2Formalization",
        "TournamentH7.GMC2NC2Capstone",
        "`DvdK1` plus `HeightWitnessSupplier`",
        "formalize THM-2067 to discharge `DvdK1`",
        "THM-2070",
        "THM-2071",
        "no hypothetical JC(2) counterexample has a quadratic member",
        "HYP-8890",
    ),
    "01-canon/ACTIVE-GUARDRAILS.md": (
        "No uniform `q<=25` theorem",
        "Uniform twelve-speed sporadic emptiness is OPEN",
        "HYP-8815 is not a disproof characterization",
        "A shared Pascal array is not a geometric bridge",
        "Braid localization does not factor every wall object",
        "Poisson rank two is not DC(2) or JC(2)",
        "A thickened safe set is not an ordinary toric complement",
        "Antisymmetry is not a whole game/torus theorem",
        "Shared kernels are not bridges",
        "MISTAKE-226/235",
        "NC2/GMC(2) is proved, not fully formalized",
        "THM-2057 closes two planes; THM-2059 joins packets",
        "Equal ranks do not identify lattices",
        "Small-prime Paley spectra assign no LRC roles",
        "THM-2058 is a carrier, not LRC(14)",
        "THM-2060/2064 give capacity, not signed ownership",
        "THM-2061 is a proved reduction, not an empty-seam theorem",
        "THM-2062 filters heredity, not phase height",
        "THM-2065 is PROVED only as a reduction",
        "THM-2066 is a bounded closure",
        "THM-2068 is menu-relative",
        "No fixed finite bank is uniform",
        "THM-2073 ends at a terminal, not at LRC(14)",
        "depth at most five",
        "Safe-child homeomorphisms transport core safe sets",
        "Do not infer terminal `chi=0`",
        "THM-2069 is a code/cogirth theorem, not code existence or LRC(14)",
        "THM-2074 is RESERVED",
        "Paper proof is not full Lean proof",
        "still assumes `DvdK1` and `HeightWitnessSupplier`",
        "THM-2067 is bare existence, not effective DvdK",
        "Support feasibility is not weighted noncancellation",
        "Fiber size is observable-relative Hartley ambiguity",
        "Clocks are not modular cusps",
        "The binary symmetric Hessian is not full JC(2)",
        "THM-1330 is a necessary atlas, not a classification",
        "THM-2071 is a pencil gate, not a cover theorem",
    ),
    "00-navigation/LRC14-PROOF-MAP.md": (
        "## 2026-07-21 current control panel",
        "bounded support-at-most-3 code of rank >=11",
        "transverse deck D_N(m)>=1/14",
        "finite Farey address; listed rays unresolved",
        "scaled AP one-tail clock/binding leaves",
        "tail-sheet capacity -> dyadic seam",
        "hereditary-primitivity CRT wheel",
        "circuit-free relation rays become finite",
        "owner words close divisor-complete max<=24",
        "fixed-bank no-go -> adaptive safe children",
        "hereditarily primitive terminal",
        "depth at most five",
        "terminal maximum `<=24`",
        "mirror-complement",
        "THM-2069 proves the code/cogirth wheel",
        "THM-2074 remains reserved",
        "### Mandatory hostile controls",
    ),
    "00-navigation/LRC-TECHNIQUE-INDEX.md": (
        "SEARCHABLE ATLAS, NOT STARTUP TRUTH",
        "## 2026-07-21 current-use overlay",
        "## LTI-532 - Rank-eleven relation-code dispatcher",
        "Transverse deck",
        "Missing-clock binding",
        "CRT packet join",
        "circuit rays + hereditary wheel",
        "owner words / safe-child tower",
        "THM-2068 is menu-relative",
        "THM-2072 proves no fixed finite bank is uniform",
        "THM-2069 proves the code/cogirth wheel",
        "THM-2074 remains reserved",
        "depth at most five",
        "terminal maximum `<=24`",
        "mirror-complement",
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
        "THM-2058 is an empty reservation",
        "THM-2060/2061 are reserved proof",
        "THM-2065 is RESERVED",
        "one `sorry` remains",
        "exactly one sorry remains",
        "adaptive depth-at-most-eight safe-child tower",
        "THM-2071 forces a square leading coefficient",
        "first quadratic-fiber survivor",
        "THM-2069 and THM-2074 are RESERVED",
    ),
    "00-navigation/START-HERE.md": (
        "every counterexample has a 2--5-term integer relation of height at most `2^20`",
        "height at most `2^21`",
        "The remaining jump is pointed Euler survival",
        "They are not root-imported",
        "THM-2058 is an unproved empty stub",
        "THM-2060/2061 are reserved proof",
        "THM-2065 is an unproved",
        "THM-2065 is reserved",
        "exactly one sorry remains",
        "depth-at-most-eight",
        "first quadratic-fiber survivor",
        "keep THM-2069/2074 reserved",
    ),
    "00-navigation/CURRENT-FRONTIER.md": (
        "but are not root-imported",
        "6. **Pointed plane transport.**",
        "The incoming THM-2054 relative-Fejer program",
        "THM-2058 supplies the primitive phase-packet",
        "THM-2058 is an empty reservation",
        "Reserved proof candidates—not theorem inputs",
        "THM-2065's proposed circuit-ray collapse",
        "THM-2065 is RESERVED",
        "The Lean formalization has one sorry",
        "THM-2071 adds the next pencil gate",
        "after at most eight levels, ending",
        "THM-2069 (code/cogirth) and THM-2074 (density one) remain reserved",
    ),
    "01-canon/ACTIVE-GUARDRAILS.md": (
        "height-`2^21` relation",
        "MISTAKE-225 corrects HYP-8865",
        "are independent reserved proof candidates",
        "THM-2065 is RESERVED",
        "The formalization has one sorry",
        "THM-2069 and THM-2074 are RESERVED",
    ),
    "00-navigation/LRC14-PROOF-MAP.md": (
        "height <=2^21",
        "THM-2069/2074 remain reserved",
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

    for relative in (
        "00-navigation/SESSION-LOG.md",
        "05-knowledge/hypotheses/INDEX.md",
    ):
        if not (REPO / relative).is_file():
            errors.append(f"missing bounded current ledger: {relative}")
            continue
        validate_surface(
            relative,
            read_required(relative, errors),
            errors,
            line_budget=LINE_BUDGETS.get(relative),
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

    for value in (230, 231, 232, 233, 234, 235, 236, 237):
        if f"## MISTAKE-{value}" not in mistakes:
            errors.append(f"MISTAKES.md: current correction MISTAKE-{value} is missing")

    phase_packets = read_required(
        "01-canon/theorems/THM-2058-primitive-phase-packets-and-deck-fan-intervals.md",
        errors,
    )
    for sentinel in (
        "PROVED from THM-2053 and settled lower-dimensional LRC",
        "S_N(m)=disjoint_union_(d|N) (N/d) S_d^prim(m)",
        "bulk/boundary/null trichotomy",
        "transports by M^{-1}",
        "one explicit",
        "script_sha256: e9da9f9f5a5fdb1dc0f35814fb72b467327e76a2b90524259041e680035f7f34",
        "output_sha256: d51938376295392390408574772e19a809d8f4d2ee42bd4526bd0e35f85abadb",
    ):
        if sentinel not in phase_packets:
            errors.append(f"THM-2058 proved carrier lacks {sentinel!r}")
    if "RESERVED / UNPROVED EMPTY STUB" in phase_packets:
        errors.append("THM-2058 regressed from proved carrier to reservation")

    for theorem, sentinels in {
        "THM-2060-crt-tail-coset-saturation.md": (
            "PROVED COROLLARY / SHARP REPACKAGING",
            "q=h/gcd(N,h)=a/gcd(a,w)",
            "q-ceil(q/7)",
            "qualitative one-tail sheet dodge",
            "result_sha256: a022b5e7bb4b6b9528365836c5546c7a977e3b7b6d0e7dc6019614e6fcc8df58",
        ),
        "THM-2061-lrc14-dyadic-two-tail-folded-seam.md": (
            "PROVED REDUCTION; NOT LRC(14)",
            "G_C subset H^o_(a,b)",
            "x,y<12R",
            "sharp measure at most 4/63",
            "result_sha256: 22fb09fa81f67418d8deaf62a5a330bc7aff3928189a81e5b3586b40203370da",
        ),
        "THM-2062-two-anchor-hereditary-primitivity-crt-wheel.md": (
            "PROVED. For nonzero coefficient rows generating Z^2",
            "|B_p|<=2",
            "positive density globally and cannot",
            "result_sha256: 16ae242af27a18ccead4865edeb56b8e67113e5ad239dd24e612714a0f1d7b4f",
        ),
        "THM-2063-one-fiber-linear-planar-keller-pairs.md": (
            "PROVED. If one nonzero member of the output pencil",
            "deg_m(sP+tQ)>=2",
            "disprove JC(2).",
            "output_sha256: b164d0266bf3ec84fa897c08947a0c3463b20e7a3f2006531a199ac7ba5ac46b",
        ),
        "THM-2064-multitail-sheet-capacity-and-dyadic-seam.md": (
            "PROVED. On one safe scaled-core clock",
            "sum_(w in W) ceil(t_w/7)/t_w < 1",
            "- LRCUpTo13",
            "exact gcd of an imprimitive eleven-speed core",
            "output_sha256: 04b5e5ede02d707210400a06033061054f7e524525e22dd852602c2a33c25026",
        ),
    }.items():
        theorem_text = read_required(f"01-canon/theorems/{theorem}", errors)
        for sentinel in sentinels:
            if sentinel not in theorem_text:
                errors.append(f"{theorem}: proved result lacks {sentinel!r}")
        if "status: >\n  RESERVED" in theorem_text:
            errors.append(f"{theorem}: proved result regressed to RESERVED")

    # This is the post-THM-2064 frontier.  Keep the scope qualifiers and frozen
    # referee hashes close to the status checks: the expensive failure mode is
    # not a missing file but a proved reduction silently becoming either a
    # reservation or a full LRC(14) claim.
    for theorem, sentinels in {
        "THM-2065-two-anchor-fejer-circuit-ray-collapse.md": (
            "PROVED REDUCTION; NOT LRC(14)",
            "template without a height-bounded persistent circuit has only finitely",
            "persistent-circuit branch",
            "remains open",
            "script_sha256: 07c7e8d3a3f04372efe52678d84b8b7b2418c1897bcb15b12860025fe185e5bd",
            "result_sha256: 54a740db647a7c92fa6168e4253f0fc7357ea3491cdeee715679f3aa1c669b07",
        ),
        "THM-2066-dyadic-seam-owner-word-crt-atlas.md": (
            "PROVED. For every",
            "bitwise complementary",
            "max(C)<=24",
            "This is not LRC(14)",
            "script_sha256: 0678aaf536c60279c9989fd86cb6cc5d02ad0fd983429dc83daefa2157ee98ad",
            "output_sha256: 5d6fba079c85c532184762bcd626f78915df219f4f2b5e70a49e6ec5447c3b85",
        ),
        "THM-2067-galois-orbit-product-closes-one-variable-dvdk.md": (
            "PROVED. Let f in C[z,z^-1]",
            "the nonzero exponents of f have one",
            "strict sign",
            "project-internal proof of exactly the bare existence input used by THM-2022",
            "does not prove DvdK's stronger critical-value/limsup theorem",
            "gives no bound on the first nonzero `m`",
            "The Lean interface named `DvdK1` is not automatically discharged",
        ),
        "THM-2068-minimal-dyadic-owner-word-clock-bank.md": (
            "PROVED by exact bitset exhaustion",
            "all 59,880 primitive divisor-complete eleven-cores",
            "{25,26,27,28,32,33,34}",
            "Minimality is only within the stated",
            "twenty-clock menu",
            "script_sha256: 0ab3b5d55a70f32876908a6bd7ea6cb5ab6cb60863286eb8fc490d5809b9ebdb",
            "output_sha256: 91c67726a5f8bfd9c7982e7f8f36f57224eba61007867815911797a3f62f6b8c",
        ),
        "THM-2069-k-deletion-code-cogirth-crt-wheel.md": (
            "PROVED. Integer rows generating Z^r give, at every prime p, an injective",
            "codeword has weight at most k",
            "initial",
            "weight distribution is the exact local wheel",
            "first failure radius",
            "is matroid cogirth",
            "CRT gives the",
            "exact vector, projective, and primitive-density products",
            "rank deficiency",
            "no-finite-wheel branch",
            "proves the Paley-e8",
            "does not solve, the [72,36,16]",
            "support-realization problem or LRC(14)",
            "script_sha256: f67736997dfba6973c077c1f1b952287ddf28281c29079f7341107e71019b1a1",
            "result_sha256: 4532964632a4f7d772d5470567dc41df1935982bc236976dfbea7e3ebb89bd3f",
        ),
        "THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation.md": (
            "PROVED. Every finite one-variable Laurent polynomial embeds exactly",
            "E[P_f^m]=(Bm)! CT(f^m)",
            "f=u^2+u+u^(-1)-u^(-2)",
            "CT(f^m)=0 for every odd m",
            "CT(f^4)=-12",
            "f(-u^(-1))=-f(u)",
            "support numerical semigroup controls feasibility only",
            "THM-2067 supplies the correct algebraic",
        ),
        "THM-2071-quadratic-fiber-square-parity-gate.md": (
            "PROVED. If one member P of a complex planar Keller pair",
            "leading fiber coefficient is constant",
            "reduced fiber degree one",
            "explicit tame normal form",
            "centered parity decomposition and a noncancellation law",
            "no hypothetical JC(2)",
            "counterexample has a quadratic member",
            "This is a pencil-direction theorem, not a statement about generic",
            "script_sha256: 3d5ce81db8601a3035db28ae63bf3d003d4f72372d304c3250098a64a8efb267",
            "output_sha256: 2b692ddf3606a0226bb22b88b0a2060a5ecff951b41ba02a8dee3af111821030",
        ),
        "THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate.md": (
            "PROVED. No finite THM-2066 owner-word clock bank fixed independently",
            "all sampled safe",
            "packets empty and every owner constraint vacuous",
            "if the closed weak-safe core contains a phase together with its half-translate",
            "certificate-strategy theorem, not LRC(14)",
        ),
        "THM-2073-lrc14-dyadic-deletion-tower.md": (
            "PROVED REDUCTION; NOT LRC(14)",
            "after at most eight levels",
            "one speed of each valuation below the terminal",
            "hereditarily primitive terminal core",
            "This does **not** prove LRC(14)",
            "script_sha256: a65a72a538897d3d9b2f6a25ce8c9099e88679fe8292b2e7213d7816e5b095e4",
            "output_sha256: c99d6dcfb2d0aae805492e2a8349ef6a1793eda4494cf54fb440aedac4d79e24",
        ),
        "THM-2075-safe-child-homeomorphism-and-wall-word-conjugacy.md": (
            "PROVED from THM-2073's unique safe-child law",
            "Doubling is a homeomorphism",
            "Component count",
            "Euler characteristic, and endpoint count are invariant",
            "every endpoint is inherited from a",
            "terminal-core speed",
            "This is a carrier theorem, not LRC(14)",
        ),
        "THM-2076-guard-capacity-terminal-rank-floor.md": (
            "PROVED by Haar measure and compact/open separation",
            "forces s>=6",
            "depth at most five",
            "measure strictly below",
            "2^(1-r)/7",
            "This is a structural reduction, not LRC(14)",
        ),
        "THM-2077-terminal-kakeya-needle-and-recursive-quarter-escape.md": (
            "PROVED. THM-2075 lifts a terminal maximizer interval",
            "uniform relative-height bounds",
            "folded",
            "owner-determinant sidecar",
            "every quotient level must fail",
            "ratio escape at every depth",
            "This is a structural reduction, not LRC(14)",
        ),
        "THM-2078-bounded-terminal-guard-containment-closure.md": (
            "PROVED by exact integer-bitset census",
            "terminal maximum at most 24",
            "4,484,931 cores",
            "exactly 30,594",
            "denominator 8192",
            "zero pairs survive",
            "not LRC(14)",
            "script_sha256: 4e0d36e38aa2fabf74bf7ba012f24c764e95cc79dd44e4a538d12814b318de53",
            "output_sha256: b542a94aebec338ca4b91253b205429f986b881b93f276ac0036237e1bd888d2",
        ),
        "THM-2079-mirror-complement-safe-child-addresses.md": (
            "PROVED from THM-2073/2075",
            "Reversal acts freely",
            "mirror components have complementary addresses",
            "a and 2^r-1-a",
            "original odd tail flips its nearest-integer owner parity",
            "covering condition is equivariant, not contradicted by mirror",
            "LRC(14) is not proved",
        ),
    }.items():
        theorem_text = read_required(f"01-canon/theorems/{theorem}", errors)
        for sentinel in sentinels:
            if sentinel not in theorem_text:
                errors.append(f"{theorem}: current scoped result lacks {sentinel!r}")
        if "status: >\n  RESERVED" in theorem_text:
            errors.append(f"{theorem}: proved result regressed to RESERVED")

    for theorem, sentinels in {
        "THM-2074-lrc14-density-one-relation-hyperplane-sieve.md": (
            "RESERVED WITH COMPLETE PROOF UNDER AUDIT",
            "This is an almost-everywhere theorem, not LRC(14)",
            "This ID reserves the counting and finite-field consequences",
            "will replace this stub after audit",
        ),
    }.items():
        theorem_text = read_required(f"01-canon/theorems/{theorem}", errors)
        for sentinel in sentinels:
            if sentinel not in theorem_text:
                errors.append(f"{theorem}: reservation lacks {sentinel!r}")
        if not re.search(r"^status:\s*>\n\s+RESERVED\b", theorem_text, re.MULTILINE):
            errors.append(f"{theorem}: reservation was promoted without a truth-surface update")

    gmc2_paper = read_required(
        "01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md",
        errors,
    )
    for sentinel in (
        "PROVED. For every finite exact support",
        "project-internal Galois orbit-product theorem THM-2067",
        "proves the specialized normalized relation is both",
        "zero and nonzero, then derives `NC2` and GMC(2) from `DvdK1` and the compact",
        "`HeightWitnessSupplier`",
        "THM-2067 gives exactly the bare existence statement needed",
        "one proved-but-unformalized proposition and one named composition interface",
    ):
        if sentinel not in gmc2_paper:
            errors.append(f"THM-2022 formalization ledger lacks {sentinel!r}")
    if "The sole internal Lean boundary" in gmc2_paper:
        errors.append(
            "THM-2022: stale singular Lean-boundary summary hides the separate DvdK1 interface"
        )

    lean_formalization_files = {
        "04-computation/lean/TournamentH7/TournamentH7/GMC2DvdKInterface.lean": (
            "def DvdK1 : Prop :=",
            "All downstream theorems take a proof of `DvdK1` explicitly",
            "theorem exists_nonzero_face_seed",
        ),
        "04-computation/lean/TournamentH7/TournamentH7/GMC2NC2.lean": (
            "structure HeightWitnessSupplier : Prop where",
            "theorem nc2_of_dvdK1_of_heightWitnessSupplier",
            "(hDvdK : GMC2DvdKInterface.DvdK1)",
            "(hHeight : HeightWitnessSupplier)",
            "theorem gmc2_of_dvdK1_of_heightWitnessSupplier",
        ),
        "04-computation/lean/TournamentH7/TournamentH7/GMC2NC2Capstone.lean": (
            "former work-in-progress theorem",
            "sorry-free compatibility",
            "It deliberately exposes both remaining inputs",
            "theorem nc2_of_dvdk1_of_heightWitnessSupplier",
            "theorem gmc2_of_dvdk1_of_heightWitnessSupplier",
        ),
        "04-computation/lean/TournamentH7/TournamentH7/GMC2DvdKTwoCharge.lean": (
            "Elementary single-character (two-charge) DvdK, DvdK-premise-free",
            "lemma balanced_unique",
            "theorem exists_nonzero_ct_pair",
            "theorem dvdk1_pair",
            "Swapped orientation",
            "theorem exists_nonzero_ct_pair'",
            "#print axioms GMC2DvdKTwoCharge.dvdk1_pair",
        ),
        "04-computation/lean/TournamentH7/TournamentH7/GMC2Formalization.lean": (
            "import TournamentH7.GMC2NC2",
            "derives both `NC2` and GMC(2) from `DvdK1`",
            "plus a compact `HeightWitnessSupplier`",
            "this interface therefore stays explicit",
        ),
    }
    for relative, sentinels in lean_formalization_files.items():
        lean_text = read_required(relative, errors)
        for sentinel in sentinels:
            if sentinel not in lean_text:
                errors.append(f"{relative}: formal boundary lacks {sentinel!r}")
        if re.search(r"^\s*(?:by\s+)?sorry\b", lean_text, re.MULTILINE):
            errors.append(f"{relative}: executable sorry returned to the checked capstone")

    for relative in (
        "README.md",
        "00-navigation/START-HERE.md",
        "00-navigation/CURRENT-FRONTIER.md",
        "01-canon/ACTIVE-GUARDRAILS.md",
        "00-navigation/SESSION-LOG.md",
        "05-knowledge/hypotheses/INDEX.md",
    ):
        formal_summary = maintained_text(relative, errors)
        for sentinel in ("DvdK1", "HeightWitnessSupplier"):
            if sentinel not in formal_summary:
                errors.append(f"{relative}: current Lean summary omits {sentinel}")
        for retired in (
            "exactly one sorry remains",
            "exactly one `sorry` remains",
            "one remaining sorry",
            "one remaining `sorry`",
        ):
            if retired.casefold() in formal_summary.casefold():
                errors.append(f"{relative}: obsolete one-sorry summary returned: {retired!r}")

    corrected_jc = read_required(
        "05-knowledge/hypotheses/"
        "HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs.md",
        errors,
    )
    for sentinel in (
        "MIXED / CORRECTED BY MISTAKE-237",
        "det Hess(P)",
        "-4 d^2(d-1)^2 A B",
        "The symmetric-Jacobian reduction used for a general planar Keller map changes",
        "Three live programs, not three equivalent forms",
        "fiber-degree descent",
    ):
        if sentinel not in corrected_jc:
            errors.append(f"HYP-8905 corrected map lacks {sentinel!r}")

    for relative in (
        "04-computation/planar_jc_is_nc2_one_sidedness_deathstar_S103.py",
        "05-knowledge/results/planar_jc_is_nc2_one_sidedness_deathstar_S103.out",
        "04-computation/jc2_via_vanishing_conjecture_and_the_cf_termination_boxeph_S225.py",
        "05-knowledge/results/jc2_via_vanishing_conjecture_and_the_cf_termination_boxeph_S225.out",
    ):
        corrected_artifact = read_required(relative, errors)
        if "MISTAKE-237" not in corrected_artifact or "HYP-8905" not in corrected_artifact:
            errors.append(f"{relative}: corrected JC artifact lacks its audit route")
        for retired in (
            "THE UNIFICATION (JC <-> GMC2",
            "JC(3) FALSE",
            "Alpoge's Keller counterexample",
            "JC(2), LRC(14), GMC(2) ALL reduce",
        ):
            if retired in corrected_artifact:
                errors.append(f"{relative}: retracted JC bridge returned: {retired!r}")

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
    if not hypotheses.startswith("> **CURRENT DIGEST"):
        errors.append("hypotheses/INDEX.md: current digest must remain first")
    if "INDEX-HISTORICAL-THROUGH-2026-07-21.md" not in hypotheses:
        errors.append("hypotheses/INDEX.md: split historical ledger is not routed")
    if "# Hypothesis Log — Index" in hypotheses:
        errors.append("hypotheses/INDEX.md: historical suffix leaked into bounded router")
    current_digest = hypotheses
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
    for value in (8878, 8879, 8885, 8890, 8895, 8900, 8905):
        if f"HYP-{value}" not in current_digest:
            errors.append(
                f"hypotheses/INDEX.md: current digest lacks incoming HYP-{value} routing"
            )
    if "CRT PACKET JOIN (PROVED; not LRC(14))" not in current_digest:
        errors.append("hypotheses/INDEX.md: current digest lacks THM-2059 routing")
    for sentinel in (
        "PRIMITIVE PACKET AND OWNER INTERVAL (PROVED; not LRC(14))",
        "SHARP TAIL COSETS (PROVED; not LRC(14))",
        "COMMON-FIBER CAPACITY (PROVED; not LRC(14))",
        "DYADIC SEAM (PROVED REDUCTION; not LRC(14))",
        "HEREDITARY CRT WHEEL (PROVED; not LRC(14))",
        "CIRCUIT-RAY COLLAPSE (PROVED REDUCTION; not LRC(14))",
        "OWNER-WORD ATLAS (PROVED / FINITE-EXACT)",
        "BOUNDED-MENU MINIMUM (PROVED / FINITE-EXACT)",
        "FIXED-BANK NO-GO (PROVED)",
        "DYADIC DELETION TOWER (PROVED REDUCTION; not LRC(14))",
        "CODE-COGIRTH WHEEL (PROVED; not LRC(14))",
        "DENSITY-ONE SIEVE (RESERVED)",
        "ONE-FIBER-LINEAR KELLER PAIRS (PROVED; JC(2) OPEN)",
        "QUADRATIC-FIBER GATE (PROVED; JC(2) OPEN)",
        "ONE-VARIABLE CT EXISTENCE (PROVED)",
        "SEMIGROUP BYPASS CORRECTION (PROVED)",
        "unbounded hereditarily primitive terminal, which remains open",
        "`DvdK1` and `HeightWitnessSupplier` interfaces remain",
        "HYP-8879-lrc-gmc-weighted-fiber-analogy-corrected.md",
        "HYP-8895-return-semigroup-reachability-not-noncancellation.md",
        "HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs.md",
    ):
        if sentinel not in current_digest:
            errors.append(f"hypotheses/INDEX.md: current digest lacks {sentinel!r}")

    session_log = read_required("00-navigation/SESSION-LOG.md", errors)
    if not session_log.startswith("> **CURRENT-TRUTH WARNING"):
        errors.append("SESSION-LOG.md: current-truth warning must remain first")
    if "SESSION-LOG-HISTORICAL-THROUGH-2026-07-21.md" not in session_log:
        errors.append("SESSION-LOG.md: split historical ledger is not routed")
    for sentinel in (
        "THM-2065 PROVED REDUCTION",
        "THM-2066 PROVED / FINITE-EXACT",
        "THM-2068 PROVED / FINITE-EXACT",
        "THM-2072 PROVED",
        "THM-2073 PROVED REDUCTION",
        "THM-2069 PROVED",
        "THM-2074 is RESERVED",
        "THM-2022 proves NC2 and GMC(2) on paper",
        "`DvdK1` plus",
        "`HeightWitnessSupplier`",
        "Do not summarize the",
        "current state as “one sorry” or as a complete formal proof",
        "THM-2070 PROVED CORRECTION",
        "THM-2071 PROVED",
    ):
        if sentinel not in session_log:
            errors.append(f"SESSION-LOG.md: current digest lacks {sentinel!r}")

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

    compact_lrc_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "LRC14", "--recent", "1", "--max-matches", "12",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    if "CURRENT-FRONTIER.md:12: ## LRC(14)" not in compact_lrc_smoke.stdout:
        errors.append("start_session.py: compact LRC14 did not route to LRC(14)")

    for relative, expected in (
        (
            "01-canon/theorems/"
            "THM-1430-explicit-symmetric-keller-counterexample-on-C6.md",
            "VERIFIED",
        ),
        (
            "01-canon/theorems/"
            "THM-1990-the-figurate-reciprocal-ladder-and-the-harmonic-edge.md",
            "PROVED",
        ),
    ):
        actual = file_status(relative)
        if actual != expected:
            errors.append(
                f"start_session.py: folded status for {relative} is "
                f"{actual}, expected {expected}"
            )

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

    jc_bridge_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "binary symmetric Hessian NC2 GMC JC2",
            "--recent", "1", "--max-matches", "12",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    if "MISTAKE-237" not in jc_bridge_smoke.stdout:
        errors.append("start_session.py: symmetric-JC topic did not surface MISTAKE-237")
    if (
        "[MIXED] 05-knowledge/hypotheses/"
        "HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs.md"
        not in jc_bridge_smoke.stdout
    ):
        errors.append("start_session.py: symmetric-JC topic lacks corrected HYP-8905")

    packet_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "THM-2058 primitive phase packets",
            "--recent", "1", "--max-matches", "16",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    packet = packet_smoke.stdout
    if (
        "[PROVED] 01-canon/theorems/THM-2058-primitive-phase-packets-and-deck-fan-intervals.md"
        not in packet
    ):
        errors.append("start_session.py: THM-2058 is not routed as proved canon")
    reservation_block = packet.partition("  Reservations (not results):")[2].partition(
        "  Hypotheses:"
    )[0]
    if "THM-2058" in reservation_block:
        errors.append("start_session.py: THM-2058 leaked into reservations")

    for theorem, topic in (
        ("THM-2060-crt-tail-coset-saturation.md", "THM-2060 tail coset saturation"),
        ("THM-2061-lrc14-dyadic-two-tail-folded-seam.md", "THM-2061 dyadic seam"),
        ("THM-2062-two-anchor-hereditary-primitivity-crt-wheel.md", "THM-2062 hereditary CRT wheel"),
        ("THM-2063-one-fiber-linear-planar-keller-pairs.md", "THM-2063 one fiber linear Keller"),
        ("THM-2064-multitail-sheet-capacity-and-dyadic-seam.md", "THM-2064 multi tail capacity"),
        ("THM-2065-two-anchor-fejer-circuit-ray-collapse.md", "THM-2065 circuit ray collapse"),
        ("THM-2066-dyadic-seam-owner-word-crt-atlas.md", "THM-2066 owner word atlas"),
        ("THM-2067-galois-orbit-product-closes-one-variable-dvdk.md", "THM-2067 Galois orbit product"),
        ("THM-2068-minimal-dyadic-owner-word-clock-bank.md", "THM-2068 minimal clock bank"),
        ("THM-2069-k-deletion-code-cogirth-crt-wheel.md", "THM-2069 deletion code cogirth wheel"),
        ("THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation.md", "THM-2070 Wick return cancellation"),
        ("THM-2071-quadratic-fiber-square-parity-gate.md", "THM-2071 quadratic fiber parity gate"),
        ("THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate.md", "THM-2072 fixed owner clock no go"),
        ("THM-2073-lrc14-dyadic-deletion-tower.md", "THM-2073 dyadic deletion tower"),
        ("THM-2075-safe-child-homeomorphism-and-wall-word-conjugacy.md", "THM-2075 safe child homeomorphism"),
        ("THM-2076-guard-capacity-terminal-rank-floor.md", "THM-2076 guard capacity terminal rank"),
        ("THM-2077-terminal-kakeya-needle-and-recursive-quarter-escape.md", "THM-2077 terminal needle quarter escape"),
        ("THM-2078-bounded-terminal-guard-containment-closure.md", "THM-2078 bounded terminal guard closure"),
        ("THM-2079-mirror-complement-safe-child-addresses.md", "THM-2079 mirror complement addresses"),
    ):
        smoke = subprocess.run(
            (
                sys.executable, str(REPO / "agents/start_session.py"),
                "--topic", topic, "--recent", "1", "--max-matches", "16",
            ),
            cwd=REPO, text=True, capture_output=True, check=False,
        )
        proposal_packet = smoke.stdout
        expected = f"[PROVED] 01-canon/theorems/{theorem}"
        if expected not in proposal_packet:
            errors.append(f"start_session.py: {theorem} is not routed as proved")
        canon_block = proposal_packet.partition("  Canon:")[2].partition(
            "  Reservations (not results):"
        )[0]
        if theorem not in canon_block:
            errors.append(f"start_session.py: {theorem} is missing from proved canon")

    for theorem, topic in (
        (
            "THM-2074-lrc14-density-one-relation-hyperplane-sieve.md",
            "THM-2074 density one hyperplane sieve",
        ),
    ):
        reserved_smoke = subprocess.run(
            (
                sys.executable, str(REPO / "agents/start_session.py"),
                "--topic", topic, "--recent", "1", "--max-matches", "16",
            ),
            cwd=REPO, text=True, capture_output=True, check=False,
        )
        expected = f"[RESERVED] 01-canon/theorems/{theorem}"
        if expected not in reserved_smoke.stdout:
            errors.append(f"start_session.py: {theorem} is not routed as reserved")
        reserved_canon = reserved_smoke.stdout.partition("  Canon:")[2].partition(
            "  Reservations (not results):"
        )[0]
        if theorem in reserved_canon:
            errors.append(f"start_session.py: {theorem} leaked into proved canon")

    hypothesis_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "HYP-8846 HYP-8871 HYP-8879 HYP-8895 HYP-8905",
            "--recent", "1", "--max-matches", "24",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    hypothesis_packet = hypothesis_smoke.stdout
    for sentinel in (
        "[OPEN] 05-knowledge/hypotheses/HYP-8846-lrc14-pointed-plane-transport.md",
        "[OPEN] 05-knowledge/hypotheses/HYP-8871-lrc14-owner-sector-klein-sail-automaton.md",
        "[PARTIAL] 05-knowledge/hypotheses/HYP-8879-lrc-gmc-weighted-fiber-analogy-corrected.md",
        "[PARTIAL] 05-knowledge/hypotheses/HYP-8895-return-semigroup-reachability-not-noncancellation.md",
        "[MIXED] 05-knowledge/hypotheses/HYP-8905-binary-symmetric-jc2-subcase-and-separate-descent-programs.md",
    ):
        if sentinel not in hypothesis_packet:
            errors.append("start_session.py: folded-status routing lacks " + repr(sentinel))

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
