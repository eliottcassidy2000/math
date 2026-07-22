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
        "THM-2058 proves the primitive",
        "THM-2059 the arbitrary-clock CRT join",
        "2064 reduce imprimitive two-tail capacity",
        "THM-2062 adds a hereditary CRT wheel",
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
        "THM-2057 closes two AP one-tail families",
        "THM-2058 splits reduced-order packets",
        "THM-2059 joins",
        "THM-2060/2064 reduce imprimitive two-tail capacity",
        "THM-2062 adds the hereditary",
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
        "THM-2065's proposed circuit-ray collapse is",
        "TournamentH7.GMC2Formalization",
        "TournamentH7.GMC2NC2Capstone",
        "unique minimal balanced channel",
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
        "THM-2065 is RESERVED",
        "Fiber size is observable-relative Hartley ambiguity",
        "Clocks are not modular cusps",
        "The binary symmetric Hessian is not full JC(2)",
        "THM-1330 is a necessary atlas, not a classification",
    ),
    "00-navigation/LRC14-PROOF-MAP.md": (
        "## 2026-07-21 current control panel",
        "bounded support-at-most-3 code of rank >=11",
        "transverse deck D_N(m)>=1/14",
        "finite Farey address; listed rays unresolved",
        "scaled AP one-tail clock/binding leaves",
        "tail-sheet capacity -> dyadic seam",
        "hereditary-primitivity CRT wheel",
        "### Mandatory hostile controls",
    ),
    "00-navigation/LRC-TECHNIQUE-INDEX.md": (
        "SEARCHABLE ATLAS, NOT STARTUP TRUTH",
        "## 2026-07-21 current-use overlay",
        "## LTI-532 - Rank-eleven relation-code dispatcher",
        "Transverse deck",
        "Missing-clock binding",
        "CRT packet join",
        "THM-2060/2064 reduce tail capacity",
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
    ),
    "00-navigation/START-HERE.md": (
        "every counterexample has a 2--5-term integer relation of height at most `2^20`",
        "height at most `2^21`",
        "The remaining jump is pointed Euler survival",
        "They are not root-imported",
        "THM-2058 is an unproved empty stub",
        "THM-2060/2061 are reserved proof",
    ),
    "00-navigation/CURRENT-FRONTIER.md": (
        "but are not root-imported",
        "6. **Pointed plane transport.**",
        "The incoming THM-2054 relative-Fejer program",
        "THM-2058 supplies the primitive phase-packet",
        "THM-2058 is an empty reservation",
        "Reserved proof candidates—not theorem inputs",
    ),
    "01-canon/ACTIVE-GUARDRAILS.md": (
        "height-`2^21` relation",
        "MISTAKE-225 corrects HYP-8865",
        "are independent reserved proof candidates",
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

    circuit_candidate = read_required(
        "01-canon/theorems/THM-2065-two-anchor-fejer-circuit-ray-collapse.md",
        errors,
    )
    for sentinel in (
        "RESERVED WITH A COMPLETE PROOF UNDER AUDIT",
        "No consumer may cite this reservation as a theorem",
        "support-three-to-five integer relation",
    ):
        if sentinel not in circuit_candidate:
            errors.append(f"THM-2065 reservation lacks {sentinel!r}")
    if re.search(r"^status:\s*>?\s*PROVED\b", circuit_candidate, re.MULTILINE):
        errors.append("THM-2065 was promoted without updating startup truth")

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
        "CIRCUIT-RAY COLLAPSE (RESERVED; PROOF UNDER AUDIT)",
        "ONE-FIBER-LINEAR KELLER PAIRS (PROVED; JC(2) OPEN)",
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

    reserved_smoke = subprocess.run(
        (
            sys.executable, str(REPO / "agents/start_session.py"),
            "--topic", "THM-2065 circuit ray collapse", "--recent", "1",
            "--max-matches", "16",
        ),
        cwd=REPO, text=True, capture_output=True, check=False,
    )
    if (
        "[RESERVED] 01-canon/theorems/THM-2065-two-anchor-fejer-circuit-ray-collapse.md"
        not in reserved_smoke.stdout
    ):
        errors.append("start_session.py: THM-2065 is not routed as reserved")

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
