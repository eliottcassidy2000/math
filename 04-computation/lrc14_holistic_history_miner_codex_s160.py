#!/usr/bin/env python3
"""Mine the repo's LRC14 history into route families and guardrails.

This is not a proof script.  It is a reproducible synthesis aid for the
current LRC14 proof search: it scans navigation, hypothesis, result-index, and
reflection files; groups recurring ideas into proof-route families; extracts
refutation/guardrail snippets; and emits a tournament analysis on route
families.  Vertices are proof carriers rather than runners or arcs.
"""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class Family:
    key: str
    label: str
    tier: int
    terms: tuple[str, ...]


FAMILIES = (
    Family(
        "quotient_tournament",
        "Unlabelled Quotients And Tournament Shadows",
        10,
        (
            "a000568",
            "tournament",
            "chamber",
            "isomorphism",
            "apex-pressure",
            "fixed-round",
            "source-fiber",
            "unlabelled",
            "marked chamber",
        ),
    ),
    Family(
        "endpoint_fold",
        "Endpoint/Fold Certificate Calculus",
        20,
        (
            "endpoint",
            "owner",
            "fold",
            "pincer",
            "boundary",
            "circuit",
            "pair-sum",
            "certificate calculus",
            "exposed cell",
        ),
    ),
    Family(
        "qdiv_farey_haar",
        "Exact q/Farey/Haar-Baire Fronts",
        30,
        (
            "qdiv",
            "q_threshold",
            "farey",
            "exact m",
            "haar",
            "baire",
            "strict-open",
            "safe mass",
            "boundary-only",
            "rational interval",
        ),
    ),
    Family(
        "c27_unital_k33",
        "C27/Unital/K33 State-Lift Labels",
        40,
        (
            "c27",
            "unital",
            "k33",
            "state-lift",
            "state lift",
            "h=7",
            "d9",
            "shell transfer",
            "unit petal",
        ),
    ),
    Family(
        "wide_gk8_relation",
        "Wide/Relation-Lattice/gK8 Moment Control",
        45,
        (
            "wide",
            "gk8",
            "l_y",
            "delsarte",
            "p0",
            "cap",
            "bonferroni",
            "survival",
            "miss-zeta",
            "residual profile",
            "doublet",
            "tornheim",
            "relation lattice",
        ),
    ),
    Family(
        "labelled_packets",
        "Labelled Packets, Fixed Margins, Gauntlets",
        55,
        (
            "packet",
            "fixed-margin",
            "gauntlet",
            "family",
            "sporadic",
            "source-spectrum",
            "labelled",
            "f6",
            "f7",
            "johnson",
        ),
    ),
    Family(
        "apex_lift_comb",
        "Moon-Core Apex/Lift/Comb Packets",
        58,
        (
            "moon",
            "apex",
            "aperture",
            "comb",
            "14z",
            "lift-packet",
            "few-apex",
            "low-multiple",
            "multiple",
        ),
    ),
    Family(
        "boundary_moment_nork",
        "Boundary-Moment, NORK, Pinch Templates",
        62,
        (
            "nork",
            "pinch",
            "boundary-moment",
            "covering moment",
            "endpoint-owner",
            "zero-open",
            "bridge",
            "gk8/l_y",
        ),
    ),
    Family(
        "dual_certificates",
        "Endpoint-Credit, Twist, Count, Toeplitz Duals",
        64,
        (
            "toeplitz",
            "fourier",
            "twist",
            "danger-count",
            "multiplicity",
            "endpoint-credit",
            "farkas",
            "polynomial dual",
            "psd",
            "potential",
            "dual certificate",
        ),
    ),
    Family(
        "formal_finite",
        "Formal/Finite-Check Proof Interfaces",
        70,
        (
            "lean",
            "finite check",
            "finite-check",
            "periodmax",
            "verified",
            "exhaustive",
            "bounded",
            "part a",
            "theorem",
            "thm-",
        ),
    ),
    Family(
        "external_carriers",
        "External Carriers And Analogy Reservoir",
        15,
        (
            "irreducible",
            "polynomial",
            "bunyakovsky",
            "faulhaber",
            "triangular",
            "pollock",
            "unit distance",
            "pentagonal",
            "euler",
            "pisano",
            "tiling",
            "tetrahedral",
            "flower",
        ),
    ),
)


GUARDRAIL_TERMS = (
    "false",
    "refut",
    "dead",
    "guardrail",
    "too coarse",
    "not enough",
    "fails",
    "fail",
    "blind",
    "invisible",
    "wrong",
    "cannot",
    "overclaim",
    "mistake",
    "insufficient",
)


RETENTION_TERMS = {
    "exact": 3,
    "label": 3,
    "labelled": 3,
    "owner": 2,
    "source-spectrum": 4,
    "fixed-margin": 3,
    "qdiv": 2,
    "farey": 2,
    "haar": 2,
    "baire": 2,
    "c27": 2,
    "unital": 2,
    "k33": 2,
    "endpoint": 2,
    "moment": 2,
    "dual": 3,
    "toeplitz": 3,
    "state-lift": 3,
    "farkas": 2,
}


def iter_markdown_files() -> list[Path]:
    files: list[Path] = []
    for name in (
        "00-navigation/CONCEPT-MAP.md",
        "00-navigation/OPEN-QUESTIONS.md",
        "00-navigation/SESSION-LOG.md",
        "00-navigation/TANGENTS.md",
        "05-knowledge/hypotheses/INDEX.md",
        "05-knowledge/results/INDEX.md",
    ):
        p = ROOT / name
        if p.exists():
            files.append(p)
    for base in (ROOT / "05-knowledge/hypotheses", ROOT / "07-reflections"):
        if not base.exists():
            continue
        for p in sorted(base.glob("*.md")):
            low = p.name.lower()
            if (
                "lrc" in low
                or "lonely-runner" in low
                or "unit-distance" in low
                or "pollock" in low
                or "tournament" in low
                or "faulhaber" in low
                or "triangular" in low
            ):
                files.append(p)
    return sorted(dict.fromkeys(files))


def norm(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


def family_hits(text: str) -> Counter[str]:
    low = text.lower()
    hits: Counter[str] = Counter()
    for fam in FAMILIES:
        for term in fam.terms:
            hits[fam.key] += low.count(term)
    return hits


def retention_score(text: str) -> int:
    low = text.lower()
    return sum(low.count(term) * weight for term, weight in RETENTION_TERMS.items())


def split_entries(path: Path, text: str) -> list[tuple[str, str]]:
    """Return rough logical entries with a stable label."""
    rel = path.relative_to(ROOT).as_posix()
    if path.name == "INDEX.md":
        entries = []
        for line in text.splitlines():
            if line.startswith("- **") or line.startswith("| **"):
                entries.append((f"{rel}:{len(entries) + 1}", line))
        return entries or [(rel, text)]
    chunks = re.split(r"\n(?=#{1,3} )", text)
    if len(chunks) <= 1:
        return [(rel, text)]
    out = []
    for i, chunk in enumerate(chunks, start=1):
        if chunk.strip():
            first = chunk.strip().splitlines()[0].lstrip("# ").strip()
            out.append((f"{rel}:{i}:{first[:70]}", chunk))
    return out


def extract_ids(text: str) -> Counter[str]:
    ids = Counter()
    for match in re.findall(r"\b(?:HYP|THM|OPEN-Q|T)\-?\d+\b", text):
        ids[match] += 1
    return ids


def top_items(counter: Counter[str], n: int = 8) -> str:
    if not counter:
        return "(none)"
    return ", ".join(f"{k}:{v}" for k, v in counter.most_common(n))


def main() -> None:
    files = iter_markdown_files()
    entry_count = 0
    family_entry_hits: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    family_file_hits: dict[str, set[str]] = defaultdict(set)
    guardrails: dict[str, list[tuple[str, str]]] = defaultdict(list)
    seen_guardrails: set[tuple[str, str]] = set()
    id_counts: Counter[str] = Counter()
    cross_terms: Counter[str] = Counter()
    co_occurrence: Counter[tuple[str, str]] = Counter()

    for path in files:
        text = path.read_text(encoding="utf-8", errors="ignore")
        id_counts.update(extract_ids(text))
        entries = split_entries(path, text)
        for label, entry in entries:
            entry_count += 1
            hits = family_hits(entry)
            active = [key for key, count in hits.items() if count > 0]
            for key in active:
                family_entry_hits[key].append((label, hits[key], retention_score(entry)))
                family_file_hits[key].add(path.relative_to(ROOT).as_posix())
            for i, left in enumerate(active):
                for right in active[i + 1 :]:
                    co_occurrence[tuple(sorted((left, right)))] += 1

            for line_no, line in enumerate(entry.splitlines(), start=1):
                low_line = line.lower()
                if not any(term in low_line for term in GUARDRAIL_TERMS):
                    continue
                line_hits = family_hits(line)
                line_active = [key for key, count in line_hits.items() if count > 0]
                if not line_active:
                    line_active = active
                clean = norm(line)
                for key in line_active:
                    marker = (key, clean)
                    if marker in seen_guardrails:
                        continue
                    seen_guardrails.add(marker)
                    if len(guardrails[key]) < 12:
                        guardrails[key].append((f"{label}:L{line_no}", clean[:260]))

            for term in FAMILIES[-1].terms:
                cross_terms[term] += entry.lower().count(term)

    print("LRC14 HOLISTIC HISTORY MINER")
    print("============================")
    print(f"repo_root: {ROOT}")
    print(f"files_scanned: {len(files)}")
    print(f"logical_entries_scanned: {entry_count}")
    print(f"id_top: {top_items(id_counts, 12)}")
    print()

    print("ROUTE FAMILY CENSUS")
    print("-------------------")
    family_scores: dict[str, float] = {}
    for fam in FAMILIES:
        entries = family_entry_hits.get(fam.key, [])
        hit_sum = sum(count for _, count, _ in entries)
        file_count = len(family_file_hits.get(fam.key, set()))
        retention = sum(score for _, _, score in entries)
        norm_score = fam.tier + (retention / math.sqrt(1 + max(file_count, 1))) / 1000
        family_scores[fam.key] = norm_score
        print(f"{fam.label}")
        print(f"  term_hits={hit_sum} files={file_count} entries={len(entries)}")
        print(f"  retention_weight={retention} tiered_score={norm_score:.3f}")
        for label, count, score in sorted(entries, key=lambda x: (-x[1], -x[2], x[0]))[:5]:
            print(f"  top_entry: hits={count} retention={score} {label}")
        print()

    print("REFUTATION / GUARDRAIL SNIPPETS BY FAMILY")
    print("-----------------------------------------")
    for fam in FAMILIES:
        snippets = guardrails.get(fam.key, [])
        print(f"{fam.label}: {len(snippets)} sampled guardrail entries")
        for label, snippet in snippets[:4]:
            print(f"  - {label}: {snippet}")
        print()

    print("CROSS-TOPIC RESERVOIR COUNTS")
    print("----------------------------")
    print(top_items(cross_terms, len(cross_terms)))
    print()

    print("CO-OCCURRENCE SPINE")
    print("-------------------")
    for (left, right), count in co_occurrence.most_common(18):
        l_lab = next(f.label for f in FAMILIES if f.key == left)
        r_lab = next(f.label for f in FAMILIES if f.key == right)
        print(f"{count:4d}  {l_lab}  <->  {r_lab}")
    print()

    print("TOURNAMENT ANALYSIS ON PROOF-ROUTE FAMILIES")
    print("-------------------------------------------")
    print("vertices: proof-route families, not runners, arcs, or raw residues")
    print(
        "pairwise_observable: which family keeps more of the LRC predicate payload "
        "(exact scale, open/boundary status, source labels, owner labels, and "
        "dual obstruction data)"
    )
    print(
        "switch/gauge: higher tiered_score dominates; ties use the fixed "
        "historical Hamiltonian path listed below"
    )
    hamiltonian = [fam.key for fam in sorted(FAMILIES, key=lambda f: f.tier)]
    print(
        "tie_hamiltonian_path: "
        + " -> ".join(next(f.label for f in FAMILIES if f.key == key) for key in hamiltonian)
    )

    scores = Counter()
    edges: list[tuple[str, str]] = []
    for i, a in enumerate(FAMILIES):
        for b in FAMILIES[i + 1 :]:
            if (family_scores[a.key], -a.tier) > (family_scores[b.key], -b.tier):
                winner, loser = a, b
            else:
                winner, loser = b, a
            edges.append((winner.key, loser.key))
            scores[winner.key] += 1
            scores.setdefault(loser.key, scores.get(loser.key, 0))

    score_hist = Counter(scores.values())
    print(f"score_histogram: {dict(sorted(score_hist.items()))}")
    print("top_scores:")
    for key, score in scores.most_common():
        print(f"  {score:2d}  {next(f.label for f in FAMILIES if f.key == key)}")

    # The tiered-score rule is acyclic by construction except exact equalities;
    # test 3-cycles anyway as a guard against accidental orientation mistakes.
    edge_set = set(edges)
    c3 = 0
    keys = [f.key for f in FAMILIES]
    for i, a in enumerate(keys):
        for j, b in enumerate(keys[i + 1 :], start=i + 1):
            for c in keys[j + 1 :]:
                if (
                    ((a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set)
                    or ((a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set)
                ):
                    c3 += 1
    print(f"directed_3_cycles: {c3}")
    print()

    print("SYNTHESIS READOUT")
    print("-----------------")
    print(
        "The mined history supports a one-way refinement pattern: scalar or "
        "unlabelled quotients are useful scouts, but every serious obstruction "
        "forces more payload to be retained.  The durable current object is a "
        "labelled packet sheaf over exact q/Farey/Haar fronts, C27/unital/K33 "
        "owner data, endpoint/boundary-moment charts, and several dual cones.  "
        "AP/Goddyn-Wong survive as common equality atoms; the live proof task is "
        "to show that no qdiv>14 zero-open packet can remain invisible to all "
        "packet labels and all dual certificates."
    )


if __name__ == "__main__":
    main()
