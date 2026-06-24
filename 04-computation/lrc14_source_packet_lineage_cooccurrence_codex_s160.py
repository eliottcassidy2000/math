#!/usr/bin/env python3
"""Mine the LRC14 research log for route families and quotient guardrails.

This is a synthesis aid for HYP-2976.  It does not recompute lonely-runner
values.  Instead it scans the persistent hypothesis/result/forum history and
asks which proof-object fibers each route keeps before scalarizing.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from itertools import combinations
from pathlib import Path
import re


ROOTS = [
    Path("05-knowledge/hypotheses"),
    Path("05-knowledge/results"),
    Path("07-reflections"),
    Path("poke-forum/posts"),
    Path("00-navigation"),
]
SELF_OUTPUT = Path("05-knowledge/results/lrc14_source_packet_lineage_cooccurrence_codex_s160.out")

ID_RE = re.compile(r"\b(?:HYP|THM|OPEN-Q)-[+A-Za-z0-9]*\d+\b")
HYP_NUM_RE = re.compile(r"\bHYP-\+?(\d+)\b")

CLUSTERS = {
    "qdiv_exact_period": {
        "keywords": [
            "qdiv",
            "q_threshold",
            "q-witness",
            "THM-523",
            "exact-period",
            "denominator",
            "lcm",
            "prime-power",
            "covering",
            "omit-prime",
        ],
        "fibers": {"scale", "divisibility", "exact_period"},
    },
    "ap_gw_tight_locus": {
        "keywords": [
            "AP",
            "Goddyn",
            "GW",
            "Jacobs",
            "doubling",
            "12 -> 24",
            "12->24",
            "apex",
            "Steinhaus",
            "one-hole",
        ],
        "fibers": {"boundary", "apex", "owners", "tight_atoms"},
    },
    "tournament_state_lifts": {
        "keywords": [
            "Tournament Analysis",
            "tournament",
            "isomorphism",
            "Hamiltonian",
            "score_hist",
            "SCC",
            "state lift",
            "TournamentStateLift",
            "H=7",
            "OCF",
        ],
        "fibers": {"state", "achievability", "conflict", "source"},
    },
    "farey_product_c27": {
        "keywords": [
            "Farey",
            "mediant",
            "3/41",
            "2/27",
            "C27",
            "unital",
            "K33",
            "excess",
            "branch",
            "product",
        ],
        "fibers": {"scale", "branch", "owners", "incidence"},
    },
    "haar_baire_boundary": {
        "keywords": [
            "Haar",
            "Baire",
            "regular-open",
            "strict-open",
            "safe mass",
            "endpoint",
            "boundary",
            "pinch",
            "bridge",
            "taut",
        ],
        "fibers": {"topology", "boundary", "endpoints", "measure"},
    },
    "packet_fixed_margin": {
        "keywords": [
            "packet",
            "fixed-margin",
            "family",
            "sporadic",
            "F6",
            "F7",
            "Johnson",
            "source-spectrum",
            "SourceSpec",
            "non-migration",
        ],
        "fibers": {"source", "owners", "fixed_margin", "families", "residual"},
    },
    "moon_induction_large_speed": {
        "keywords": [
            "Moon",
            "THM-571",
            "scale-separated",
            "remove-large",
            "one-large",
            "multi-large",
            "comb",
            "induction",
            "top-balanced",
            "apex-majority",
        ],
        "fibers": {"induction", "scale", "boundary_state", "descent"},
    },
    "moment_dual_certificates": {
        "keywords": [
            "dual",
            "moment",
            "Toeplitz",
            "PSD",
            "Fourier",
            "multiplicity",
            "danger-count",
            "twist",
            "Farkas",
            "eigen",
        ],
        "fibers": {"dual", "moments", "harmonics", "certificates"},
    },
    "wide_gk8_decorrelation": {
        "keywords": [
            "wide",
            "gK8",
            "L_y",
            "Delsarte",
            "Fejer",
            "Bonferroni",
            "resonance",
            "equidistribution",
            "decorrelation",
            "cap",
        ],
        "fibers": {"analytic", "moments", "wide_tail", "cap"},
    },
    "scalar_quotient_guardrails": {
        "keywords": [
            "scalar",
            "residue liar",
            "magnitude-blind",
            "too coarse",
            "refuted",
            "false",
            "not enough",
            "destroys",
            "lie",
            "dead",
        ],
        "fibers": {"guardrail"},
    },
}

TIE_PATH = [
    "packet_fixed_margin",
    "haar_baire_boundary",
    "moment_dual_certificates",
    "farey_product_c27",
    "qdiv_exact_period",
    "moon_induction_large_speed",
    "tournament_state_lifts",
    "wide_gk8_decorrelation",
    "ap_gw_tight_locus",
    "scalar_quotient_guardrails",
]

GUARDRAIL_PATTERNS = [
    "refut",
    "false",
    "too coarse",
    "not enough",
    "magnitude-blind",
    "residue liar",
    "destroy",
    "guardrail",
    "dead",
    "weak",
    "cannot",
    "no fixed",
    "lie",
    "fails",
]


def iter_docs() -> list[Path]:
    docs: list[Path] = []
    for root in ROOTS:
        if not root.exists():
            continue
        docs.extend(
            p
            for p in root.rglob("*")
            if p.suffix in {".md", ".out"} and p != SELF_OUTPUT
        )
    return sorted(docs)


def read(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def doc_id(path: Path, text: str) -> str:
    match = re.search(r"^id:\s*(\S+)", text, re.MULTILINE)
    if match:
        return match.group(1)
    match = ID_RE.search(path.name)
    if match:
        return match.group(0)
    if "poke-forum/posts" in str(path):
        return path.parent.name
    return str(path)


def status(text: str) -> str:
    match = re.search(r"^status:\s*(.+)$", text, re.MULTILINE)
    if match:
        return match.group(1).strip()
    match = re.search(r"\*\*Status:\*\*\s*(.+)$", text, re.MULTILINE)
    if match:
        return match.group(1).strip()
    return "unknown"


def hyp_number(docid: str, path: Path) -> int | None:
    match = HYP_NUM_RE.search(docid) or HYP_NUM_RE.search(path.name)
    if not match:
        return None
    return int(match.group(1))


def score_clusters(text: str) -> dict[str, int]:
    lower = text.lower()
    out: dict[str, int] = {}
    for name, spec in CLUSTERS.items():
        count = 0
        for kw in spec["keywords"]:
            count += lower.count(kw.lower())
        if count:
            out[name] = count
    return out


def is_lrc_relevant(text: str, scores: dict[str, int]) -> bool:
    lower = text.lower()
    return (
        "lrc14" in lower
        or "lrc(14" in lower
        or "open-q-108" in lower
        or sum(scores.values()) >= 5
    )


def short_line(line: str, limit: int = 190) -> str:
    clean = " ".join(line.strip().split())
    if len(clean) <= limit:
        return clean
    return clean[: limit - 3] + "..."


def directed_3cycles(vertices: list[str], orient: dict[tuple[str, str], str]) -> int:
    total = 0
    for a, b, c in combinations(vertices, 3):
        ab = orient[(a, b)] == a
        bc = orient[(b, c)] == b
        ca = orient[(c, a)] == c
        ac = orient[(a, c)] == a
        cb = orient[(c, b)] == c
        ba = orient[(b, a)] == b
        if (ab and bc and ca) or (ac and cb and ba):
            total += 1
    return total


def scc_sizes(vertices: list[str], orient: dict[tuple[str, str], str]) -> list[int]:
    graph = {v: set() for v in vertices}
    rev = {v: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        winner = orient[(a, b)]
        loser = b if winner == a else a
        graph[winner].add(loser)
        rev[loser].add(winner)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str) -> int:
        seen.add(v)
        size = 1
        for w in rev[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(vertices: list[str], orient: dict[tuple[str, str], str]) -> int:
    # Held-Karp dynamic count of directed Hamiltonian paths for n<=10.
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    pred_mask = [0] * n
    for a, b in combinations(vertices, 2):
        winner = orient[(a, b)]
        loser = b if winner == a else a
        pred_mask[idx[loser]] |= 1 << idx[winner]
    dp: dict[tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if pred_mask[nxt] & (1 << last):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def main() -> None:
    docs = []
    cluster_counts: Counter[str] = Counter()
    primary_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    cluster_docs: dict[str, list[tuple[int | None, str, str, Path]]] = defaultdict(list)
    co_edges: Counter[tuple[str, str]] = Counter()
    guardrails: list[tuple[int | None, str, str, str]] = []
    related_edges: Counter[tuple[str, str]] = Counter()

    for path in iter_docs():
        text = read(path)
        scores = score_clusters(text)
        if not is_lrc_relevant(text, scores):
            continue
        docid = doc_id(path, text)
        hnum = hyp_number(docid, path)
        stat = status(text)
        docs.append((path, text, docid, hnum, stat, scores))
        status_counts[stat.split("/")[0].strip()] += 1
        assigned = [name for name, score in scores.items() if score >= 2 or name in {"scalar_quotient_guardrails", "tournament_state_lifts"}]
        if not assigned and scores:
            assigned = [max(scores, key=scores.get)]
        for name in assigned:
            cluster_counts[name] += 1
            cluster_docs[name].append((hnum, docid, stat, path))
        if assigned:
            primary = max(assigned, key=lambda name: (scores.get(name, 0), len(CLUSTERS[name]["fibers"])))
            primary_counts[primary] += 1
        for a, b in combinations(sorted(set(assigned)), 2):
            co_edges[(a, b)] += 1
        ids = sorted(set(ID_RE.findall(text)))
        for rid in ids:
            if rid != docid and docid.startswith("HYP-") and rid.startswith(("HYP-", "THM-", "OPEN-Q-")):
                related_edges[(docid, rid)] += 1
        for line in text.splitlines():
            lower = line.lower()
            if any(pat in lower for pat in GUARDRAIL_PATTERNS):
                if "lrc" in lower or "ap" in lower or "gw" in lower or "tournament" in lower or "denominator" in lower or "scalar" in lower:
                    guardrails.append((hnum, docid, str(path), short_line(line)))

    print("# LRC14 Holistic Lineage Synthesis Mine")
    print()
    print("Scope:")
    print(f"- scanned_docs={sum(1 for _ in iter_docs())}")
    print(f"- lrc_relevant_docs={len(docs)}")
    print(f"- cluster_count={len(CLUSTERS)}")
    print("- roots=" + ", ".join(str(p) for p in ROOTS if p.exists()))
    print()

    print("## Cluster Counts")
    for name, count in cluster_counts.most_common():
        fibers = ",".join(sorted(CLUSTERS[name]["fibers"]))
        print(f"- {name}: docs={count}, primary={primary_counts[name]}, fibers={fibers}")
    print()

    print("## First/Last Hypothesis Touches By Cluster")
    for name in TIE_PATH:
        rows = cluster_docs[name]
        nums = [h for h, _, _, _ in rows if h is not None]
        if not nums:
            continue
        first = min(nums)
        last = max(nums)
        first_docs = sorted((h, d) for h, d, _, _ in rows if h == first)[:3]
        last_docs = sorted((h, d) for h, d, _, _ in rows if h == last)[:3]
        print(f"- {name}: HYP-{first} -> HYP-{last}; first={first_docs}; last={last_docs}")
    print()

    print("## Status Topline")
    for stat, count in status_counts.most_common(16):
        print(f"- {stat or 'unknown'}: {count}")
    print()

    print("## High-Signal Co-Occurrence Edges")
    for (a, b), count in co_edges.most_common(24):
        print(f"- {a} <-> {b}: {count}")
    print()

    print("## Guardrails, Refutations, And Corrections")
    # Prefer recent and unique snippets, but keep a few older corrections.
    seen_lines: set[str] = set()
    selected = []
    for hnum, docid, path, line in sorted(guardrails, key=lambda row: (row[0] is None, -(row[0] or 0), row[1])):
        key = re.sub(r"\d+", "#", line.lower())
        if key in seen_lines:
            continue
        seen_lines.add(key)
        selected.append((hnum, docid, path, line))
        if len(selected) >= 42:
            break
    for hnum, docid, path, line in selected:
        prefix = f"HYP-{hnum}" if hnum is not None else docid
        print(f"- {prefix} [{docid}] {line}")
    print()

    print("## Proof-Carrier Tournament")
    vertices = TIE_PATH[:]
    retained_score = {
        name: len(CLUSTERS[name]["fibers"]) + min(5, cluster_counts[name] // 25)
        for name in vertices
    }
    tie_rank = {name: len(vertices) - i for i, name in enumerate(TIE_PATH)}
    orient: dict[tuple[str, str], str] = {}
    scores: Counter[str] = Counter()
    for a, b in combinations(vertices, 2):
        ka = (retained_score[a], tie_rank[a])
        kb = (retained_score[b], tie_rank[b])
        winner = a if ka >= kb else b
        orient[(a, b)] = winner
        orient[(b, a)] = winner
        scores[winner] += 1
    score_hist = Counter(scores[v] for v in vertices)
    sorted_vertices = sorted(vertices, key=lambda v: (-scores[v], TIE_PATH.index(v)))
    print("Vertices are proof carriers, not runners.")
    print("Pair observable = retained source-spectrum fibers plus corpus support.")
    print("Switch/gauge = larger retained-fiber score; ties follow the Hamiltonian path.")
    print("Tie Hamiltonian path:")
    print("  " + " > ".join(TIE_PATH))
    print("Score order:")
    for v in sorted_vertices:
        print(f"- {v}: outscore={scores[v]}, retained_score={retained_score[v]}, docs={cluster_counts[v]}")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={directed_3cycles(vertices, orient)}")
    print(f"scc_sizes={scc_sizes(vertices, orient)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(vertices, orient)}")
    print("Achievable proof-carrier isomorphism class under this gauge: transitive tournament.")
    print()

    print("## Assumption Challenge")
    print("- Runner vertices alone are not selected: residue and apex tournaments repeatedly collapsed AP with magnitude liars or loose near-misses.")
    print("- Alternate vertices considered by the scanned corpus include gaps, fixed circle sections, endpoints, wall crossings, cover arcs, Fourier modes, packet families, exact-period sectors, state-lift packets, and proof obligations.")
    print("- Preserved predicate in the chosen quotient: whether the object still carries the observer-source LRC witness, exact scale, regular-open/boundary status, owner labels, and a named discharge or obstruction route.")
    print("- Destroyed information: raw row identity and many runner labels after they are converted into packet owners; this is acceptable only when exact scale and boundary/source fibers are retained.")
    print()

    print("## Synthesis Readout")
    print("- The historical flow is not many unrelated attacks.  It is a sequence of quotient tests around one proof object: SourceSpec / labelled packet / endpoint-dual carrier.")
    print("- Successful finite audits keep exact scale, boundary topology, endpoint owners, and packet labels until the last step.")
    print("- Failed or corrected routes are useful because each names a forbidden quotient: residue-only, scalar-energy-only, fixed-denominator-only, raw tournament-isomorphism-only, raw divergence, and arbitrary-digraph H=7.")
    print("- Current theorem target: every primitive 13-speed row either exits through qdiv/equidistribution/descent, is AP/GW boundary, has a positive regular-open or dual-certified safe interval, or constructs the K33/H=7 TournamentStateLift.  A counterexample must be a zero-open, qdiv>14, fixed-margin source-spectrum kernel avoiding all these exits.")


if __name__ == "__main__":
    main()
