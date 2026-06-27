#!/usr/bin/env python3
"""S160: holistic lineage synthesis for the LRC14 proof workspace.

This is a metadata computation, not an LRC14 proof.  It mines the local
research records for how the proof object changed across sessions and records a
Tournament Analysis over proof carriers rather than runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]


FILES = {
    "hyp_index": ROOT / "05-knowledge/hypotheses/INDEX.md",
    "results_index": ROOT / "05-knowledge/results/INDEX.md",
    "tangents": ROOT / "00-navigation/TANGENTS.md",
    "forum_index": ROOT / "poke-forum/index.md",
    "session_log": ROOT / "00-navigation/SESSION-LOG.md",
    "concept_map": ROOT / "00-navigation/CONCEPT-MAP.md",
    "open_questions": ROOT / "00-navigation/OPEN-QUESTIONS.md",
}


KEYWORDS: dict[str, tuple[str, ...]] = {
    "observer-source/tournament": (
        "observer",
        "source",
        "tournament",
        "A000568",
        "rooted",
        "marked",
    ),
    "q-threshold/Farey scale": (
        "qdiv",
        "q_threshold",
        "q-threshold",
        "Farey",
        "denominator",
        "exact M",
    ),
    "Haar-Baire topology": (
        "Haar",
        "Baire",
        "Borel",
        "regular-open",
        "strict-open",
        "boundary-only",
    ),
    "endpoint ownership": (
        "endpoint",
        "owner",
        "pinch",
        "bridge",
        "aperture",
        "credit",
        "taut",
    ),
    "C27/unital/K33 state": (
        "C27",
        "unital",
        "K33",
        "state-lift",
        "H=7",
        "THM-572",
    ),
    "fixed-margin packet families": (
        "fixed-margin",
        "labelled-packet",
        "packet",
        "family",
        "sporadic",
        "classifier",
        "gauntlet",
    ),
    "wide/gK8 analytic tail": (
        "wide",
        "gK8",
        "Bonferroni",
        "decorrelation",
        "tail",
        "Fejer",
        "sector",
    ),
    "dual certificates": (
        "dual",
        "moment",
        "Toeplitz",
        "Fourier",
        "PSD",
        "twist",
        "polynomial",
        "barrier",
    ),
    "Ramanujan/divisor quotient": (
        "Ramanujan",
        "divisor",
        "sigma_k",
        "Dirichlet",
        "cyclotomic",
        "multiplicative",
        "quotient-admissibility",
        "exact-period",
        "primitive-unit",
        "projector",
    ),
    "AP/Goddyn-Wong equality": (
        "AP",
        "GW",
        "Goddyn",
        "Jacobsthal",
        "equality",
        "boundary atom",
    ),
    "refuted scalar quotient": (
        "refut",
        "false",
        "too weak",
        "magnitude-blind",
        "not a proof",
        "cannot",
        "guardrail",
        "negative signal",
    ),
}


PHASES = [
    ("P0 observer/tournament quotient era", 1800, 2299),
    ("P1 invariant/sheaf and finite endpoint era", 2300, 2599),
    ("P2 AP-window, sector, and far-tail era", 2600, 2829),
    ("P3 wide/gK8 and analytic residual era", 2830, 2908),
    ("P4 source-spectrum/Farey/C27/unital era", 2909, 2949),
    ("P5 packet classifier and boundary-moment era", 2950, 2969),
    ("P6 dual-certificate cluster", 2970, 2975),
    ("P7 holistic proof-object convergence", 2976, 2999),
]


ROUTE_MUTATIONS = [
    (
        "endpoint full-cover obstruction",
        "marked endpoint-owner hypergraph",
        "kept",
        "A counterexample must cover every witness interval, but endpoint data "
        "must be labelled by owner and side.",
        ("HYP-1802", "HYP-1841", "HYP-2965", "HYP-2970", "HYP-2975"),
    ),
    (
        "raw A000568 / tournament shadow",
        "observer-source and source-deleted fibers",
        "refined",
        "Unmarked and even observer-pointed tournament classes mix safe and "
        "unsafe rows; threshold/source labels are required.",
        ("HYP-1977", "HYP-2486", "HYP-2953", "HYP-2954"),
    ),
    (
        "direct finite endpoint brute force",
        "structural finite-core reductions",
        "refuted as strategy",
        "The finite endpoint space below the useful bounds is too large, and "
        "q=14V endpoint checks already fail on small rows.",
        ("HYP-2596", "HYP-2955", "HYP-2963", "HYP-2966"),
    ),
    (
        "scalar moment / additive-energy routes",
        "labelled packet and residual-leak routes",
        "refined",
        "Scalar moments see pressure but lose packet ownership; same-frequency "
        "and gK8 tails need labelled depth and residual packets.",
        ("HYP-2890", "HYP-2903", "HYP-2969", "HYP-2971", "HYP-2973"),
    ),
    (
        "Haar/Baire boundary fronts",
        "strict-open vs boundary-only packet split",
        "kept",
        "AP/GW are zero-open boundary atoms; most near rows are positive-open "
        "and need exact interval/owner ledgers rather than analogy.",
        ("HYP-2948", "HYP-2949", "HYP-2951", "HYP-2955"),
    ),
    (
        "C27/unital/K33 analogies",
        "state-lift owner labels",
        "refined",
        "Design and graph carriers are useful only as typed labels attached to "
        "exact M/Farey/Haar packets.",
        ("HYP-2891", "HYP-2892", "HYP-2937", "HYP-2954", "HYP-2908"),
    ),
    (
        "fixed-margin family idea",
        "labelled packet theorem",
        "promoted",
        "Families are fixed-margin labelled packet classes; sporadics are "
        "bounded singleton packets after family coordinates are removed.",
        ("HYP-2956", "HYP-2961", "HYP-2962", "HYP-2963"),
    ),
    (
        "boundary-moment bridge",
        "multi-chart feasible-region packet",
        "refined",
        "One all-covered denominator chart is not an obstruction; covering rows "
        "can be positive Haar-open.",
        ("HYP-2965", "HYP-2968", "HYP-2969", "HYP-2975"),
    ),
    (
        "dual certificate scouts",
        "projections of one cover/noncover sheaf",
        "merged",
        "Endpoint potentials, multiplicity barriers, twists, danger counts, and "
        "Toeplitz negativity are complementary Farkas-like shadows.",
        ("HYP-2970", "HYP-2971", "HYP-2972", "HYP-2973", "HYP-2974"),
    ),
    (
        "scalar divisor signature",
        "Ramanujan/cyclotomic packet quotient",
        "reserved as guardrail",
        "Divisor data are admissible only when the retained labels distinguish "
        "AP/GW boundary atoms, Toeplitz/Ramanujan exits, and K33/state-lift debts.",
        ("HYP-2979", "HYP-2978", "HYP-2977", "HYP-2976", "THM-406", "THM-572"),
    ),
]


COUNTEREXAMPLE_SIEVE = [
    ("primitive normalization", "work with 13 nonzero speeds at total n=14"),
    ("q-clock witness", "qdiv < 14 gives a direct strict phase"),
    ("AP/GW boundary atoms", "qdiv = 14 and zero-open routes to the known equality atoms"),
    ("positive Haar/Baire open front", "any strict interval front discharges the row"),
    ("unit petal / two-block discharge", "C27 unit-petal packets are positive-open"),
    ("K33 / state-lift obligation", "nonunit K33 packets are not scalar rows; route to THM-572"),
    ("THM-571 apex-majority discharge", "|14Z cap S| >= 7 leaves the live Moon core"),
    ("few-apex lift packet", "1 <= |14Z cap S| <= 6 must have positive lift mass or state-lift"),
    ("NORK pinch atlas", "non-AP/GW zero-open residuals have not appeared in the large local bank"),
    ("Ramanujan exact-period packet labels", "scalar divisor signatures must refine to primitive cyclotomic/Fourier endpoint labels"),
    ("boundary-moment ledger", "covering packets need a labelled multi-chart moment bridge"),
    ("dual certificate wall", "moment, twist, endpoint, and Toeplitz duals certify all positive audited hard rows"),
]


@dataclass(frozen=True)
class HypEntry:
    hyp: str
    num: int
    source: str
    title: str
    status: str
    body: str
    refs: tuple[str, ...]


@dataclass(frozen=True)
class Route:
    name: str
    vector: tuple[int, ...]
    preserves: str
    destroys: str


def read(name: str) -> str:
    return FILES[name].read_text(errors="replace")


def hyp_num(hyp: str) -> int:
    m = re.search(r"HYP-\+?(\d+)", hyp)
    return int(m.group(1)) if m else -1


def parse_hyp_index(text: str) -> list[HypEntry]:
    entries: list[HypEntry] = []
    for line in text.splitlines():
        if not line.startswith("- **HYP-"):
            continue
        m = re.match(
            r"- \*\*(HYP-\+?\d+) / ([^(]+) \(([^)]+)\):\*\* ([^.]+)\.  (.*)",
            line,
        )
        if not m:
            continue
        hyp, source, title, status, body = m.groups()
        refs = tuple(dict.fromkeys(re.findall(r"HYP-\+?\d+|THM-\d+|OPEN-Q-\d+", body)))
        entries.append(
            HypEntry(
                hyp=hyp,
                num=hyp_num(hyp),
                source=source.strip(),
                title=title.strip(),
                status=status.strip(),
                body=body.strip(),
                refs=refs,
            )
        )
    return entries


def parse_hyp_files() -> list[HypEntry]:
    entries: list[HypEntry] = []
    for path in sorted((ROOT / "05-knowledge/hypotheses").glob("HYP-*.md")):
        text = path.read_text(errors="replace")
        m = re.match(r"(HYP-\+?\d+)", path.name)
        if not m:
            continue
        hyp = m.group(1)
        title = path.stem
        heading = re.search(r"^#\s+(.*)$", text, re.M)
        if heading:
            title = heading.group(1).replace(hyp + ":", "").strip()
        status = "detail file"
        sm = re.search(r"(?im)^\*\*Status:\*\*\s*(.*)$|^status:\s*(.*)$", text)
        if sm:
            status = next(g for g in sm.groups() if g).strip()
        refs = tuple(dict.fromkeys(re.findall(r"HYP-\+?\d+|THM-\d+|OPEN-Q-\d+", text)))
        entries.append(
            HypEntry(
                hyp=hyp,
                num=hyp_num(hyp),
                source="detail-file",
                title=title,
                status=status,
                body=text.replace("\n", " "),
                refs=refs,
            )
        )
    return entries


def phase_for(num: int) -> str:
    for name, lo, hi in PHASES:
        if lo <= num <= hi:
            return name
    return "other"


def keyword_hits(text: str) -> Counter[str]:
    lower = text.lower()
    counts: Counter[str] = Counter()
    for label, words in KEYWORDS.items():
        for word in words:
            if word.lower() in lower:
                counts[label] += 1
    return counts


def lrc_records(entries: list[HypEntry]) -> list[HypEntry]:
    out = []
    for entry in entries:
        blob = f"{entry.title} {entry.status} {entry.body}".lower()
        if (
            "lrc14" in blob
            or "lrc(14)" in blob
            or "n=14" in blob
            or "lonely runner" in blob
            or re.search(r"\blrc\b", blob)
        ):
            out.append(entry)
    return out


def top_refs(entries: list[HypEntry]) -> Counter[str]:
    c: Counter[str] = Counter()
    for entry in entries:
        c.update(entry.refs)
    return c


def hamiltonian_path_count(routes: list[Route]) -> tuple[int, list[str]]:
    def beats(a: Route, b: Route) -> bool:
        if a.vector != b.vector:
            return a.vector > b.vector
        return a.name < b.name

    count = 0
    first_path: list[str] = []
    for perm in permutations(routes):
        if all(beats(perm[i], perm[i + 1]) for i in range(len(perm) - 1)):
            count += 1
            if not first_path:
                first_path = [r.name for r in perm]
    return count, first_path


def route_tournament() -> tuple[dict[int, int], int, int, list[str]]:
    routes = [
        Route(
            "proof-object sheaf",
            (9, 9, 9, 9, 9, 9, 8, 9, 9),
            "strict predicate, exact scale, topology, owners, packets, state-lift, duals",
            "nothing intentionally; this is the synthesis object",
        ),
        Route(
            "dual-certificate cluster",
            (9, 8, 8, 5, 5, 5, 9, 8, 8),
            "cover/noncover certificates through endpoint, moment, twist, and Fourier duals",
            "some packet ownership unless reattached",
        ),
        Route(
            "Ramanujan/cyclotomic packet quotient",
            (9, 8, 8, 6, 8, 7, 6, 8, 9),
            "divisor/cyclotomic phase labels with explicit quotient-admissibility guards",
            "endpoint geometry unless reattached to Farey/Haar owners",
        ),
        Route(
            "labelled packet classifier",
            (9, 8, 8, 8, 9, 9, 5, 9, 9),
            "families, sporadics, endpoint owners, C27/K33 labels",
            "some analytic harmonic detail",
        ),
        Route(
            "boundary endpoint bridge",
            (9, 8, 9, 9, 7, 7, 5, 8, 8),
            "open-vs-boundary status, safe bridges, endpoint owners",
            "global fixed-margin and Fourier detail",
        ),
        Route(
            "source-spectrum observer pullback",
            (9, 8, 7, 6, 7, 7, 4, 7, 8),
            "observer-source predicate with Farey/time movie",
            "fine endpoint ownership unless labelled",
        ),
        Route(
            "wide/gK8 analytic route",
            (8, 6, 5, 4, 5, 5, 5, 8, 7),
            "positive mass/cap bounds for wide and decorrelated branches",
            "low-frontier packet identities",
        ),
        Route(
            "AP/GW boundary skeleton",
            (7, 7, 9, 8, 7, 6, 3, 7, 7),
            "tight equality atoms and denominator-14 boundary support",
            "strict positive rows",
        ),
        Route(
            "raw tournament shadow",
            (4, 2, 2, 1, 1, 1, 1, 3, 2),
            "some pairwise/cyclic pressure",
            "magnitude, topology, endpoints, labels",
        ),
        Route(
            "raw scalar invariant",
            (3, 3, 1, 1, 1, 1, 1, 2, 1),
            "one-dimensional stress signals",
            "most proof fibers",
        ),
    ]
    scores: Counter[str] = Counter()
    c3 = 0
    for a, b in combinations(routes, 2):
        winner = a if (a.vector, -ord(a.name[0])) > (b.vector, -ord(b.name[0])) else b
        scores[winner.name] += 1
        scores.setdefault((b if winner is a else a).name, 0)
    by_score = Counter(scores.values())

    def beats(a: Route, b: Route) -> bool:
        if a.vector != b.vector:
            return a.vector > b.vector
        return a.name < b.name

    for a, b, c in combinations(routes, 3):
        if beats(a, b) and beats(b, c) and beats(c, a):
            c3 += 1
        if beats(a, c) and beats(c, b) and beats(b, a):
            c3 += 1

    hp_count, hp = hamiltonian_path_count(routes)
    return dict(sorted(by_score.items())), c3, hp_count, hp


def format_refs(counter: Counter[str], n: int = 12) -> str:
    return ", ".join(f"{k}:{v}" for k, v in counter.most_common(n))


def main() -> None:
    hyp_text = read("hyp_index")
    index_entries = parse_hyp_index(hyp_text)
    file_entries = parse_hyp_files()
    by_hyp: dict[str, HypEntry] = {entry.hyp: entry for entry in file_entries}
    by_hyp.update({entry.hyp: entry for entry in index_entries})
    entries = sorted(by_hyp.values(), key=lambda e: (e.num, e.hyp, e.title))
    lrc_entries = lrc_records(entries)
    all_text = "\n".join(read(name) for name in FILES)

    print("S160 LRC14 HOLISTIC LINEAGE SYNTHESIS")
    print("=" * 78)
    print("[0] Corpus")
    print(f"  parsed hypothesis-index rows: {len(index_entries)}")
    print(f"  parsed hypothesis detail files: {len(file_entries)}")
    print(f"  merged unique hypothesis rows: {len(entries)}")
    print(f"  LRC/LRC14/N=14 hypothesis rows: {len(lrc_entries)}")
    for name, path in FILES.items():
        text = path.read_text(errors="replace")
        print(
            f"  {name:14s}: lines={text.count(chr(10))+1:5d} "
            f"lrc14={text.lower().count('lrc14'):4d} hyp={text.count('HYP-'):5d}"
        )
    print()

    print("[1] Phase ledger")
    phase_entries: dict[str, list[HypEntry]] = defaultdict(list)
    for entry in lrc_entries:
        phase_entries[phase_for(entry.num)].append(entry)
    for phase, lo, hi in PHASES:
        rows = sorted(phase_entries.get(phase, []), key=lambda e: e.num)
        if not rows:
            continue
        print(f"  {phase}: {len(rows)} rows")
        reps = rows[:2] + rows[-3:]
        seen = set()
        for row in reps:
            if row.hyp in seen:
                continue
            seen.add(row.hyp)
            print(f"    {row.hyp:9s} {row.title[:66]}")
    print()

    print("[2] Route-family keyword load")
    family_counts: Counter[str] = Counter()
    phase_family: dict[str, Counter[str]] = defaultdict(Counter)
    for entry in lrc_entries:
        hits = keyword_hits(f"{entry.title} {entry.status} {entry.body}")
        family_counts.update(hits)
        phase_family[phase_for(entry.num)].update(hits)
    for label, count in family_counts.most_common():
        print(f"  {label:30s} {count:4d}")
    print()
    print("  phase-leading families:")
    for phase, _, _ in PHASES:
        if phase_family.get(phase):
            print(f"    {phase}: {format_refs(phase_family[phase], 6)}")
    print()

    print("[3] Refuted/refined quotient signals")
    refute_lines = []
    for entry in lrc_entries:
        blob = f"{entry.status}. {entry.body}"
        if keyword_hits(blob)["refuted scalar quotient"]:
            refute_lines.append(entry)
    print(f"  rows with guardrail/refutation language: {len(refute_lines)}")
    for entry in refute_lines[:18]:
        fragment = entry.body[:190].replace("`", "")
        print(f"    {entry.hyp:9s} {entry.title[:44]:44s} :: {fragment}...")
    print()

    print("[4] Dependency hubs in LRC14 rows")
    refs = top_refs(lrc_entries)
    print(f"  top refs: {format_refs(refs, 16)}")
    print()

    print("[5] Current convergence cluster")
    cluster_ids = [
        "HYP-2953",
        "HYP-2954",
        "HYP-2956",
        "HYP-2961",
        "HYP-2962",
        "HYP-2963",
        "HYP-2964",
        "HYP-2965",
        "HYP-2966",
        "HYP-2967",
        "HYP-2968",
        "HYP-2969",
        "HYP-2970",
        "HYP-2971",
        "HYP-2972",
        "HYP-2973",
        "HYP-2974",
        "HYP-2975",
        "HYP-2976",
        "HYP-2977",
        "HYP-2978",
        "HYP-2979",
    ]
    by_id = {e.hyp: e for e in entries}
    for hid in cluster_ids:
        if hid in by_id:
            e = by_id[hid]
            hits = ", ".join(keyword_hits(f"{e.title} {e.body}").keys())
            print(f"  {hid:9s} {e.title[:54]:54s} | {hits}")
        elif hid == "HYP-2975":
            print("  HYP-2975 LRC14 taut bridge graph curvature              | endpoint ownership, dual certificates")
        elif hid == "HYP-2976":
            print("  HYP-2976 current holistic synthesis                     | proof-object sheaf")
    print()

    print("[6] Route mutations: tested, refuted, refined")
    for source, target, verdict, lesson, support in ROUTE_MUTATIONS:
        suffix = f" support={','.join(support)}" if support else ""
        print(f"  {verdict.upper():18s} {source} -> {target}{suffix}")
        print(f"    {lesson}")
    print()

    print("[7] Counterexample-family squeeze")
    for i, (gate, readout) in enumerate(COUNTEREXAMPLE_SIEVE, 1):
        print(f"  {i:02d}. {gate:34s} :: {readout}")
    print()

    print("[8] Tournament Analysis over proof carriers")
    print("  vertices considered:")
    print("    runners, raw arcs, endpoints, safe intervals, fixed-margin packets,")
    print("    C27/unital/K33 labels, state-lift obligations, exact-period ledgers,")
    print("    Ramanujan/divisor packets, multiplicity distributions, twist ladders,")
    print("    Fourier modes, proof routes.")
    print("  chosen vertices:")
    print("    proof carriers / quotient maps, because this session asks which")
    print("    structures preserve the LRC predicate through the history.")
    print("  pair observable:")
    print("    retention vector = (strict predicate, scale, topology, endpoints,")
    print("    packet labels, state-lift, dual certificate, finite atlas,")
    print("    anti-scalar guardrail).")
    print("  switch/gauge:")
    print("    lexicographic retention; ties use carrier name.")
    hist, c3, hp_count, hp = route_tournament()
    print(f"  fingerprint: score_hist={hist} c3={c3} hp={hp_count}")
    print(f"  Hamiltonian path: {' > '.join(hp)}")
    print()

    print("[9] Synthesis theorem target")
    print("  Finite proof-object sheaf theorem for LRC14:")
    print("    Every primitive 13-speed row reaches one of the named sheaf exits:")
    print("    q-witness, AP/GW boundary equality, unit-petal discharge,")
    print("    K33/state-lift obligation, positive boundary/lift packet,")
    print("    fixed-margin labelled family discharge, Ramanujan/cyclotomic")
    print("    exact-period quotient-admissibility packet, or a dual certificate from endpoint")
    print("    potential, multiplicity moment, danger-count moment, rational twist,")
    print("    or Fourier-Toeplitz negativity.")
    print("  A proposed counterexample must be invisible to every current quotient;")
    print("  the remaining work is to prove that invisibility forces AP/GW or a")
    print("  named state-lift packet, not a new scalar extremizer.")


if __name__ == "__main__":
    main()
