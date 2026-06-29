#!/usr/bin/env python3
"""
Mine the repository history for niche mathematical topic clusters.

This is a lightweight archaeology tool, not a theorem prover.  It scans:

1. commit subjects across all refs;
2. current tracked text files in the research directories;
3. all historical touched paths and subjects.

The output is meant to support synthesis notes by showing which strange doors
have been repeatedly opened and where to look next.  It deliberately avoids a
full `git log -G` scan by default because this repository has enough history
that content-pickaxe passes are too slow for routine close-out sessions.
"""

from __future__ import annotations

import collections
import datetime as _dt
import pathlib
import re
import subprocess
from dataclasses import dataclass


ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "05-knowledge/results/repo_history_niche_topic_miner_codex_20260629.out"

RESEARCH_PREFIXES = (
    "00-navigation/",
    "01-canon/",
    "04-computation/",
    "05-knowledge/",
    "07-reflections/",
)

CURRENT_CONTENT_PREFIXES = (
    "00-navigation/",
    "01-canon/",
    "05-knowledge/hypotheses/",
    "05-knowledge/variables/",
    "07-reflections/",
)

SKIP_SUFFIXES = (
    ".png",
    ".jpg",
    ".jpeg",
    ".gif",
    ".pdf",
    ".sqlite",
    ".out",
    ".pyc",
)


@dataclass(frozen=True)
class Cluster:
    key: str
    regex: str
    proof_pull: str


CLUSTERS = [
    Cluster(
        "p-adic-local-lifting",
        r"\b(2-adic|p-adic|Hensel|Krasner|Kummer|Dedekind|Monna|Morita|Monsky|Hasse-Arf|local[- ]to[- ]global|non[- ]archimedean|Skolem[- ]Mahler[- ]Lech)\b",
        "Stability and lift legality: track which quotient coordinates survive small local perturbations.",
    ),
    Cluster(
        "discrepancy-haar-fejer",
        r"\b(Erd[oő]s[-– ]Tur[aá]n|Roth[-– ]Hal[aá]sz|Roth[-– ]Vaughan|Beck[-– ]Fiala|van der Corput|Haar|Fej[eé]r|Kaczynski|large sieve|circle method|discrepancy)\b",
        "Certificate backend: turn interval/fiber imbalance into signed Haar/Fejer or sieve-clock inequalities.",
    ),
    Cluster(
        "analytic-prime-sieve",
        r"\b(Barban[-– ]Davenport[-– ]Halberstam|Brun[-– ]Titchmarsh|Mertens|M[eö]bius|Mobius|Ramanujan sum|divisor function|Goldbach|ternary Goldbach|parabolic cylinder|saddle[- ]point|explicit formula|upper[- ]bound quadratic sieve)\b",
        "Smoothing dispatcher: separate admissible packet families before applying one global estimate.",
    ),
    Cluster(
        "cyclotomic-galois-q7",
        r"\b(Q\(sqrt\(-?7\)\)|Q\(sqrt\(-7\)\)|sqrt-?7|zeta_?7|ζ_?7|cyclotomic|Gauss sum|Legendre|C3|C6|Chebyshev|de Moivre|Joukowski|Hermite[-– ]Biehler|Lindemann|Weierstrass)\b",
        "Owner-charge and field sidecar: split residue, magnitude, and monodromy before quotienting.",
    ),
    Cluster(
        "automata-lacunary-sequences",
        r"\b(Moser[-– ]de Bruijn|fibbinary|Ostrowski[-– ]Hadamard|gap theorem|finite automaton|Zeckendorf|Thue[-– ]Morse|Stern[-– ]Brocot|Farey|three[- ]gap|continued fraction|Sturmian)\b",
        "Finite-state language for legal fibers, gap packets, and carry/bypass words.",
    ),
    Cluster(
        "tournament-quotient-guardrails",
        r"\b(A000568|Worpitzky|Eulerian|Lee[-– ]Yang|resolvent|Bring radical|Abel[-– ]Ruffini|canary|filler|metagraph|Hamiltonian path|score sequence|tournament)\b",
        "Controlled forgetting: every quotient must declare its retained predicate and destroyed coordinates.",
    ),
    Cluster(
        "topology-path-cuts",
        r"\b(Borel|Baire|Cech|Čech|Borsuk[-– ]Ulam|Menger|Green current|normal fan|barcode|persistent|Theta\*|ANYA|CWave|Field D\*|Schwarz[-– ]Christoffel|cut|seam|bypass)\b",
        "Relative topology: seam complements, punctures, low-rank escapes, and path-lift certificates.",
    ),
    Cluster(
        "polygonal-faulhaber-pollock",
        r"\b(Fermat polygonal|polygonal number|Pollock|Faulhaber|pentagonal|Euler partition|square pyramidal|triangular tower|odd Faulhaber|boundary moment|Bernoulli)\b",
        "Boundary-moment carrier: odd moments and additive tilings as obstruction ledgers.",
    ),
    Cluster(
        "unit-distance-codes-designs",
        r"\b(unit distance|self[- ]dual|\[72, ?36, ?16\]|Golay|unital|BIBD|block design|q=3 unital|Steiner|Paley design)\b",
        "Design view: convert local packets into incidence/block constraints with dual certificates.",
    ),
    Cluster(
        "irreducibility-polynomial-tiling",
        r"\b(Bunyakovsky|irreducib|Newton polytope|factorization|Sophie Germain|Fermat[-– ]Catalan|convolution|tiling lift|integer polynomial|Mahler measure)\b",
        "Lift obstruction: reducibility as hidden convolution lift; LRC quotients must prove no hidden lift is lost.",
    ),
    Cluster(
        "constants-thresholds-decoys",
        r"\b(Ramanujan[-– ]Soldner|Meissel[-– ]Mertens|Skewes|Graham's constant|Kempner|Champernowne|Lehmer|Lucas[-– ]Lehmer[-– ]Riesel|Brun constant)\b",
        "Canary constants: useful only after becoming exact thresholds, clocks, or named counterexamples.",
    ),
    Cluster(
        "physics-spectral-statmech",
        r"\b(octonion|G_2|Calabi[-– ]Yau|Ising|spin[- ]glass|phi\^?4|Potts|protein folding|Golay|Lee[-– ]Yang|partition function|ferromagnetic|antiferromagnetic)\b",
        "Spectral/partition lens: check whether a sidecar is a true energy certificate or only metaphor.",
    ),
]


def run_git(args: list[str], timeout: int = 45) -> str:
    return subprocess.check_output(
        ["git", *args],
        cwd=ROOT,
        text=True,
        stderr=subprocess.DEVNULL,
        timeout=timeout,
    )


def tracked_text_files() -> list[pathlib.Path]:
    files = []
    for line in run_git(["ls-files"]).splitlines():
        if not line.startswith(CURRENT_CONTENT_PREFIXES):
            continue
        if line.endswith(SKIP_SUFFIXES):
            continue
        path = ROOT / line
        if path.is_file():
            files.append(path)
    return files


def safe_read(path: pathlib.Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        try:
            return path.read_text(encoding="latin-1")
        except Exception:
            return ""
    except Exception:
        return ""


def first_matching_lines(text: str, regex: re.Pattern[str], limit: int = 4) -> list[str]:
    hits = []
    for line in text.splitlines():
        if regex.search(line):
            clean = re.sub(r"\s+", " ", line).strip()
            if clean and clean not in hits:
                hits.append(clean[:220])
            if len(hits) >= limit:
                break
    return hits


def commit_subjects() -> list[tuple[str, str, str]]:
    out = run_git(["log", "--all", "--format=%h%x09%ad%x09%s", "--date=short"])
    rows = []
    for line in out.splitlines():
        parts = line.split("\t", 2)
        if len(parts) == 3:
            rows.append((parts[0], parts[1], parts[2]))
    return rows


def history_subject_path_index() -> list[tuple[str, str, str, list[str]]]:
    out = run_git(
        ["log", "--all", "--name-only", "--format=COMMIT%x09%h%x09%ad%x09%s", "--date=short"],
        timeout=90,
    )
    rows: list[tuple[str, str, str, list[str]]] = []
    current: tuple[str, str, str] | None = None
    paths: list[str] = []
    for line in out.splitlines():
        if line.startswith("COMMIT\t"):
            if current is not None:
                rows.append((*current, paths))
            _, h, d, s = line.split("\t", 3)
            current = (h, d, s)
            paths = []
        elif line.startswith(RESEARCH_PREFIXES):
            paths.append(line)
    if current is not None:
        rows.append((*current, paths))
    return rows


def main() -> None:
    now = _dt.datetime.now().isoformat(timespec="seconds")
    files = tracked_text_files()
    subjects = commit_subjects()
    history_index = history_subject_path_index()
    all_text = {path: safe_read(path) for path in files}

    lines: list[str] = []
    lines.append("Repo History Niche Topic Miner")
    lines.append("=" * 80)
    lines.append(f"generated={now}")
    lines.append(f"tracked_text_files_scanned={len(files)}")
    lines.append(f"commit_subjects_scanned={len(subjects)}")
    lines.append(f"history_commits_with_paths_scanned={len(history_index)}")
    lines.append("")

    rows = []
    for cluster in CLUSTERS:
        regex = re.compile(cluster.regex, re.IGNORECASE)
        current_paths = []
        current_hits = 0
        samples = []
        path_counts = collections.Counter()
        for path, text in all_text.items():
            matches = regex.findall(text)
            if matches:
                current_paths.append(str(path.relative_to(ROOT)))
                current_hits += len(matches)
                path_counts[str(path.relative_to(ROOT))] = len(matches)
                if len(samples) < 8:
                    for sample in first_matching_lines(text, regex):
                        samples.append(f"{path.relative_to(ROOT)}: {sample}")
                        if len(samples) >= 8:
                            break

        subject_hits = [
            f"{h} {d} {s}" for h, d, s in subjects if regex.search(s)
        ]
        hist_rows = []
        for h, d, s, paths in history_index:
            if regex.search(s) or any(regex.search(path) for path in paths):
                hist_rows.append(f"{h} {d} {s}")
        rows.append(
            {
                "cluster": cluster,
                "current_paths": current_paths,
                "current_hits": current_hits,
                "subject_hits": subject_hits,
                "history_rows": hist_rows,
                "samples": samples,
                "path_counts": path_counts,
            }
        )

    rows.sort(
        key=lambda r: (
            len(r["history_rows"]),
            len(r["subject_hits"]),
            r["current_hits"],
            len(r["current_paths"]),
        ),
        reverse=True,
    )

    lines.append("Ranked Cluster Scorecard")
    lines.append("-" * 80)
    lines.append(
        "cluster | current_hits | current_paths | commit_subject_hits | history_subject_path_commits"
    )
    for r in rows:
        lines.append(
            f"{r['cluster'].key} | {r['current_hits']} | {len(r['current_paths'])} | "
            f"{len(r['subject_hits'])} | {len(r['history_rows'])}"
        )
    lines.append("")

    for r in rows:
        c = r["cluster"]
        lines.append(f"Cluster: {c.key}")
        lines.append("-" * 80)
        lines.append(f"regex={c.regex}")
        lines.append(f"proof_pull={c.proof_pull}")
        lines.append(f"current_hits={r['current_hits']}")
        lines.append(f"current_paths={len(r['current_paths'])}")
        lines.append(f"commit_subject_hits={len(r['subject_hits'])}")
        lines.append(f"history_subject_path_commits={len(r['history_rows'])}")
        lines.append("top_current_paths:")
        for path, count in r["path_counts"].most_common(12):
            lines.append(f"  {count:4d}  {path}")
        if r["subject_hits"]:
            lines.append("recent_subject_hits:")
            for hit in r["subject_hits"][:12]:
                lines.append(f"  {hit}")
        if r["history_rows"]:
            lines.append("history_subject_path_hits:")
            for hit in r["history_rows"][:12]:
                lines.append(f"  {hit}")
        if r["samples"]:
            lines.append("sample_current_lines:")
            for sample in r["samples"][:8]:
                lines.append(f"  {sample}")
        lines.append("")

    lines.append("Synthesis Notes")
    lines.append("-" * 80)
    lines.append(
        "The highest-value recurring pattern is not any single proper noun. It is "
        "the repeated warning that a quotient is legal only when it keeps the "
        "sidecar needed by the target predicate."
    )
    lines.append(
        "For the current random031 seam-complement frontier, the useful old topics "
        "are those that can certify one of: local lift stability, signed discrepancy, "
        "relative path reachability, finite-state gap language, design incidence, "
        "or controlled quotient forgetting."
    )
    lines.append(
        "Decoy constants should stay in the inspiration bucket until they become "
        "exact clocks, thresholds, obstruction labels, or failed-route canaries."
    )

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
