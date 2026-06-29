#!/usr/bin/env python3
"""HYP-3512: random031 niche-history connection scanner.

This is a repo-history archaeology scout, not a proof.  It searches current
archive files plus full commit subject/path history for older niche carriers,
then scores them against the current random031 seam-complement receiver.

Tournament Analysis declaration:
  vertices: niche proof carriers, not runners, gates, rows, residues, or raw
            git commits;
  pairwise observable: preserved random031 receiver payload across nine axes;
  switch/gauge: higher weighted receiver score wins, with evidence count only
                as a weak tie-breaker;
  tie Hamiltonian path: topology/quotient carriers before scalar analogies.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
import math
import re
import subprocess
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
SEARCH_ROOTS = (
    "00-navigation",
    "05-knowledge/hypotheses",
    "05-knowledge/results",
    "07-reflections",
    "04-computation",
)

AXES = (
    "seam_complement_topology",
    "bypass_phase_flow",
    "private_label_firewall",
    "owner_monodromy_charge",
    "quotient_guardrail",
    "two_adic_recursion",
    "rank2_escape_dual",
    "formal_receiver_readiness",
    "false_positive_control",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    aliases: tuple[str, ...]
    anchors: tuple[str, ...]
    scores: tuple[int, ...]
    connection: str
    guardrail: str
    next_experiment: str

    @property
    def base_score(self) -> int:
        return sum(self.scores)


def carrier(
    name: str,
    aliases: Iterable[str],
    anchors: Iterable[str],
    scores: tuple[int, ...],
    connection: str,
    guardrail: str,
    next_experiment: str,
) -> Carrier:
    if len(scores) != len(AXES):
        raise ValueError(f"{name} has {len(scores)} scores for {len(AXES)} axes")
    return Carrier(
        name=name,
        aliases=tuple(aliases),
        anchors=tuple(anchors),
        scores=scores,
        connection=connection,
        guardrail=guardrail,
        next_experiment=next_experiment,
    )


CARRIERS = (
    carrier(
        "coordinate_resurrection_sheaf",
        ("coordinate resurrection", "coordinate_resurrection", "repair-cover", "repair cover"),
        ("HYP-3118", "HYP-3120", "LRC-LENS-MAP"),
        (5, 4, 4, 4, 5, 4, 4, 5, 5),
        "The vertical half-turn and any seam-complement quotient must name the first destroyed coordinate: free-hole status, branch sheet, bypass purity, or endpoint rank.",
        "Do not use coordinate resurrection as a slogan; require cover rank, live section, and terminal exit/debt.",
        "Add a random031 quotient-audit table with destroyed fields and minimal repair sidecars for mirror, vertical half-turn, owner-word, and u-fiber quotients.",
    ),
    carrier(
        "normal_fan_cech_barcode",
        ("normal-fan", "normal fan", "cech", "barcode", "Cech"),
        ("HYP-3034", "HYP-3096", "HYP-3101", "LTI-345"),
        (5, 4, 3, 4, 5, 3, 4, 4, 5),
        "The four dead islands and 14 free-hole packets want a relative Cech/barcode test on the seam complement, not another projection-edge search.",
        "Topology alone forgets owner labels and branch legality unless barcode components retain owner and mirror words.",
        "Build a relative Cech/H1 complex whose boundary is rank-2 escape plus forbidden seam, then test owner-deletion persistence on bypass/free-hole components.",
    ),
    carrier(
        "finite_address_phi_tuple",
        ("finite-address", "finite_address", "Phi tuple", "phi tuple", "finite address"),
        ("HYP-3083", "HYP-3091", "HYP-3119", "HYP-3120"),
        (4, 4, 4, 4, 5, 4, 4, 5, 4),
        "Random031 needs a receiver tuple: seam status, u-index, branch, mirror mate, hit rank, owner word, firewall status, and terminal exit.",
        "A finite address without the endpoint/owner payload becomes another row hash.",
        "Define `Random031SeamComplementPacket` fields and check that HYP-3486 cells instantiate them without raw witness-count shortcuts.",
    ),
    carrier(
        "fiber_pgf_sheet_count",
        ("fiber-PGF", "fiber PGF", "sheet-count", "sheet count", "Q-masked"),
        ("HYP-3140", "HYP-3137", "HYP-3136", "HYP-3485"),
        (4, 5, 2, 3, 5, 5, 3, 4, 5),
        "HYP-3486's 258 occupied fibers and 24 double fibers are exactly the kind of sheet-count packet HYP-3140 warned not to scalarize.",
        "Raw `282` or raw `12` hides which fibers are bypass, free-hole, or ordinary and whether vertical gluing is legal.",
        "Compute a class-conditioned fiber PGF for random031: bypass/free-hole/ordinary by branch, mirror orbit, endpoint rank, and owner word.",
    ),
    carrier(
        "observer_cut_payload_orbit",
        ("observer-cut", "observer cut", "observer payload", "observer_gluing", "source-perspective"),
        ("HYP-3054", "HYP-3056", "HYP-3102", "HYP-3119"),
        (4, 4, 4, 5, 5, 3, 4, 5, 5),
        "The seven-owner seam is an observer-cut payload orbit: deleting the seam is legal only if the complement remembers the owner charge crossing it.",
        "A visible phase route without the observer payload cannot account for seam-only owners `(45,147,169,173)`.",
        "Attach observer-cut payload orbits to the pure 12-cell bypass component and compare bypass owners `(23,93,113)` with seam-only owners.",
    ),
    carrier(
        "q27_q31_resource_descent",
        ("Q27", "Q31", "q27", "q31", "C=27", "C27", "resource descent"),
        ("HYP-2470", "HYP-2471", "HYP-2480", "HYP-3483", "HYP-3484"),
        (4, 4, 3, 5, 4, 5, 3, 4, 4),
        "`C=2n-1=27` already sits at the seam: additive owner debt is n+2-like while phase flow is n*2-like.",
        "Resource descent is not a proof unless it names the exact owner/bypass fields it normalizes.",
        "Compare the seven-owner seam word to the C=27 carry/doubling split and ask which owners are stationary boundary debt versus moving phase carriers.",
    ),
    carrier(
        "endpoint_credit_farkas",
        ("endpoint-credit", "endpoint credit", "Farkas", "dual certificate", "Menger", "Green"),
        ("HYP-3437", "HYP-3451", "HYP-3417", "HYP-3120"),
        (3, 4, 4, 4, 4, 3, 5, 4, 5),
        "Rank-2 seam-complement routing should become a dual certificate: every routed cell has one rank-2 endpoint gate, and free holes are named exceptions.",
        "A max-flow/min-cut shadow is illegal if it forgets private labels or owner monodromy.",
        "Turn the HYP-3486 rank-2 exits into an endpoint-credit/Farkas ledger over components, with separate free-hole and bypass rows.",
    ),
    carrier(
        "et_hensel_fiber_zipper",
        ("ET/Hensel", "Hensel", "fiber zipper", "coarse-to-exact", "ET_unit"),
        ("HYP-3020", "HYP-3024", "HYP-3119"),
        (3, 4, 3, 3, 5, 4, 3, 4, 4),
        "The two-adic cylinder coordinate is a coarse-to-exact zipper: u-fibers split by branch and legality, not by raw residue alone.",
        "Hensel language is only useful if it distinguishes legal mirror from illegal vertical half-turn gluing.",
        "Test whether u-fiber classes form a zipper: singleton/double fibers, branch sheet, hit rank, free-hole status, and mirror mate.",
    ),
    carrier(
        "proof_circuit_missing_input",
        ("proof circuit", "missing-input", "missing input", "circuit certificate", "proof-circuit"),
        ("HYP-3116", "HYP-3117", "HYP-3120", "LTI-512"),
        (3, 3, 4, 4, 4, 3, 4, 5, 5),
        "Every tempting random031 shortcut has a missing-input vector: branch, u, mirror, owner, firewall, rank, and free-hole status.",
        "A small Boolean closure over raw flags will overfit unless all essential sidecars are marked.",
        "Run a missing-input audit for candidate random031 lemmas: rank-2 route, free-hole packet, pure bypass, and private firewall.",
    ),
    carrier(
        "chiral_stalk_cech_guard",
        ("chiral-stalk", "chiral stalk", "mirror/converse", "orientation payload", "state_lift_sign"),
        ("HYP-3123", "LRC-LENS-MAP", "HYP-3486"),
        (4, 4, 3, 3, 5, 4, 3, 4, 4),
        "The mirror-paired bypass blocks are oriented: mirror is legal and class-preserving, vertical half-turn is not.",
        "A mirror quotient that ignores chirality can mix free-hole and ordinary cells.",
        "Add a chiral guard word to each random031 component: mirror mate, branch orientation, vertical-halfturn legality, and class preservation.",
    ),
    carrier(
        "exact_period_totient_packet",
        ("exact-period", "exact period", "Mobius", "totient", "CRT curvature"),
        ("HYP-2886", "HYP-2899", "HYP-2900", "HYP-3311"),
        (3, 3, 3, 3, 4, 5, 3, 3, 4),
        "The q=14V phase pullback and q=2422 fiber graph are exact-period ledgers; period data should be packetized before analytic shortcuts.",
        "Totient/CRT counts are false friends if they forget branch filters and owner charges.",
        "Record the exact-period packet for the 12 bypass cells and 40 free-hole cells: residue class, u-fiber, branch, mirror orbit, and hit rank.",
    ),
    carrier(
        "resolvent_middle_layer_packet",
        ("resolvent", "middle-layer", "SPEC", "De Moivre", "bounded-core"),
        ("HYP-3135", "HYP-3139", "HYP-3486"),
        (3, 4, 2, 4, 4, 3, 3, 3, 4),
        "The bypass component may be a finite middle-layer payload: SPEC/owner-current debt should be named only after bypass purity is retained.",
        "Raw root or SPEC constants do not know the punctures or private-label firewall.",
        "For the pure 12-cell bypass, compute the owner-current/SPEC sidecar before declaring seam-only owners as residual debt.",
    ),
    carrier(
        "signed_polymer_dirichlet_network",
        ("signed polymer", "Dirichlet", "F7", "state-lift", "Ising domain"),
        ("HYP-3120", "THM-572", "HYP-3115", "HYP-2908"),
        (2, 3, 3, 3, 3, 3, 3, 3, 4),
        "If the seam-complement packet emits named debt, F7/state-lift should distinguish real terminal obstruction from phantom boundary atom.",
        "Energy/conductance analogies are late-stage tests, not carriers for the random031 topology itself.",
        "Use only after rank-2/free-hole/bypass lemmas fail: test whether the remaining owner debt forms an F7/state-lift exit.",
    ),
    carrier(
        "vitali_antipoisson_width",
        ("Vitali", "anti-Poisson", "anti Poisson", "coimage"),
        ("THM-406", "HYP-2152", "HYP-2936"),
        (3, 2, 2, 3, 4, 2, 2, 2, 4),
        "Warns that low-order moments can miss all-order cancellation; random031's scalar counts are shadows of richer packets.",
        "Does not directly see branch legality, private labels, or owner seam charges.",
        "Use as a guardrail against replacing the fiber graph by raw counts `282`, `242`, `40`, or `12`.",
    ),
    carrier(
        "ocf_forbidden_h7_packet",
        ("OCF", "forbidden H7", "H=7", "hard-core partition"),
        ("HYP-2908", "THM-572", "HYP-2936"),
        (2, 2, 2, 3, 3, 2, 2, 3, 4),
        "A terminal fallback if owner-current/two-adic/SPEC debt becomes a real forbidden state-lift packet.",
        "Raw H-values are too late and too coarse for the seam-complement geometry.",
        "Keep as named residual boundary only after packet closure fails.",
    ),
    carrier(
        "zeckendorf_partial_cube_shadow",
        ("Zeckendorf", "Fibonacci-cube", "partial-cube", "fibbinary", "Theta class"),
        ("HYP-3146", "LRC-LENS-MAP", "SESSION-LOG"),
        (2, 2, 2, 2, 3, 3, 2, 2, 3),
        "Partial-cube cut language is a possible telemetry shadow for forbidden adjacency and mirror packet cuts.",
        "Sequence automata do not know random031 owner labels unless lifted into a cut carrier.",
        "Use only as a visualization/control vocabulary for bypass/free-hole adjacency masks.",
    ),
    carrier(
        "sexy_prime_residue_sieve_echo",
        ("sexy-prime", "sexy prime", "prime residue", "residue-sieve"),
        ("HYP-3120", "SESSION-LOG"),
        (1, 1, 1, 1, 2, 1, 1, 1, 5),
        "A useful negative control: local residue sieves can look structured while lacking the analytic sidecar needed for closure.",
        "No direct random031 proof content without analytic distribution input.",
        "Use as a false-friend detector only.",
    ),
)


def decode_bytes(raw: bytes) -> str:
    if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        try:
            return raw.decode("utf-16")
        except UnicodeDecodeError:
            pass
    return raw.decode("utf-8", errors="ignore")


def git(*args: str, timeout: int = 30) -> str:
    proc = subprocess.run(
        ("git", *args),
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=True,
    )
    return proc.stdout


def tracked_files() -> list[Path]:
    out = git("ls-files", *SEARCH_ROOTS)
    files: list[Path] = []
    for line in out.splitlines():
        path = ROOT / line
        if path.suffix.lower() in {".md", ".txt", ".out", ".py", ".lean"}:
            files.append(path)
    return files


def text_suffix(path: str) -> bool:
    return Path(path).suffix.lower() in {".md", ".txt", ".out", ".py", ".lean"}


def alias_hits(text: str, aliases: tuple[str, ...]) -> int:
    lower = text.lower()
    return sum(lower.count(alias.lower()) for alias in aliases)


def line_samples(path: Path, text: str, aliases: tuple[str, ...], limit: int = 3) -> list[str]:
    samples: list[str] = []
    lowered_aliases = tuple(alias.lower() for alias in aliases)
    for number, line in enumerate(text.splitlines(), 1):
        low = line.lower()
        if any(alias in low for alias in lowered_aliases):
            rel = path.relative_to(ROOT)
            snippet = re.sub(r"\s+", " ", line.strip())
            samples.append(f"{rel}:{number}: {snippet[:180]}")
            if len(samples) >= limit:
                break
    return samples


def current_archive_evidence(files: list[Path]) -> dict[str, dict[str, object]]:
    evidence: dict[str, dict[str, object]] = {}
    texts: list[tuple[Path, str]] = []
    for path in files:
        try:
            text = decode_bytes(path.read_bytes())
        except OSError:
            continue
        texts.append((path, text))

    for item in CARRIERS:
        hits = 0
        by_file: Counter[str] = Counter()
        samples: list[str] = []
        for path, text in texts:
            count = alias_hits(text, item.aliases)
            if count:
                hits += count
                by_file[str(path.relative_to(ROOT))] += count
                if len(samples) < 5:
                    samples.extend(line_samples(path, text, item.aliases, limit=5 - len(samples)))
        evidence[item.name] = {
            "hits": hits,
            "files": len(by_file),
            "top_files": by_file.most_common(5),
            "samples": samples,
        }
    return evidence


def git_history_evidence() -> dict[str, dict[str, object]]:
    log = git(
        "log",
        "--all",
        "--format=commit%x09%h%x09%cs%x09%s",
        "--name-only",
        *("--", *SEARCH_ROOTS),
        timeout=60,
    )
    subject_lines: list[str] = []
    path_lines: list[str] = []
    commits_seen: set[str] = set()
    for line in log.splitlines():
        if line.startswith("commit\t"):
            parts = line.split("\t", 3)
            if len(parts) == 4:
                commits_seen.add(parts[1])
                subject_lines.append("\t".join(parts[1:]))
            continue
        path = line.strip()
        if path and text_suffix(path):
            path_lines.append(path)

    evidence: dict[str, dict[str, object]] = {}
    for item in CARRIERS:
        needles = item.aliases + item.anchors
        subject_hits = [
            line
            for line in subject_lines
            if any(alias.lower() in line.lower() for alias in needles)
        ][:8]
        path_matches = [
            path
            for path in path_lines
            if any(alias.lower() in path.lower() for alias in needles)
        ]
        unique_paths = sorted(set(path_matches))
        evidence[item.name] = {
            "subject_hits": subject_hits,
            "path_hits": len(path_matches),
            "path_files": len(unique_paths),
            "path_samples": unique_paths[:8],
            "commits_scanned": len(commits_seen),
            "path_events_scanned": len(path_lines),
            "unique_paths_scanned": len(set(path_lines)),
        }
    return evidence


def tournament_rank(
    archive: dict[str, dict[str, object]],
    history: dict[str, dict[str, object]],
) -> list[tuple[Carrier, int, int, int, float]]:
    rows = []
    for item in CARRIERS:
        current_hits = int(archive[item.name]["hits"])
        subject_hits = len(history[item.name]["subject_hits"])
        path_hits = int(history[item.name]["path_hits"])
        evidence_bonus = (
            min(7, int(math.log2(current_hits + 1)))
            + min(3, subject_hits)
            + min(2, int(math.log2(path_hits + 1)))
        )
        total = item.base_score + evidence_bonus
        rows.append((item, item.base_score, evidence_bonus, total, current_hits))
    return sorted(rows, key=lambda row: (-row[3], -row[1], row[0].name))


def tournament_fingerprint(ranked: list[tuple[Carrier, int, int, int, float]]) -> dict[str, object]:
    totals = {item.name: total for item, _, _, total, _ in ranked}
    score_hist = Counter(totals.values())
    names = [item.name for item, *_ in ranked]
    edge_flips_against_tie_path = 0
    directed_3cycles = 0
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            if totals[a] < totals[b]:
                edge_flips_against_tie_path += 1
    # Higher total gives a transitive tournament except exact-score ties; the
    # fixed ranked order orients ties, so there are no directed cycles.
    return {
        "vertices": len(names),
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": directed_3cycles,
        "scc_sizes": [1 for _ in names],
        "hamiltonian_path_count_under_fixed_tie_path": 1,
        "edge_flips_against_tie_path": edge_flips_against_tie_path,
        "hamiltonian_path": names,
    }


def print_section(title: str) -> None:
    print(f"\n## {title}")


def main() -> None:
    files = tracked_files()
    archive = current_archive_evidence(files)
    history = git_history_evidence()
    ranked = tournament_rank(archive, history)
    fp = tournament_fingerprint(ranked)

    print("HYP-3512 RANDOM031 NICHE-HISTORY CONNECTION SCANNER")
    print("status=EVIDENCE / repo-history synthesis scout; not an LRC14 proof")
    print(f"tracked_files_scanned={len(files)}")
    any_history = next(iter(history.values()))
    print(f"history_commits_scanned={any_history['commits_scanned']}")
    print(f"history_path_events_scanned={any_history['path_events_scanned']}")
    print(f"history_unique_paths_scanned={any_history['unique_paths_scanned']}")
    print(f"carrier_count={len(CARRIERS)}")
    print(f"receiver_axes={','.join(AXES)}")

    print_section("Top carrier bridges")
    print(
        "rank | carrier | base | evidence_bonus | total | "
        "archive_hits | archive_files | history_subject_hits | history_path_hits"
    )
    for rank, (item, base, bonus, total, current_hits) in enumerate(ranked[:12], 1):
        files_hit = archive[item.name]["files"]
        subject_hits = len(history[item.name]["subject_hits"])
        path_hits = history[item.name]["path_hits"]
        print(
            f"{rank:02d} | {item.name} | {base} | {bonus} | {total} | "
            f"{current_hits} | {files_hit} | {subject_hits} | {path_hits}"
        )

    print_section("Tournament Analysis")
    print(f"vertices={fp['vertices']}")
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"scc_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_path_count_under_fixed_tie_path={fp['hamiltonian_path_count_under_fixed_tie_path']}")
    print(f"edge_flips_against_tie_path={fp['edge_flips_against_tie_path']}")
    print("hamiltonian_path=" + " -> ".join(fp["hamiltonian_path"]))

    print_section("Top bridge details")
    for rank, (item, base, bonus, total, current_hits) in enumerate(ranked[:10], 1):
        ev = archive[item.name]
        hist = history[item.name]
        print(f"\n### {rank:02d}. {item.name}")
        print(f"anchors={','.join(item.anchors)}")
        print(f"score_vector={dict(zip(AXES, item.scores))}")
        print(f"connection={item.connection}")
        print(f"guardrail={item.guardrail}")
        print(f"next_experiment={item.next_experiment}")
        print(f"current_hits={current_hits} current_files={ev['files']}")
        print(f"top_files={ev['top_files']}")
        print(
            f"history_subject_hits={len(hist['subject_hits'])} "
            f"history_path_hits={hist['path_hits']} "
            f"history_path_files={hist['path_files']}"
        )
        if ev["samples"]:
            print("samples:")
            for sample in ev["samples"]:
                print(f"  - {sample}")
        if hist["subject_hits"]:
            print("commit_subject_hits:")
            for sample in hist["subject_hits"][:4]:
                print(f"  - {sample}")
        if hist["path_samples"]:
            print("history_path_samples:")
            for sample in hist["path_samples"][:4]:
                print(f"  - {sample}")

    print_section("False friends and late-stage carriers")
    for item, base, bonus, total, current_hits in ranked[-5:]:
        print(f"- {item.name}: {item.guardrail}")

    print_section("Assumption Challenge")
    print(
        "considered_vertices=runners,gaps,fixed_circle_sections,section_boundaries,"
        "wall_crossing_events,residues,cover_arcs,Fourier_modes,matroid_circuits,"
        "u_fibers,mirror_pairs,free_hole_packets,owner_labels,proof_obligations"
    )
    print("chosen_vertices=niche_proof_carriers")
    print(
        "preserved_predicate=random031 terminal discharge after hard-seam deletion,"
        " with private-label firewall and seam-owner debt retained"
    )
    print(
        "destroyed_information=raw runner order, raw row names, raw witness counts,"
        " old analogy labels, and scalar maturity unless reconstructed by receiver fields"
    )
    print(
        "challenged_assumption=older niche topics should not be ranked by age,"
        " popularity, or scalar elegance; they should be ranked by which random031"
        " seam-complement coordinate they can legally restore"
    )


if __name__ == "__main__":
    main()
