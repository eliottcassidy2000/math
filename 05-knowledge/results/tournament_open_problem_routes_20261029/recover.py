"""Reproducible lexical discovery; hit counts are not mathematical validation."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PATTERNS = {
    "erdos_hajnal": r"Erd[oöőős]*[–-]Hajnal|Erd.s.Hajnal",
    "sumner_tree_embedding": r"Sumner|Havet.{0,3}Thomass|unavoidab.*tree",
    "hamiltonian_spectrum": r"H.spectrum|Hamiltonian.path.*spectrum|two.prime.lane|Chan.{0,3}Pak",
    "hamiltonian_extremum": r"Szele|maximum.{0,30}Hamilton|A038375",
    "cycle_packing": r"Bermond|Lichiardopol|Erd.s.P.sa|Tuza",
    "seymour_neighborhood": r"Seymour.{0,35}(neigh|conject)|second.neigh",
    "caccetta_haggkvist": r"Caccetta|H.ggkvist",
    "ramsey": r"tournament.Ramsey|Ramsey.tournament|Erd.s.Moser",
    "reconstruction": r"Stockmeyer|tournament.reconstruction|reconstruct.*tournament",
    "symmetric_functions": r"R.dei.Berge|Stanley.Stembridge|noncommut.*deletion",
    "density_decidability": r"tournament.{0,40}(undecid|quasirandom)|density.inequalit",
    "domination": r"Sch.tte|domination.number|heroes.in.tournament",
    "hadamard_designs": r"skew.Hadamard|doubly.regular.tournament|conference.matrices",
    "determinant_bound": r"H\s*(?:>=|≥)\s*disc|H.ge.disc|rooted.Pfaffian",
    "metagraph_closure": r"selected.path.*(conject|recover)|black.core.connect|joint.port.*(carrier|compos)",
    "coloring": r"Hadwiger.Nelson|Alon.Tarsi",
    "graceful_square_sum": r"graceful.tree|square.sum.conject|square.sum.*Hamilton",
    "lonely_runner_transfer": r"improvement.tournament|winding.tournament|tournament.*lonel|lonel.*tournament",
    "arithmetic_dynamics": r"tournament.{0,50}(Collatz|Mahler|Turing)|(?:Collatz|Mahler|Turing).{0,50}tournament",
    "jacobian_transfer": r"tournament.{0,60}Jacobian|Jacobian.{0,60}tournament",
}

def main():
    args = ["rg", "--files", "00-navigation", "01-canon", "05-knowledge", "07-reflections", "-g", "*.md"]
    paths = subprocess.check_output(args, cwd=ROOT, text=True, encoding="utf-8").splitlines()
    # Frozen nested checkout copies duplicate discovery and are not current truth.
    excluded = {"inputs", "clean-repo", "session65-repo", "node_modules", ".git"}
    compiled = {k: re.compile(v, re.I) for k, v in PATTERNS.items()}
    rows = {k: [] for k in compiled}
    digest_seen = set()
    scanned = duplicates = 0
    for name in sorted(paths):
        p = Path(name)
        if excluded.intersection(p.parts) or p.parent == HERE.relative_to(ROOT):
            continue
        absolute = ROOT / p
        if os.name == "nt":
            absolute = Path("\\\\?\\" + str(absolute))
        raw = absolute.read_bytes()
        sha = hashlib.sha256(raw).hexdigest()
        if sha in digest_seen:
            duplicates += 1
            continue
        digest_seen.add(sha)
        scanned += 1
        lines = raw.decode("utf-8", errors="replace").splitlines()
        for key, pat in compiled.items():
            hits = [(i, line.strip()[:350]) for i, line in enumerate(lines, 1) if pat.search(line)]
            if hits:
                rows[key].append({"path": p.as_posix(), "sha256": sha, "hit_count": len(hits), "first_hits": hits[:3]})
    data = {"base_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "scope": "Distinct Markdown contents in navigation, canon, knowledge and reflections; lexical candidate discovery, not exhaustive semantic reading.",
            "patterns": PATTERNS, "scanned_files": scanned, "duplicate_contents_skipped": duplicates,
            "categories": rows}
    (HERE / "discovery.json").write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"scanned_files": scanned, "duplicates": duplicates, "candidate_files": {k: len(v) for k, v in rows.items()}}, indent=2))

if __name__ == "__main__":
    main()
