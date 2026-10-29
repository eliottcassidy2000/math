"""Refresh explicit source pins after reviewing every intended change."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INHERITED = [
    "05-knowledge/results/relation_tower_frobenius_ports_20261028/check_root_switch.py",
    "01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md",
    "05-knowledge/results/session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md",
    "05-knowledge/results/session89_consecutive_tournament_pairs_20261028/counting/LOW_BITS.md",
    "05-knowledge/results/session90_even_graphs_odd_descent_20261029/even/DECORATED_EVEN_GRAPH_DUALITY.md",
]
paths = sorted([p for p in HERE.iterdir() if p.suffix in (".py", ".md", ".json") and p.name != "manifest.json"]
               + [ROOT / relative for relative in INHERITED])
manifest = {"normalization": "SHA256 after CRLF -> LF; no other normalization",
            "source_snapshot_before_concurrent_integration": "3eb49fa63d7ff4e58742b1ae84944f832e87a118",
            "sha256_lf": {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                          for p in paths}}
(HERE / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"pinned_files": len(paths)}))
