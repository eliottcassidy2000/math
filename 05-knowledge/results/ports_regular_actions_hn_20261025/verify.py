"""Read-only source/receipt validation. --refresh explicitly replaces the pins."""
import ast
import hashlib
import json
from pathlib import Path
import re
import sys
from replay import JOBS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INHERITED = [
    "01-canon/theorems/THM-516-burnside-core-kernel-phi-reframe.md",
    "01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md",
    "05-knowledge/results/trijection_hourglass_ports_20261025/endpoint/GENERAL_TAIL.md",
    "05-knowledge/results/trijection_hourglass_ports_20261025/trijection/README.md",
    "05-knowledge/results/defect_boundary_20261020/boundary.md",
    "05-knowledge/results/tournament_enumeration_20261021a/README.md",
    "05-knowledge/results/circle_ports_polygon_palettes_20261021/README.md",
    "05-knowledge/results/session51_structural_errors_relational_packets_20261024/prime-offsets/THEOREM.md",
    "05-knowledge/results/six_seven_geometry_20261025/fano_codes/RESULT.md",
    "05-knowledge/results/six_seven_geometry_20261025/root_deletion/DELETION_AND_CUBE.md",
    "05-knowledge/results/six_seven_gale_cube_20261025/BUDGET_AND_GALE.md",
    "05-knowledge/results/session61_hadwiger_six_seven_20261025/torus-fano/TORUS_FANO_CAMPEDELLI.md",
    "05-knowledge/results/session61_hadwiger_six_seven_20261025/SOURCE_STATUS.md",
]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def digest(path):
    # Source documents are text; LF normalization makes checkout EOL immaterial.
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def main():
    excluded = {"manifest.json", "verification.json", "verification-optimized.json"}
    files = sorted(p for p in HERE.rglob("*") if p.is_file()
                   and "__pycache__" not in p.parts and p.name not in excluded
                   and (p.suffix in {".md", ".py", ".json", ".out"} or p.name == ".gitattributes"))
    manifest = HERE / "manifest.json"
    if "--refresh" in sys.argv:
        data = {"hash_semantics": "SHA256 of UTF-8 source bytes with CRLF normalized to LF",
                "packet": {p.relative_to(HERE).as_posix(): digest(p) for p in files},
                "inherited": {name: digest(ROOT / name) for name in INHERITED}}
        manifest.write_bytes((json.dumps(data, indent=2, sort_keys=True) + "\n").encode())
    data = json.loads(manifest.read_text(encoding="utf-8"))
    need(set(data["packet"]) == {p.relative_to(HERE).as_posix() for p in files}, "packet file set changed")
    for name, expected in data["packet"].items():
        need(digest(HERE / name) == expected, ("packet hash", name))
    need(set(data["inherited"]) == set(INHERITED), "inherited file set changed")
    for name, expected in data["inherited"].items():
        need(digest(ROOT / name) == expected, ("inherited hash", name))
    for source, normal, optimized in JOBS:
        a = (HERE / normal).read_text(encoding="utf-8")
        b = (HERE / optimized).read_text(encoding="utf-8")
        need(a == b, ("receipt mismatch", source))
        if normal.endswith(".json"):
            need(json.loads(a).get("status") == "PASS", ("receipt failed", source))
        else:
            need(a.startswith("PASS"), ("receipt failed", source))
    links = parsed = 0
    for path in files:
        if path.suffix == ".py":
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            parsed += 1
        if path.suffix != ".md":
            continue
        for url in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if "://" in url or url.startswith("#"):
                continue
            target = url.split("#", 1)[0].strip("<>")
            need((path.parent / target).exists(), ("broken link", str(path), url))
            links += 1
    print(json.dumps({"status": "PASS", "packet_files": len(files),
                      "inherited_files": len(INHERITED), "receipt_pairs": len(JOBS),
                      "python_files": parsed, "local_links": links}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
