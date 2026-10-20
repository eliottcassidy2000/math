"""Pin this packet, its imported owner machinery, and exact replay receipts."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def normalized(data):
    return data.replace(b"\r\n", b"\n")


def inputs():
    own = [p for p in HERE.iterdir() if p.is_file() and
           p.suffix in (".py", ".md", ".out", ".json") and
           p.name not in ("manifest.json", "verify.out")]
    deps = [ROOT / ("05-knowledge/results/" + p) for p in (
        "rotation_owners_spectrum_20261020e/owners_check.py",
        "rotation_owners_spectrum_20261020e/owners.md",
        "deficiency_three_20261020/ownership_check.py",
        "deficiency_three_20261020/ownership.md",
        "defect_boundary_20261020/sparse_count.py",
        "defect_boundary_20261020/sparse_theory.md",
        "metagraph_growth_20261020/enumeration.py",
        "metagraph_growth_20261020/growth_check.py",
        "spectrum_basins_20261020e/outer.md",
        "spectrum_basins_20261020e/pure_hubs.md",
        "complement_path_games_20261020/README.md",
        "complement_path_games_20261020/AUDIT_AND_PROVENANCE.md",
    )]
    return sorted(set(own + deps))


def digest(path):
    return hashlib.sha256(normalized(path.read_bytes())).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    manifest = HERE / "manifest.json"
    current = {p.relative_to(ROOT).as_posix(): digest(p) for p in inputs()}
    if args.record:
        manifest.write_text(json.dumps({"normalization": "CRLF to LF",
                                       "files": current}, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    saved = json.loads(manifest.read_text(encoding="utf-8"))
    require(current == saved["files"], "manifest changed: inspect before rerecording")
    print("SHA256 source/input/receipt coverage:", len(current), "PASS", flush=True)
    scripts = sorted(p for p in HERE.glob("*.py") if p.name != "verify.py")
    for path in scripts:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
                f"removable assertion in {path.name}")
        require(path.with_suffix(".out").exists(), f"missing receipt for {path.name}")
    links = 0
    for path in HERE.glob("*.md"):
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if "://" in target or target.startswith("#"):
                continue
            require((path.parent / target.split("#")[0]).resolve().exists(),
                    f"broken local link in {path.name}: {target}")
            links += 1
    print("Scripts parse without removable assertions:", len(scripts), flush=True)
    print("Local Markdown link targets:", links, "PASS", flush=True)
    if args.replay:
        for path in scripts:
            expected = normalized(path.with_suffix(".out").read_bytes())
            for flags in ([], ["-O"]):
                result = subprocess.run([sys.executable, *flags, str(path)], cwd=ROOT,
                                        capture_output=True, check=False)
                require(result.returncode == 0,
                        f"{path.name} {flags}: {result.stderr.decode(errors='replace')}")
                require(normalized(result.stdout) == expected,
                        f"receipt mismatch: {path.name} {flags}")
            print("Normal/optimized replay:", path.name, "PASS", flush=True)
    print("PASS", flush=True)


if __name__ == "__main__":
    main()
