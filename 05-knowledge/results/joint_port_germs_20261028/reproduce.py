"""Reproduce this packet in an isolated temporary layout, including python -O."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
REL = HERE.relative_to(ROOT)
ORACLE = Path("05-knowledge/results/relation_tower_frobenius_ports_20261028/check_root_switch.py")
RUNS = (
    ("check_ports.py", "ports.json"),
    ("check_strong_tags.py", "strong_tags.json"),
    ("check_lrc_germs.py", "lrc_germs.json"),
)

def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()

def main():
    manifest = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["sha256_lf"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError("Source pin changed; review before resealing: " + name)
    print(json.dumps({"source_pins": len(manifest["sha256_lf"]), "status": "PASS"}), flush=True)
    with tempfile.TemporaryDirectory(prefix="math-port-germs-") as temporary:
        scratch = Path(temporary)
        packet = scratch / REL
        packet.mkdir(parents=True)
        for name in ["port_core.py"] + [a for a, b in RUNS]:
            shutil.copyfile(HERE / name, packet / name)
        (scratch / ORACLE).parent.mkdir(parents=True)
        shutil.copyfile(ROOT / ORACLE, scratch / ORACLE)
        for optimize in (False, True):
            checks = 0
            for script, receipt in RUNS:
                command = [sys.executable, "-B"]
                if optimize:
                    command.append("-O")
                command.append(str(packet / script))
                result = subprocess.run(command, cwd=scratch, capture_output=True,
                                        text=True, encoding="utf-8", timeout=300)
                if result.returncode:
                    raise RuntimeError(script + "\n" + result.stdout + "\n" + result.stderr)
                actual = json.loads((packet / receipt).read_text(encoding="utf-8"))
                expected = json.loads((HERE / receipt).read_text(encoding="utf-8"))
                if actual != expected:
                    raise RuntimeError("Receipt mismatch: " + receipt)
                checks += actual["checks"]
                print(json.dumps({"script": script, "optimized": optimize,
                                  "checks": actual["checks"], "status": "PASS"}), flush=True)
            print(json.dumps({"optimized": optimize, "total_checks": checks,
                              "receipts_equal": True}), flush=True)
    print("PASS: source pins and both exact replays")

if __name__ == "__main__":
    main()
