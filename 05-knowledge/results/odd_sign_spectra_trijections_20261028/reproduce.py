"""Reproduce all receipts without modifying the repository; normal and -O."""
from pathlib import Path
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NAMES = ("sign_spectrum", "paired_oddness", "lrc_reflection", "unit_distance")

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def normalized_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()

def comparable(obj):
    result = dict(obj)
    # Library version and a float diagnostic are not mathematical receipt fields.
    result.pop("sympy_version", None)
    result.pop("nonedge_squared_distance_gap_float_diagnostic", None)
    return result

def main():
    manifest = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    for relative, digest in manifest["sha256_lf"].items():
        need(normalized_hash(ROOT / relative) == digest, "source pin changed: " + relative)
    expected = {name: json.loads((HERE / (name + ".json")).read_text(encoding="utf-8")) for name in NAMES}
    runs = []
    with tempfile.TemporaryDirectory(prefix="odd-sign-spectra-") as directory:
        temp = Path(directory)
        for name in NAMES:
            shutil.copyfile(HERE / (name + ".py"), temp / (name + ".py"))
        for optimized in (False, True):
            for name in NAMES:
                command = [sys.executable, "-B"] + (["-O"] if optimized else []) + [str(temp / (name + ".py"))]
                result = subprocess.run(command, cwd=temp, capture_output=True, text=True, encoding="utf-8")
                need(result.returncode == 0, name + " failed: " + result.stderr + result.stdout)
                actual = json.loads((temp / (name + ".json")).read_text(encoding="utf-8"))
                need(comparable(actual) == comparable(expected[name]), "receipt differs: " + name)
                runs.append({"checker": name, "optimized": optimized, "checks": actual["checks"]})
                print(json.dumps(runs[-1]), flush=True)
    hostile_path = ROOT / "05-knowledge/results/relation_tower_frobenius_ports_20261028/check_root_switch.py"
    spec = importlib.util.spec_from_file_location("inherited_root_switch", hostile_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    witness = module.hostile()
    need(witness["automorphism_orders"] == [81, 27], "81/27 hostile groups")
    need(witness["labelled_cyclic_triples"] == 67, "81/27 common labelled triples")
    need(witness["identity_order_legal"] == [True, False], "81/27 named path separator")
    print(json.dumps({"status": "PASS", "pins": len(manifest["sha256_lf"]),
                      "runs": runs, "hostile_81_27": witness}, indent=2))

if __name__ == "__main__":
    main()
