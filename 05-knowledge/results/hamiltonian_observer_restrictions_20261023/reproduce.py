"""Replay scoped exact receipts and inherited-source hashes; requires SymPy."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    provenance = json.loads((HERE / "SOURCE_PROVENANCE.json").read_text(encoding="utf-8"))
    for name, expected in provenance["inherited_sha256"].items():
        if digest(ROOT / name) != expected:
            raise RuntimeError(f"Inherited input changed; inspect before replay: {name}")
    jobs = [
        ("coupling/check.py", "coupling/receipt.json"),
        ("coupling/check_native_descent.py", "coupling/native_descent_receipt.json"),
        ("compiler/check.py", "compiler/receipt.json"),
        ("history/check_gain_lift.py", "history/gain_lift_receipt.json"),
        ("trace/check_trace.py", "trace/trace_receipt.json"),
    ]
    mode = ["-O"] if "--optimized" in sys.argv else []
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONIOENCODING="utf-8")
    for program, receipt in jobs:
        before = digest(HERE / receipt)
        result = subprocess.run([sys.executable, *mode, "-X", "utf8", "-B", str(HERE / program)],
                                cwd=ROOT, env=env, text=True, encoding="utf-8",
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            print(result.stdout)
            raise RuntimeError(f"Checker failed: {program}")
        if digest(HERE / receipt) != before:
            raise RuntimeError(f"Receipt changed during replay: {receipt}")
        print(f"PASS {program}: frozen receipt unchanged", flush=True)
    print(f"PASS {len(jobs)} exact checkers; inherited sources unchanged; mode={mode or ['normal']}")


if __name__ == "__main__":
    main()
