"""Run normal/-O checks and a relocated replay; no third-party packages."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
SCRIPTS = ("verify.py", "verify_consequences.py", "verify_examples.py")
SOURCES = ("cut_response.py",) + SCRIPTS


def run(script, optimized=False, cwd=None):
    command = [sys.executable] + (["-O"] if optimized else []) + [str(script)]
    result = subprocess.run(command, cwd=cwd, capture_output=True, check=True)
    receipt = json.loads(result.stdout)
    if receipt.get("status") != "PASS":
        raise RuntimeError(("nonpassing receipt", script, receipt))
    return result.stdout


def main():
    destination = HERE / "receipts"
    destination.mkdir(exist_ok=True)
    normal = {}
    outputs = {}
    for script in SCRIPTS:
        stem = Path(script).stem
        normal[script] = run(HERE / script)
        optimized = run(HERE / script, optimized=True)
        if normal[script] != optimized:
            raise RuntimeError(("normal/-O mismatch", script))
        for mode, payload in [("normal", normal[script]), ("optimized", optimized)]:
            name = f"{stem}.{mode}.json"
            (destination / name).write_bytes(payload)
            outputs[name] = hashlib.sha256(payload).hexdigest()
    with tempfile.TemporaryDirectory(prefix="cut-response-portability-") as temp:
        relocated = Path(temp) / "package"
        relocated.mkdir()
        for source in SOURCES:
            shutil.copyfile(HERE / source, relocated / source)
        for script in SCRIPTS:
            payload = run(relocated / script, cwd=temp)
            if payload != normal[script]:
                raise RuntimeError(("relocation mismatch", script))
    summary = {"status": "PASS", "normal_optimized_identical": True,
               "unrelated_working_directory_replay_identical": True,
               "scripts": list(SCRIPTS), "receipt_sha256": outputs}
    (destination / "verification.json").write_text(json.dumps(summary, sort_keys=True, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
