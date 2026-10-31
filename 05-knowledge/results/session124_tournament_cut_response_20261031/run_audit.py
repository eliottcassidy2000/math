"""Replay both independent auditors normally, under -O, and after relocation."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
AUDITS = [("independent_check.py", "normal.json", "optimized.json"),
          ("independent_consequences.py", "consequences.normal.json", "consequences.optimized.json")]


def run(path, optimized=False, cwd=None):
    command = [sys.executable] + (["-O"] if optimized else []) + [str(path)]
    return subprocess.run(command, cwd=cwd, check=True, capture_output=True).stdout


def main():
    normal = {}
    hashes = {}
    for script, normal_name, optimized_name in AUDITS:
        payload = run(HERE / "audit" / script)
        optimized = run(HERE / "audit" / script, optimized=True)
        if payload != optimized:
            raise RuntimeError(("independent -O mismatch", script))
        for name, data in [(normal_name, payload), (optimized_name, optimized)]:
            (HERE / "audit" / name).write_bytes(data)
            hashes[name] = hashlib.sha256(data).hexdigest()
        normal[script] = payload
    with tempfile.TemporaryDirectory(prefix="cut-response-independent-") as temp:
        relocated = Path(temp) / "packet"
        (relocated / "audit").mkdir(parents=True)
        shutil.copyfile(HERE / "cut_response.py", relocated / "cut_response.py")
        for script, _, _ in AUDITS:
            shutil.copyfile(HERE / "audit" / script, relocated / "audit" / script)
        for script, _, _ in AUDITS:
            payload = run(relocated / "audit" / script, cwd=temp)
            if payload != normal[script]:
                raise RuntimeError(("independent relocation mismatch", script))
    summary = {"status": "PASS", "normal_optimized_identical": True,
               "unrelated_working_directory_replay_identical": True,
               "receipt_sha256": hashes}
    (HERE / "audit" / "verification.json").write_text(json.dumps(summary, sort_keys=True, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
