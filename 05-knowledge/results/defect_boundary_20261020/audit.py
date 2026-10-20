"""Verify the frozen packet; --replay also reruns every normal/-O control."""
import ast
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def need(value, message):
    if not value:
        raise RuntimeError(message)


def digest(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    manifest = json.loads((HERE / "manifest.json").read_text(encoding="utf-8"))
    for relative, expected in manifest["files"].items():
        need(digest(REPO / relative) == expected, "changed file: " + relative)
    print("Normalized SHA-256 checks:", len(manifest["files"]), "files; PASS",
          flush=True)
    for name in manifest["scripts"]:
        source = HERE / (name + ".py")
        tree = ast.parse(source.read_text(encoding="utf-8"))
        need(not any(isinstance(n, ast.Assert) for n in ast.walk(tree)),
             "disabled assertion: " + name)
        need((HERE / (name + ".out")).exists(), "missing output: " + name)
        if "--replay" in sys.argv:
            outputs = []
            for flags in ([], ["-O"]):
                run = subprocess.run([sys.executable, *flags, str(source)],
                                     capture_output=True)
                need(run.returncode == 0, name + ": " +
                     run.stdout.decode("utf-8", errors="replace") +
                     run.stderr.decode("utf-8", errors="replace"))
                need(not run.stderr, name + ": stderr")
                outputs.append(run.stdout.replace(b"\r\n", b"\n"))
            retained = (HERE / (name + ".out")).read_bytes().replace(b"\r\n", b"\n")
            need(outputs[0] == outputs[1] == retained,
                 "normal/-O/retained mismatch: " + name)
            print(name + ": normal = -O = retained; PASS", flush=True)
    print("All control scripts use active exception checks; PASS")
    print("PACKET AUDIT PASS")


if __name__ == "__main__":
    main()
