"""Replay all packet controls normally and with -O; compare literal stdout.

Run with --refresh only to replace the retained outputs after both runs agree.
"""
import ast
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def need(value, message):
    if not value:
        raise RuntimeError(message)


def main():
    for name in ("growth_check", "core_count", "allocation_check"):
        path = ROOT / (name + ".py")
        tree = ast.parse(path.read_text(encoding="utf-8"))
        need(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
             name + ": disabled assertion")
        results = []
        for flags in ([], ["-O"]):
            run = subprocess.run([sys.executable, *flags, str(path)],
                                 capture_output=True, text=True, encoding="utf-8")
            need(run.returncode == 0, name + ": " + run.stdout + run.stderr)
            need(not run.stderr, name + ": unexpected stderr " + run.stderr)
            results.append(run.stdout.replace("\r\n", "\n"))
        need(results[0] == results[1], name + ": normal/-O output differs")
        output = ROOT / (name + ".out")
        if "--refresh" in sys.argv:
            output.write_text(results[0], encoding="utf-8", newline="\n")
        need(output.read_text(encoding="utf-8").replace("\r\n", "\n") == results[0],
             name + ": retained output differs")
        print(name + ": normal = -O = retained stdout; no assertions; PASS", flush=True)
    print("PACKET REPLAY PASS")


if __name__ == "__main__":
    main()
