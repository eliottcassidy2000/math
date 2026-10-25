"""Reproduce saved mathematical receipts; --refresh writes them explicitly."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
JOBS = [
    ("geometry/check.py", "geometry/normal.json", "geometry/optimized.json"),
    ("ports/check.py", "ports/normal.json", "ports/optimized.json"),
    ("ports/boundary_check.py", "ports/boundary-normal.json", "ports/boundary-optimized.json"),
    ("hn/check.py", "hn/normal.json", "hn/optimized.json"),
    ("counting/check.py", "counting/normal.json", "counting/optimized.json"),
    ("counting/check_excess.py", "counting/excess-normal.json", "counting/excess-optimized.json"),
    ("counting/check_colored_inventory.py", "counting/colored-normal.json", "counting/colored-optimized.json"),
    ("counting/check_clone_involution.py", "counting/clone-normal.json", "counting/clone-optimized.json"),
    ("ports/INDEPENDENT_REVIEW.md", "ports/review-normal.json", "ports/review-optimized.json"),
    ("counting/INDEPENDENT_REVIEW.md", "counting/review-normal.out", "counting/review-optimized.out"),
    ("hn/RELATIONAL_LIFT.md", "hn/relation-normal.json", "hn/relation-optimized.json"),
]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def run_job(job, refresh):
    source, normal, optimized = job
    path = HERE / source
    if path.suffix == ".md":
        code = path.read_text(encoding="utf-8").split("```python\n")[1].split("```", 1)[0]
        args = ["-c", code]
        cwd = HERE.parents[2]
    else:
        args = [str(path)]
        cwd = path.parent
    observations = []
    for flags, receipt in (([], normal), (["-O"], optimized)):
        result = subprocess.run([sys.executable, "-B", *flags, *args], cwd=cwd,
                                capture_output=True, text=True, timeout=240)
        need(result.returncode == 0, (source, "child failed", result.returncode, result.stderr[-4000:]))
        output = result.stdout.replace("\r\n", "\n")
        if receipt.endswith(".json"):
            value = json.loads(output)
            need(value.get("status") == "PASS", (source, "failed mathematical receipt"))
            output = json.dumps(value, indent=2, sort_keys=True) + "\n"
        else:
            need(output.startswith("PASS"), (source, "missing PASS"))
        observations.append(output)
        target = HERE / receipt
        if refresh:
            target.write_bytes(output.encode("utf-8"))
        else:
            need(target.read_text(encoding="utf-8") == output, (source, "receipt changed", receipt))
    need(observations[0] == observations[1], (source, "normal/optimized mismatch"))
    return source


def main():
    refresh = "--refresh" in sys.argv
    with ThreadPoolExecutor(max_workers=3) as pool:
        for source in pool.map(lambda job: run_job(job, refresh), JOBS):
            print("PASS", source, flush=True)
    print(json.dumps({"status": "PASS", "paired_replays": len(JOBS), "refresh": refresh}))


if __name__ == "__main__":
    main()
