"""Reproduce exact mathematical receipts; only --refresh writes files."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
JOBS = [
    ("ports/check.py", "ports/normal.json", "ports/optimized.json"),
    ("ports/endcap_check.py", "ports/endcap-normal.json", "ports/endcap-optimized.json"),
    ("counting/check.py", "counting/normal.json", "counting/optimized.json"),
    ("counting/check_positive_jet.py", "counting/jet-normal.json", "counting/jet-optimized.json"),
    ("bridges/check.py", "bridges/normal.json", "bridges/optimized.json"),
    ("valuation/check.py", "valuation/normal.json", "valuation/optimized.json"),
    ("valuation/INDEPENDENT_REVIEW.md", "valuation/review-normal.json", "valuation/review-optimized.json"),
]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def run_job(job, refresh):
    source, normal, optimized = job
    path = HERE / source
    if path.suffix == ".md":
        code = path.read_text(encoding="utf-8").split("```python\n")[1].split("```", 1)[0]
        args, cwd = ["-c", code], HERE.parents[2]
    else:
        args, cwd = [str(path)], path.parent
    observations = []
    for flags, receipt in (([], normal), (["-O"], optimized)):
        result = subprocess.run([sys.executable, "-B", *flags, *args], cwd=cwd,
                                capture_output=True, text=True, timeout=240)
        need(result.returncode == 0, (source, result.returncode, result.stderr[-4000:]))
        value = json.loads(result.stdout)
        need(value.get("status") == "PASS", (source, "failed mathematical receipt"))
        output = json.dumps(value, indent=2, sort_keys=True) + "\n"
        observations.append(output)
        target = HERE / receipt
        if refresh:
            target.write_bytes(output.encode("utf-8"))
        else:
            need(target.read_text(encoding="utf-8") == output, (source, "receipt changed"))
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
