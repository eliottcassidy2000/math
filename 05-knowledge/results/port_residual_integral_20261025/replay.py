"""Replay exact finite controls. Only --refresh writes receipts."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
JOBS = [
    ("h3/check.py", "h3/normal.json", "h3/optimized.json"),
    ("h4/check.py", "h4/normal.json", "h4/optimized.json"),
    ("h3/INDEPENDENT_REVIEW.md", "h3/review-normal.json", "h3/review-optimized.json"),
    ("PORT_FRONTIER.md", "port-frontier-normal.json", "port-frontier-optimized.json"),
    ("low_orders.py", "low_orders.normal.json", "low_orders.optimized.json"),
    ("LOW_ORDERS_REVIEW.md", "low-orders-review-normal.json", "low-orders-review-optimized.json"),
    ("counting/check.py", "counting/normal.json", "counting/optimized.json"),
    ("counting/INDEPENDENT_REVIEW.md", "counting/review-normal.json", "counting/review-optimized.json"),
    ("counting/modular.py", "counting/modular-normal.json", "counting/modular-optimized.json"),
    ("counting/MODULAR_REVIEW.md", "counting/modular-review-normal.json", "counting/modular-review-optimized.json"),
    ("counting/REVERSE_POTENTIAL.md", "counting/reverse-normal.json", "counting/reverse-optimized.json"),
    ("parity/check.py", "parity/normal.json", "parity/optimized.json"),
    ("parity/INDEPENDENT_REVIEW.md", "parity/review-normal.json", "parity/review-optimized.json"),
    ("precision/check.py", "precision/normal.json", "precision/optimized.json"),
    ("precision/INDEPENDENT_REVIEW.md", "precision/review-normal.json", "precision/review-optimized.json"),
]


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def run_job(job, refresh):
    source, normal, optimized = job
    path = HERE / source
    if path.suffix == ".md":
        blocks = re.findall(r"^```python\n(.*?)^```", path.read_text(encoding="utf-8"), re.M | re.S)
        need(len(blocks) == 1, (source, "expected exactly one executable Python block"))
        args = ["-c", blocks[0]]
    else:
        args = [str(path)]
    results = []
    for flags, receipt in (([], normal), (["-O"], optimized)):
        result = subprocess.run([sys.executable, "-B", *flags, *args], cwd=ROOT,
                                capture_output=True, text=True, timeout=300)
        need(result.returncode == 0, (source, result.returncode, result.stderr[-4000:]))
        if source == "counting/INDEPENDENT_REVIEW.md":
            need(result.stdout.strip() == "PASS: independent 105-cycle ceiling and mask-deletion controls",
                 (source, "unexpected independent control result"))
            value = {"status": "PASS", "control": result.stdout.strip()}
        else:
            value = json.loads(result.stdout)
        need(value.get("status") == "PASS", (source, "failed mathematical receipt"))
        output = json.dumps(value, indent=2, sort_keys=True) + "\n"
        results.append(output)
        target = HERE / receipt
        if refresh:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(output.encode("utf-8"))
        else:
            need(target.read_text(encoding="utf-8") == output, (source, "receipt changed"))
    need(results[0] == results[1], (source, "normal/optimized mismatch"))
    return source


def main():
    refresh = "--refresh" in sys.argv
    with ThreadPoolExecutor(max_workers=3) as pool:
        for source in pool.map(lambda job: run_job(job, refresh), JOBS):
            print("PASS", source, flush=True)
    print(json.dumps({"status": "PASS", "paired_replays": len(JOBS), "refresh": refresh}))


if __name__ == "__main__":
    main()
