"""Exact isomorphism catalog and theorem-bounded attachment searches."""
import argparse
import json
from collections import Counter
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=7)
    parser.add_argument("--slope", type=int, default=63)
    parser.add_argument("--intercept", type=int)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="tournament_catalog_") as tmp:
        exe = Path(tmp) / "catalog.exe"
        subprocess.run(["g++", "-O3", "-std=c++20", str(HERE / "catalog.cpp"), "-o", str(exe)], check=True)
        command = [str(exe), str(args.order), str(args.slope)]
        if args.intercept is not None:
            command.append(str(args.intercept))
        data = json.loads(subprocess.check_output(command, text=True))
    expected = [1, 1, 2, 4, 12, 56, 456, 6880]
    if data["counts"] != expected[:args.order]:
        raise RuntimeError(("Unexpected isomorphism census", data["counts"]))
    summaries = {}
    for target_h in (5, 29, 35):
        rows = [r for r in data["representatives"] if r["h"] == target_h and r["n"] <= 7]
        strong = [r for r in rows if r["strong"]]
        def odd_part(n):
            return n // (n & -n)
        summaries[str(target_h)] = {
            "all_count": len(rows), "strong_count": len(strong),
            "all_p2": sorted({r["p2"] for r in rows}),
            "strong_p2": sorted({r["p2"] for r in strong}),
            "strong_p2_odd_parts": sorted({odd_part(r["p2"]) for r in strong}),
            "strong_p2_odd_parts_dividing_63": sorted({odd_part(r["p2"]) for r in strong if 63 % odd_part(r["p2"]) == 0}),
        }
    data["prefix_summaries"] = summaries
    suffix = "" if args.intercept is None else f"_c{args.intercept}"
    destination = HERE / f"catalog{args.order}_f{args.slope}{suffix}.json"
    destination.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(destination), "counts": data["counts"], "hit_counts": dict(Counter(r["intercept"] for r in data["hits"])), "prefix_summaries": summaries}, indent=2))


if __name__ == "__main__":
    main()
