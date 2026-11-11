#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: the predictions printed in window_split_prereg_27_24.out (columns 'E pred',
'f_div pred', 'E_nd pred') against the registered ones (prereg_window_27_24.json and
prereg_window_27_24.predictions.out), and a recount of the test items and of the pass/fail of every interval from
the measured columns of the output (5 decimals).

Usage: python3 check_prereg_consistency.py DIR   (DIR = the session's experiment directory, read-only)
"""
import json
import os
import re
import sys

d = sys.argv[1]
spec = json.load(open(os.path.join(d, "prereg_window_27_24.json")))
print("json keys:", sorted(spec.keys()))
print("window:", spec["window"], "| predicted triples per odd N:", spec["predicted_triples_per_odd_N"])
pred_json = spec["predictions"]
# predictions.out
pred_out = {}
for line in open(os.path.join(d, "prereg_window_27_24.predictions.out")):
    m = re.match(r"\s+p =\s+(\d+) first divisible h =\s+(-?\d+)\s+s = 0: E ([+-][\d.]+)\s+f_div ([\d.]+)\s+E_nd ([+-][\d.]+)\s+s >= 1: E ([+-][\d.]+)\s+f_div ([\d.]+)\s+E_nd ([+-][\d.]+)", line)
    if m:
        p = int(m.group(1))
        pred_out[(p, 0)] = (float(m.group(3)), float(m.group(4)), float(m.group(5)))
        pred_out[(p, 1)] = (float(m.group(6)), float(m.group(7)), float(m.group(8)))
# run output
meas = {}
pred_run = {}
for line in open(os.path.join(d, "window_split_prereg_27_24.out")):
    f = line.split()
    if len(f) == 10 and f[0].isdigit() and f[1] in ("s=0", "s>=1"):
        p = int(f[0])
        sh = 0 if f[1] == "s=0" else 1
        meas[(p, sh)] = (float(f[2]), float(f[4]), float(f[6]), float(f[9]))
        pred_run[(p, sh)] = (float(f[3]), float(f[5]), float(f[7]))
print("primes in json:", len(pred_json) // 2, " in predictions.out:", len(pred_out) // 2, " in the run output:", len(pred_run) // 2)
bad = 0
for (p, sh), v in sorted(pred_run.items()):
    j = pred_json[f"{p}:{sh}"]
    a = (round(j["E"], 5), round(j["fdiv"], 6), round(j["End"], 5))
    b = pred_out[(p, sh)]
    if any(abs(x - y) > 1.01e-5 for x, y in zip(a, v)) or any(abs(x - y) > 1.01e-5 for x, y in zip(b, v)):
        bad += 1
        print("  MISMATCH", p, sh, "json", a, "predictions.out", b, "run", v)
print("predictions: run output vs json vs predictions.out, mismatches:", bad)
# tests
tests = spec["tests"]
print("tests in the json:", {k: len(v["items"]) if isinstance(v, dict) and "items" in v else type(v).__name__ for k, v in tests.items()} if isinstance(tests, dict) else type(tests))
tot = 0
fails = []
names = {}
if isinstance(tests, dict):
    it = tests.items()
else:
    it = [(t.get("name", str(i)), t) for i, t in enumerate(tests)]
for name, t in it:
    items = t["items"]
    tot += len(items)
    nchecked = 0
    for item in items:
        key = item["key"]
        parts = key.split(":")
        val = None
        if parts[0] == "class":
            p, sh, what = int(parts[1]), int(parts[2]), parts[3]
            if (p, sh) in meas:
                val = {"E": meas[(p, sh)][0], "fdiv": meas[(p, sh)][1], "End": meas[(p, sh)][2]}[what]
        if val is None:
            continue
        nchecked += 1
        if not (item["lo"] - 6e-6 <= val <= item["hi"] + 6e-6):
            fails.append((name, key, val, item["lo"], item["hi"]))
    names[name] = (len(items), nchecked)
print("items per test (total, re-checked here from the printed measured columns):", names)
print("total items:", tot, " failures among the re-checked items:", fails)
# the tightest items: margin relative to the half-width
rows = []
for name, t in it:
    for item in t["items"]:
        parts = item["key"].split(":")
        if parts[0] == "class":
            p, sh, what = int(parts[1]), int(parts[2]), parts[3]
            if (p, sh) in meas:
                val = {"E": meas[(p, sh)][0], "fdiv": meas[(p, sh)][1], "End": meas[(p, sh)][2]}[what]
                lo, hi = item["lo"], item["hi"]
                if hi > lo and hi < 9:
                    mid, half = (lo + hi) / 2, (hi - lo) / 2
                    rows.append((abs(val - mid) / half, name, item["key"], val, lo, hi))
rows.sort(reverse=True)
print("tightest items (|value - centre| / half-width):")
for r in rows[:8]:
    print("   %.3f  %s  %s  value %.5f in [%.5f, %.5f]" % r)
