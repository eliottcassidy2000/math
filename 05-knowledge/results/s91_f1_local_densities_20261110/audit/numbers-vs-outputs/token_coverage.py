#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: blunt coverage check.  Every decimal number (at least two decimals) and every integer
of at least three digits in the audited documents is searched in the committed text outputs of the experiment
directory: a document number with k decimals is 'found' if some number printed in an output rounds or truncates to
it at k decimals.  Unmatched tokens are listed for manual review (they may be derived values, citations, etc.).

Usage: python3 coverage.py EXPERIMENT_DIR SNAPSHOT_DIR
"""
import os
import re
import subprocess
import sys

exp, snap = sys.argv[1], sys.argv[2]
tracked = [os.path.relpath(os.path.join(dp, f), exp) for dp, dn, fn in os.walk(exp) for f in fn]   # exp = git archive of HEAD
outs = [f for f in tracked if f.endswith((".out", ".log", ".txt"))]
num_re = re.compile(r"(?<![\w.])[+-]?\d+\.\d+(?:[eE][+-]?\d+)?|(?<![\w.])\d+(?![\w.]*\d)")
floats = []
ints = set()
for f in outs:
    try:
        txt = open(os.path.join(exp, f), encoding="utf-8", errors="replace").read()
    except OSError:
        continue
    for m in re.finditer(r"[+-]?\d+\.\d+(?:[eE][+-]?\d+)?", txt):
        try:
            floats.append(abs(float(m.group(0))))
        except ValueError:
            pass
    for m in re.finditer(r"\d+", txt):
        if len(m.group(0)) <= 12:
            ints.add(int(m.group(0)))
floats = sorted(set(floats))
import bisect


def found_float(tok):
    t = abs(float(tok))
    k = len(tok.split(".")[1])
    w = 10 ** (-k)
    lo = bisect.bisect_left(floats, t - w)
    hi = bisect.bisect_right(floats, t + w)
    for v in floats[lo:hi]:
        if abs(round(v, k) - t) < 1e-12 or abs(int(v * 10 ** k + 1e-9) / 10 ** k - t) < 1e-12:
            return True
    return False


for doc, sections in (("THM-5129.md", None), ("HYP-9660.md", None), ("HYP-9661.md", None), ("THM-5128.md", None), ("README.md", (33, 304))):
    lines = open(os.path.join(snap, doc), encoding="utf-8").read().split("\n")
    if sections:
        lines = lines[sections[0] - 1:sections[1]]
        base = sections[0]
    else:
        base = 1
    miss_f = []
    miss_i = []
    nf = ni = 0
    for i, line in enumerate(lines):
        for m in re.finditer(r"(?<![\w.^/-])[+-]?\d+\.\d\d+(?![\w])", line):
            tok = m.group(0).lstrip("+-")
            nf += 1
            if not found_float(tok):
                miss_f.append((base + i, m.group(0)))
        for m in re.finditer(r"(?<![\w.^-])\d{1,3}(?:,\d{3})+(?![\w.])|(?<![\w.^-])\d{3,}(?![\w.,]\d|\w)", line):
            tok = m.group(0).replace(",", "")
            ni += 1
            if int(tok) not in ints:
                miss_i.append((base + i, m.group(0)))
    print(f"== {doc}: {nf} decimal tokens, {len(miss_f)} not matched; {ni} integer tokens (>= 3 digits), {len(miss_i)} not matched")
    print("   decimals not matched (line, token):", miss_f)
    print("   integers not matched (line, token):", miss_i)
