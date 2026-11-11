#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: the table of README section 1 and the list of THM-5128 (G) against
window_split_prereg_27_24.out, twisted_spectrum.out and divisible_states.out (committed copies, read-only)."""
import re
import sys

exp = sys.argv[1]
snap = sys.argv[2]
meas = {}
for line in open(exp + "/window_split_prereg_27_24.out"):
    f = line.split()
    if len(f) == 10 and f[0].isdigit() and f[1] in ("s=0", "s>=1"):
        meas[(int(f[0]), 0 if f[1] == "s=0" else 1)] = (float(f[2]), float(f[3]), float(f[9]))
rates = {}
for line in open(exp + "/twisted_spectrum.out"):
    f = line.split()
    if len(f) > 4 and f[0].isdigit():
        rates[int(f[0])] = float(f[1])
first = {}
for line in open(exp + "/divisible_states.out"):
    f = line.split()
    if len(f) >= 9 and f[0].isdigit():
        first[int(f[0])] = int(f[1])

def close(doc, val):
    """doc is a decimal string; val the 5-decimal output value: doc must be a correct rounding of some number
    within 5e-6 of val"""
    nd = len(doc.split(".")[1])
    d = float(doc)
    return abs(d - val) <= 0.5 * 10 ** (-nd) + 5.1e-6

bad = 0
n = 0
for line in open(snap + "/README.md"):
    m = re.match(r"\| (\d+) \| (\d+) \| ([\d.]+|-) \| ([+-][\d.]+) / ([+-][\d.]+) \| ([+-][\d.]+) / ([+-][\d.]+) \| ([\d.]+) \|", line)
    if not m:
        continue
    p = int(m.group(1))
    n += 1
    msgs = []
    fs = int(m.group(2))
    ref_first = first.get(p, {65537: 32768}.get(p))
    if ref_first != fs:
        msgs.append(f"first divisible state {fs} vs output {ref_first}")
    if m.group(3) != "-" and not close(m.group(3), rates[p]):
        msgs.append(f"rate {m.group(3)} vs {rates[p]}")
    for sh, (a, b) in ((0, (m.group(4), m.group(5))), (1, (m.group(6), m.group(7)))):
        if not close(a, meas[(p, sh)][0]):
            msgs.append(f"sheet {sh} measured {a} vs {meas[(p, sh)][0]}")
        if not close(b, meas[(p, sh)][1]):
            msgs.append(f"sheet {sh} predicted {b} vs {meas[(p, sh)][1]}")
    g = p / (p * p - 1.0)
    nd = len(m.group(8).split(".")[1])
    if abs(float(m.group(8)) - g) > 0.5 * 10 ** (-nd) + 1e-12:
        msgs.append(f"generic limit {m.group(8)} vs p/(p^2-1) = {g:.6f}")
    print(f"README table p = {p}: " + ("ok" if not msgs else "; ".join(msgs)))
    bad += len(msgs)
print("rows:", n, "mismatches:", bad)
# THM-5128 (G)
txt = open(snap + "/THM-5128.md").read()
g = txt[txt.index("**(G) Windows"):txt.index("## Proofs")]
items = re.findall(r"`(?:p = )?(\d+)`: `([+-][\d.]+) / ([+-][\d.]+)`, `([+-][\d.]+) / ([+-][\d.]+)`", g)
bad2 = 0
for p, a, b, c, d in items:
    p = int(p)
    msgs = []
    for sh, (x, y) in ((0, (a, b)), (1, (c, d))):
        if not close(x, meas[(p, sh)][0]):
            msgs.append(f"sheet {sh} measured {x} vs {meas[(p, sh)][0]}")
        if not close(y, meas[(p, sh)][1]):
            msgs.append(f"sheet {sh} predicted {y} vs {meas[(p, sh)][1]}")
    print(f"THM-5128 (G) p = {p}: " + ("ok" if not msgs else "; ".join(msgs)))
    bad2 += len(msgs)
print("items:", len(items), "mismatches:", bad2)
