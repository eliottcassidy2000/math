#!/usr/bin/env python3
"""kind-pasteur-S128c91 -- HYP-8000 analysis: parse the extinction-hunt output,
stratify c(p) by residue class mod 14, fit within-class growth, and summarize
the extinction frontier + the k=12 retrodiction."""
import re, sys, io

FILES = ["05-knowledge/results/lrc14_extinction_hunt_kps_S128c91.out",
         "05-knowledge/results/lrc14_extinction_exact_kps_S128c91.out",
         "05-knowledge/results/lrc14_extinction_decide_kps_S128c91.out",
         "05-knowledge/results/lrc14_extinction_control_kps_S128c91.out"]
s = ""
for f in FILES:
    try: s += io.open(f, encoding='utf-8', errors='replace').read() + "
"
    except OSError: pass

rows = []
for m in re.finditer(r"p=(\d+) dk=(\d+) bound=(\d+) c=(\d+)", s):
    rows.append((int(m.group(1)), int(m.group(2)), int(m.group(4))))
print(f"exact c(p) rows: {len(rows)}")

# also fold in the HYP-7955 small-p values
small = [(29,2,7),(43,3,9),(61,4,9),(71,5,9),(101,7,10),(113,8,10),(127,9,10),(151,10,11),(173,12,11)]
allrows = sorted(set(small + rows))

print("\nresidue stratification (p mod 14 : rows p->c):")
strat = {}
for p, dk, c in allrows:
    strat.setdefault(p % 14, []).append((p, c))
for r in sorted(strat):
    seq = "  ".join(f"{p}->{c}" for p, c in strat[r])
    print(f"  p == {r:2d} (h/dk -> {7*14/(14 - (14 - r) % 14) if r else 0:.2f}-ish): {seq}")

print("\njump points (first p with each c value, overall):")
seen = {}
for p, dk, c in allrows:
    if c not in seen: seen[c] = p
for c in sorted(seen): print(f"  c = {c:2d} first at p = {seen[c]}")

ext = re.findall(r"p=(\d+) EXTINCT \(no <=13", s)
if ext:
    print(f"\nk=13 EXTINCT primes: {', '.join(ext[:20])}{'...' if len(ext) > 20 else ''}")
    print(f"FIRST EXTINCT: {ext[0]}")
else:
    print("\nno k=13 extinct primes found (yet) in output")

unk = re.findall(r"p=(\d+) (?:dk=\d+ bound=\d+ c=UNKNOWN|UNKNOWN)", s)
if unk: print(f"UNKNOWN (budget) at: {', '.join(unk[:15])}")

ext12 = re.findall(r"p=(\d+) EXTINCT for k=12", s)
print(f"\nk=12 control extinct primes in [167,733]: {len(ext12)}"
      f"{': ' + ', '.join(ext12[:20]) if ext12 else ' (none -- consistent with S-T nonempty grind)'}")
