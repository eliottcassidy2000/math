#!/usr/bin/env python3
"""corpus_archaeology_kps_S128c136.py -- kind-pasteur-2026-07-21-S128c136

MACHINE ARCHAEOLOGY of the 1392-theorem / 2554-reflection corpus: find FORGOTTEN threads
(theorems that were proved and then never built upon), map the topic zoo, and locate the gaps.

Builds the citation graph (who references whom via THM-XXXX mentions across every file),
computes in-degree, and flags:
  * FORGOTTEN: a theorem cited by <= 1 other file, i.e. a dead-end leaf -- a result the repo
    proved and forgot;
  * ORPHAN TOPICS: title keywords that appear in very few theorems (isolated ideas);
  * ERA GAPS: THM-number ranges that are dense vs sparse (abandoned eras).
"""
import os
import re
import sys
from collections import Counter, defaultdict

ROOT = "."
THMDIR = "01-canon/theorems"

# ---- 1. gather theorem ids + titles ----
titles = {}
for fn in os.listdir(THMDIR):
    if not fn.endswith(".md"):
        continue
    m = re.match(r"(THM|LEM|CONJ)-(\d+)", fn)
    if not m:
        continue
    tid = "%s-%s" % (m.group(1), m.group(2).lstrip("0") or "0")
    path = os.path.join(THMDIR, fn)
    try:
        with open(path, encoding="utf-8", errors="ignore") as f:
            head = f.read(1200)
    except Exception:
        continue
    tm = re.search(r'title:\s*"?(.+)', head)
    titles[fn] = (tid, (tm.group(1)[:90] if tm else fn))

# canonical id per file
fileid = {}
for fn, (tid, _) in titles.items():
    fileid[fn] = tid

# ---- 2. citation counts: how many DISTINCT files mention each THM-N ----
ref_files = defaultdict(set)
scan_dirs = ["01-canon/theorems", "07-reflections", "00-navigation"]
nfiles = 0
for d in scan_dirs:
    for fn in os.listdir(d):
        if not fn.endswith(".md"):
            continue
        nfiles += 1
        path = os.path.join(d, fn)
        try:
            with open(path, encoding="utf-8", errors="ignore") as f:
                txt = f.read()
        except Exception:
            continue
        for num in set(re.findall(r"THM-(\d+)", txt)):
            key = "THM-%s" % (num.lstrip("0") or "0")
            ref_files[key].add((d, fn))

print("scanned %d files for THM references" % nfiles)

# ---- 3. FORGOTTEN theorems: cited by <=1 file other than themselves ----
def own_num(fn):
    m = re.match(r"THM-(\d+)", fn)
    return "THM-%s" % (m.group(1).lstrip("0") or "0") if m else None

forgotten = []
for fn, (tid, title) in sorted(titles.items()):
    if not tid.startswith("THM"):
        continue
    num = int(re.match(r"THM-(\d+)", tid).group(1)) if re.match(r"THM-(\d+)", tid) else 0
    cites = ref_files.get(tid, set())
    # remove self-file
    external = {(d, f) for (d, f) in cites if f != fn}
    if len(external) <= 1 and num < 1700:      # old enough to have been buildable-on
        forgotten.append((len(external), num, tid, title))
forgotten.sort()
print()
print("=" * 92)
print("FORGOTTEN THEOREMS (cited by <=1 external file, THM# < 1700) -- proved then abandoned")
print("=" * 92)
print("  total flagged: %d" % len(forgotten))
for cnt, num, tid, title in forgotten[:40]:
    print("  [%d refs] %-9s %s" % (cnt, tid, title))

# ---- 4. TOPIC ZOO: keyword frequency in titles ----
KW = ["Paley", "arborescence", "Pfaffian", "homology", "Betti", "circulant", "Tutte",
      "eigenvalue", "spectral", "score", "domination", "king", "Hamiltonian", "3-cycle",
      "OCF", "Redei", "switching", "even graph", "Seidel", "nullcone", "GMC", "TNC",
      "LRC", "Jacobian", "Dixmier", "Kakeya", "moment", "binary form", "GIT", "Gauss sum",
      "quadratic residue", "BIBD", "Fibonacci", "Catalan", "Bessel", "Hermite", "Legendre",
      "Chebyshev", "Vandermonde", "permanent", "feedback", "median order", "condensation",
      "staircase", "tiling", "wiggly", "waggly", "blue", "metagraph", "Farey", "Smith",
      "Cayley", "octonion", "sedenion", "Lyapunov", "Ramanujan", "zeta", "modular"]
kwcount = Counter()
for fn, (tid, title) in titles.items():
    tl = title.lower()
    for k in KW:
        if k.lower() in tl:
            kwcount[k] += 1
print()
print("=" * 92)
print("TOPIC ZOO (title-keyword frequency) -- big = mainline, small = orphan/underexplored")
print("=" * 92)
for k, c in kwcount.most_common():
    bar = "#" * min(c, 40)
    print("  %-16s %3d  %s" % (k, c, bar))
orphans = [k for k in KW if kwcount[k] <= 2]
print()
print("  ORPHAN/near-orphan topics (<=2 theorems): %s" % ", ".join(orphans))

# ---- 5. ERA density: THM-number histogram in bins of 100 ----
print()
print("=" * 92)
print("ERA DENSITY (theorems per 100-number bin) -- sparse bins = abandoned eras / gaps")
print("=" * 92)
bins = Counter()
for fn, (tid, _) in titles.items():
    m = re.match(r"THM-(\d+)", tid)
    if m:
        bins[int(m.group(1)) // 100] += 1
for b in range(0, 19):
    c = bins.get(b, 0)
    print("  THM-%02d00..%02d99 : %3d  %s" % (b, b, c, "#" * min(c, 60)))
