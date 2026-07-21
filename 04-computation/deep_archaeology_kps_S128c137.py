#!/usr/bin/env python3
"""deep_archaeology_kps_S128c137.py -- kind-pasteur-2026-07-21-S128c137

DEEP archaeology, below the theorem layer, for NICHE forgotten ideas.  c136 found forgotten
theorems by citation.  This mines the layers underneath, where the idea-richest and most-lost
material lives:
  (A) OEIS SEQUENCES: every A-number ever mentioned, with its occurrence count.  A-numbers
      mentioned ONCE are concrete, checkable, niche forgotten connections (someone computed a
      sequence, matched it, and never returned).
  (B) RARE NAMED OBJECTS: specific polynomials/constants/named theorems mentioned once or twice.
  (C) CONFIRMED-BUT-UNCANONIZED HYPS: hypotheses marked CONFIRMED/PROVED whose id never appears
      in any theorem file -- proven ideas that never became theorems.
"""
import os
import re
import sys
from collections import Counter, defaultdict

DIRS = ["01-canon/theorems", "07-reflections", "00-navigation", "05-knowledge/hypotheses"]

# gather all text once
docs = {}
for d in DIRS:
    if not os.path.isdir(d):
        continue
    for fn in os.listdir(d):
        if fn.endswith(".md"):
            try:
                with open(os.path.join(d, fn), encoding="utf-8", errors="ignore") as f:
                    docs[(d, fn)] = f.read()
            except Exception:
                pass
print("loaded %d markdown docs" % len(docs))

# ---------- (A) OEIS sequences ----------
oeis_files = defaultdict(set)
oeis_ctx = {}
for (d, fn), txt in docs.items():
    for m in re.finditer(r"A(\d{6})", txt):
        a = "A" + m.group(1)
        oeis_files[a].add(fn)
        if a not in oeis_ctx:
            s = max(0, m.start() - 45)
            oeis_ctx[a] = txt[s:m.end() + 25].replace("\n", " ")
print()
print("=" * 92)
print("(A) OEIS SEQUENCES -- singletons are niche forgotten numerical connections")
print("=" * 92)
print("  total distinct A-numbers referenced: %d" % len(oeis_files))
singles = sorted([a for a, fs in oeis_files.items() if len(fs) == 1])
print("  mentioned in exactly ONE file (niche/forgotten): %d" % len(singles))
for a in singles[:34]:
    print("    %-9s  ...%s..." % (a, oeis_ctx[a].strip()[:70]))

# ---------- (B) rare named objects ----------
# specific patterns: named polynomials W_5, P(z)=..., special constants, "conjecture" names
NAMED = [
    r"Kelvin", r"Watson", r"Airy", r"Bessel", r"Weierstrass", r"Eisenstein", r"Dedekind",
    r"Bernoulli", r"Euler-Mascheroni", r"Liouville", r"Mahler", r"Salem", r"Pisot",
    r"Markoff", r"Apery", r"Apéry", r"Somos", r"Stern-Brocot", r"Thue-Morse",
    r"Beatty", r"Sturmian", r"continued fraction", r"Collatz", r"Sylvester",
    r"Cayley-Menger", r"Gram matrix", r"Coxeter", r"Weyl", r"root system", r"E8", r"Leech",
    r"quaternion", r"octonion", r"sedenion", r"Hurwitz", r"Clifford", r"Grassmann",
    r"tropical", r"amoeba", r"Newton polygon", r"Puiseux", r"resurgence", r"Stokes",
    r"Borel", r"hypergeometric", r"Meijer", r"Selberg", r"Macdonald", r"Jack polynomial",
    r"Schur", r"symmetric function", r"quasisymmetric", r"shuffle", r"Rota-Baxter",
    r"umbral", r"Sheffer", r"Appell", r"Riordan",
]
print()
print("=" * 92)
print("(B) RARE NAMED OBJECTS -- appearing in few files (thin threads worth reviving)")
print("=" * 92)
namecount = Counter()
namefiles = defaultdict(set)
for (d, fn), txt in docs.items():
    for pat in NAMED:
        if re.search(pat, txt, re.IGNORECASE):
            namecount[pat] += 1
            namefiles[pat].add(fn)
for pat, c in sorted(namecount.items(), key=lambda kv: kv[1]):
    if c <= 4:
        print("    %-22s %d files" % (pat.replace("\\", ""), c))

# ---------- (C) confirmed HYPs never canonized ----------
print()
print("=" * 92)
print("(C) CONFIRMED/PROVED HYPS whose id never appears in any THEOREM file (uncanonized)")
print("=" * 92)
hyp_txt = ""
for (d, fn), txt in docs.items():
    if d == "05-knowledge/hypotheses":
        hyp_txt += txt
# find HYP blocks with a status of CONFIRMED/PROVED
thm_txt = "".join(txt for (d, fn), txt in docs.items() if d == "01-canon/theorems")
thm_ids = set(re.findall(r"HYP-(\d+)", thm_txt))
uncanon = []
for m in re.finditer(r"HYP-(\d+)[^\n]*", hyp_txt):
    hid = m.group(1)
    seg = hyp_txt[m.start():m.start() + 260]
    if re.search(r"CONFIRMED|PROVED", seg) and hid not in thm_ids:
        # keep the headline
        head = seg.split("\n")[0][:96]
        uncanon.append((hid, head))
# dedup by hid
seen = set()
shown = 0
for hid, head in uncanon:
    if hid in seen:
        continue
    seen.add(hid)
    print("    HYP-%s : %s" % (hid, head[head.find("—") + 1:].strip()[:78] if "—" in head else head[:78]))
    shown += 1
    if shown >= 26:
        break
print("  (total distinct confirmed-uncanonized HYPs flagged: %d)" % len(seen))
