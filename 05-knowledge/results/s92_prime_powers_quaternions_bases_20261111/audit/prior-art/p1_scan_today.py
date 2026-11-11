#!/usr/bin/env python3
"""Auditor prior-art (S92).  Scan the text files changed on origin/main since 2026-11-10 20:00 (list produced by a
read-only git log) for phrases that would be prior statements of the things S92 types NEW.
Usage: python3 p1_scan_today.py <file-with-paths> [root]
Prints, per pattern, the files that match and up to 4 matching lines each.  Reads only."""
import io
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
listing = sys.argv[1]
root = sys.argv[2] if len(sys.argv) > 2 else "C:/tmp/wt-s92"
EXCLUDE = ("s92_prime_powers", "THM-5134", "THM-5135", "HYP-9665", "HYP-9666", "agents/")
PATTERNS = {
    "four-square / Jacobi r_4": r"four[- ]square|r_4\b|8 ?sigma\(n\)|8σ",
    "Ingham": r"Ingham",
    "quaternion": r"quaternion",
    "E8 / octonion order": r"\bE8\b|E_8|integral octonion",
    "prime-power factor": r"sigma\(p\^|prime[- ]power (law|factor|modul)|odd prime powers",
    "chi_4 / signed label": r"chi_4|\(m ?- ?s\)/4|\|4n ?\+ ?1\||signed (integer|edge|source)",
    "self-converse strong": r"self-converse",
    "strong components as primes": r"as primes|strong components? (are|as) (the )?(primes|irreducible)",
    "arc usage": r"arc[- ]usage|edge[- ]usage|using the pair consecutively|consecutive(ly)? in",
    "abundancy / sigma(n)/n": r"abundanc|sigma\(n\)/n",
    "additive square / convolution of table": r"additive (square|convolution)|N_m \* N_m|sum of two (table )?entries|two products",
    "bicomplex / C4xC2": r"bicomplex|tessarine|C_?4 ?x ?C_?2|Z_?4 ?x ?Z_?2|split[- ]complex",
    "fan primes / alternation": r"alternation prime|fan prime|2\^t ?\+ ?1",
    "29 arcs / five cycles": r"29 arcs|five (known )?cycles|18 (arcs|labels)",
    "generating function strong": r"1/\(1 ?- ?\(x|6 ?x\^5|sequences of strong|free monoid",
    "digit bound / numeral": r"numeral|digit bound|no carry|carry-free",
    "balanced state base 3": r"X ?= ?y|balanced state",
}
paths = []
with io.open(listing, encoding="utf-8", errors="replace") as fh:
    for line in fh:
        p = line.strip()
        if p and not any(e in p for e in EXCLUDE):
            paths.append(p)
texts = {}
for p in paths:
    full = os.path.join(root, p)
    if os.path.isfile(full):
        try:
            with io.open(full, encoding="utf-8", errors="replace") as fh:
                texts[p] = fh.read().split("\n")
        except OSError:
            pass
print("files scanned:", len(texts), "of", len(paths))
for name, pat in PATTERNS.items():
    rx = re.compile(pat, re.I)
    print("=" * 100)
    print("PATTERN:", name, "   regex:", pat)
    nfile = 0
    for p, lines in texts.items():
        hits = [(i + 1, l) for i, l in enumerate(lines) if rx.search(l)]
        if hits:
            nfile += 1
            if nfile > 40:
                continue
            print("  %s  (%d lines)" % (p, len(hits)))
            for i, l in hits[:4]:
                print("     %5d: %s" % (i, l.strip()[:230]))
    print("  -> files:", nfile)
