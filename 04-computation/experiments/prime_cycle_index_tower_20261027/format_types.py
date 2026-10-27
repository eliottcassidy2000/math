"""Markdown tables of the prime cycle-index types t{n}.txt written by tower.py, and the non-rigid asymptotics.
Usage: python3 format_types.py n [n ...]   (reads out/t{n}.txt)"""
import sys, re
from math import factorial, comb

def parse_ms(s):
    res = []
    for tok in s.split():
        ct, mult = tok.rsplit(':', 1); res.append((ct.strip(','), int(mult)))
    return res

def fmt(ms):
    # identity first, then by decreasing number of fixed points / cycle structure
    def key(x):
        parts = [tuple(map(int, p.split('^'))) for p in x[0].split(',') if p]
        fixed = sum(e for l, e in parts if l == 1)
        return (-fixed, parts)
    return ", ".join(f"`{ct.replace(',', ' ')}` ({m})" for ct, m in sorted(ms, key=key))

for n in map(int, sys.argv[1:]):
    rows = []
    for line in open(f"out/t{n}.txt"):
        m = re.match(r'TYPE count=(\d+) k=(\d+) aut=(\d+) AUT:(.*?) SC=(\d) ANTI:(.*)$', line.strip())
        cnt, G, A, sc, anti = int(m.group(1)), int(m.group(3)), parse_ms(m.group(4)), int(m.group(5)), parse_ms(m.group(6))
        rows.append((G, sc, A, anti, cnt))
    rows.sort(key=lambda r: (r[0], fmt(r[2]), r[1], fmt(r[3])))
    print(f"\n### order {n}: {len(rows)} types, {sum(r[4] for r in rows):,} prime classes\n")
    print("| \\|Aut\\| | Aut cycle types (elements) | self-converse | anti-automorphism cycle types (elements) | classes |")
    print("|---|---|---|---|---|")
    for G, sc, A, anti, cnt in rows:
        print(f"| {G} | {fmt(A)} | {'yes' if sc else 'no'} | {fmt(anti) if sc else '–'} | {cnt:,} |")
    nr = sum(r[4] for r in rows if r[0] > 1); nrsc = sum(r[4] for r in rows if r[0] > 1 and r[1])
    hist = {}
    for G, sc, A, anti, cnt in rows: hist[G] = hist.get(G, 0) + cnt
    print(f"\nnon-rigid {nr:,} (self-converse {nrsc:,}); |Aut| histogram {dict(sorted(hist.items()))}")
    f = n - 6
    if f >= 0:
        lead = 2 ** (comb(f, 2) + 2 * f + 5) / (12 * factorial(f))
        print(f"leading term 2^(C(n-6,2)+2(n-6)+5)/(12 (n-6)!) = {lead:,.0f} (times the prime fraction of the invariant tournaments)")
