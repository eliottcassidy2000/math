#!/usr/bin/env python3
"""
Dormancy census for the ATLAS OF ATLASES (mac-mini-2026-07-20-S124)
===================================================================
DEFINITION: an idea is DORMANT if its identifier appears nowhere outside the log
that raised it (raised once, never picked up). Items newer than CUT are excluded
as FRESH -- new is not abandoned (this exclusion matters: today's T1551-T1555
would otherwise read as a graveyard).

RESULTS (2026-07-20): tangents 50 dormant / 572 parsed entries (3.6%);
hypotheses 0 dormant-and-open of 2707 (all 46 index-only entries are CLOSED).
The ~28-tangent "abandoned Lean sprint" cluster was then checked directly and
found to be LIVE CODE (7/8 modules imported in the TournamentH7 root) with stale
tangent labels -- so genuinely abandoned mathematical content is ~nil.

BLIND SPOT (stated in the atlas): this detects dropped THREADS, not stalled
PROBLEMS. An idea cited often but never advanced looks healthy here.
"""
import re, os, collections, datetime

ROOTS = ['00-navigation','01-canon','05-knowledge','07-reflections','03-artifacts','02-court','04-computation']
CUT = datetime.date(2026, 7, 17)

def citation_index(pat):
    where = collections.defaultdict(set)
    for r in ROOTS:
        for dp, dn, fn in os.walk(r):
            for f in fn:
                if not (f.endswith('.md') or f.endswith('.py')): continue
                p = os.path.join(dp, f)
                try: txt = open(p, encoding='utf-8', errors='ignore').read()
                except Exception: continue
                for m in set(pat.findall(txt)): where[m].add(p)
    return where

Tw = citation_index(re.compile(r'\bT(\d{3,4})\b'))
tl = open('00-navigation/TANGENTS.md', encoding='utf-8', errors='ignore').read()
ents = re.findall(r'\*\*T(\d{3,4})\*\*\s*\[([^\]]*)\]\s*([^\n]{0,170})', tl)
dorm = []
for num, src, rest in ents:
    others = {f for f in Tw.get(num, set()) if 'TANGENTS.md' not in f and 'HISTORIAN' not in f}
    if others: continue
    d = re.search(r'(\d{4})-(\d{2})-(\d{2})', src)
    dt = datetime.date(*map(int, d.groups())) if d else None
    if dt and dt >= CUT: continue           # fresh, not abandoned
    dorm.append((int(num), src.split(';')[0].strip(), rest.strip()))
nmark = len(re.findall(r'\*\*T\d{3,4}\*\*', tl))
print(f"tangent markers: {nmark}, parsed entries: {len(ents)}")
print(f"DORMANT tangents (pre-{CUT}, never cited outside the log): {len(dorm)}")
print("  clusters:", collections.Counter(s for _, s, _ in dorm).most_common(8))

Hw = citation_index(re.compile(r'\bHYP-(\d{3,5})\b'))
idx = open('05-knowledge/hypotheses/INDEX.md', encoding='utf-8', errors='ignore').read()
stat = {}
for num, st in re.findall(r'\*\*HYP-(\d{3,5})[^*]*\*\*[:\s]*([A-Z\-]+)?', idx):
    stat.setdefault(num, (st or '').strip())
oo = oc = 0
for num, st in stat.items():
    if {f for f in Hw.get(num, set()) if 'hypotheses/INDEX' not in f}: continue
    if st.startswith('OPEN'): oo += 1
    else: oc += 1
print(f"hypothesis entries: {len(stat)}; index-only: OPEN={oo} CLOSED={oc}")
print("  => an index-only hypothesis with a CLOSED status is finished, not abandoned.")
