#!/usr/bin/env python3
"""constants_atlas_20260926_context.py -- for chosen constants, print the sentence in which
each recurring value occurs, one or two per thread, from the census JSON written by
constants_atlas_20260926_mine.py.  Used to type each recurrence as STRUCTURAL / TRANSPORTED /
NUMEROLOGY by reading the actual mechanism in each thread.
Usage: python3 constants_atlas_20260926_context.py <repo> key1 key2 ...   (keys like i:139, d:0.6309, f:11/8)
"""
import json, os, re, sys

root = sys.argv[1]
keys = sys.argv[2:]
J = json.load(open(os.path.join(root, '05-knowledge/results/constants_atlas_20260926_mine.json')))

def token_of(key):
    kind, val = key.split(':', 1)
    if kind == 'i':
        return re.compile(r'(?<![\w.^])' + re.escape(val) + r'(?![\w.])')
    if kind == 'f':
        return re.compile(r'(?<![\w.^])' + re.escape(val) + r'(?![\w.])')
    return re.compile(re.escape(val[:5]))

for key in keys:
    if key not in J:
        print("== %s: not in census (needs >= 4 threads)" % key)
        continue
    rx = token_of(key)
    print("== %s ==" % key)
    for thread, files in sorted(J[key].items(), key=lambda kv: -len(kv[1])):
        shown = 0
        for f in files:
            try:
                lines = open(os.path.join(root, f), encoding='utf-8', errors='ignore').read().split('\n')
            except Exception:
                continue
            for ln in lines:
                if rx.search(ln):
                    s = ln.strip()
                    i = rx.search(s).start()
                    print("   [%s] %s: ...%s..." % (thread, f.split('/')[-1][:55], s[max(0, i - 90):i + 110].replace('\n', ' ')))
                    shown += 1
                    break
            if shown >= 2:
                break
