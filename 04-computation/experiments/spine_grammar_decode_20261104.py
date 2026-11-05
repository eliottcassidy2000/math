#!/usr/bin/env python3
"""spine_grammar_decode_20261104.py -- opus-2026-11-04-S79

Decode the witness words printed by spine_grammar_census_20261104.cpp (letter indices per chain vertex, indices into
lettersFor(c) = all (T, b, A) with T over the c pending slots, b new slots, A apex subset, in that nested order) into
explicit tournaments: chain length m, reversed arcs (tail > head), apex incidences; and summarise the witness
structure (number of reversed arcs, chain length, cutwidth profile) for every recorded value.
Usage: python3 spine_grammar_decode_20261104.py k w census_file [limit]
"""
import sys, re

def letters_for(c, w, k):
    L = []
    for T in range(1 << c):
        t = bin(T).count('1')
        for b in range(0, w - (c - t) + 1):
            for A in range(1 << k):
                L.append((T, b, A))
    return L

def decode(word, w, k):
    pending = []; arcs = []; apex = []
    for v, li in enumerate(word):
        c = len(pending); L = letters_for(c, w, k); T, b, A = L[li]
        newp = []
        for i in range(c):
            if T >> i & 1: arcs.append((v, pending[i]))
            else: newp.append(pending[i])
        for j in range(b): newp.append(v)
        pending = newp; apex.append(A)
    assert not pending
    return len(word), arcs, apex

def main():
    k, w, path = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    limit = int(sys.argv[4]) if len(sys.argv) > 4 else 10**9
    rows = []
    for line in open(path):
        m = re.match(r'\s+(\d+) = ([\d,]*)', line)
        if not m: continue
        N = int(m.group(1))
        if N > limit: continue
        word = [int(x) for x in m.group(2).split(',') if x != '']
        if not word: rows.append((N, 0, [], [])); continue
        mlen, arcs, apex = decode(word, w, k)
        rows.append((N, mlen, arcs, apex))
    print(f"{len(rows)} witnesses decoded from {path}")
    for (N, mlen, arcs, apex) in rows:
        spans = sorted(t - h for (t, h) in arcs)
        ap = ''.join(str(a) for a in apex) if k else ''
        print(f"  {N}: m={mlen} arcs={len(arcs)} spans={spans} arcs={arcs}" + (f" apex={ap}" if k else ''))
    import collections
    byarcs = collections.Counter(len(arcs) for (_, _, arcs, _) in rows)
    print("values by number of reversed arcs:", sorted(byarcs.items()))
    bym = collections.Counter(mlen for (_, mlen, _, _) in rows)
    print("values by chain length:", sorted(bym.items()))

main()
