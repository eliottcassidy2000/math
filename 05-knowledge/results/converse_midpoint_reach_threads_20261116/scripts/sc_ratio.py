#!/usr/bin/env python3
"""sc_ratio.py: self-converse tournament classes SC(n) for n <= 9 (all classes and strong classes) by canonical
labelling (nauty labelg on digraph6) of each class and of its converse; the ratio SC/(converse pairs)."""
import subprocess, sys
def d6_decode(s):
    s = s.strip(); assert s[0] == '&'; n = ord(s[1]) - 63; bits = ''.join(format(ord(c) - 63, '06b') for c in s[2:])
    return n, [[int(bits[i*n+j]) for j in range(n)] for i in range(n)]
def d6_encode(n, M):
    bits = ''.join(str(M[i][j]) for i in range(n) for j in range(n)); bits += '0' * ((-len(bits)) % 6)
    return '&' + chr(n + 63) + ''.join(chr(int(bits[i:i+6], 2) + 63) for i in range(0, len(bits), 6))
def canon(lines):
    out = subprocess.run(['labelg', '-q'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split()
    assert len(out) == len(lines), (len(out), len(lines)); return out
print("n  classes  SC  pairs  SC/pairs | strong  SC_strong  pairs_strong  ratio_strong")
for n in range(3, 10):
    res = []
    for flag in ([], ['-c']):
        raw = subprocess.run(['gentourng', '-q', '-z'] + flag + [str(n)], capture_output=True, text=True).stdout.split()
        rev = []
        for s in raw:
            m, M = d6_decode(s); rev.append(d6_encode(m, [[M[j][i] for j in range(m)] for i in range(m)]))
        c1 = canon(raw); c2 = canon(rev); sc = sum(1 for a, b in zip(c1, c2) if a == b); tot = len(raw); pairs = (tot - sc) // 2
        res.append((tot, sc, pairs, (sc / pairs) if pairs else float('inf')))
    (t, s, p, r), (ts, ss, ps, rs) = res
    print(f"{n}  {t}  {s}  {p}  {r:.4f} | {ts}  {ss}  {ps}  {rs:.4f}")
