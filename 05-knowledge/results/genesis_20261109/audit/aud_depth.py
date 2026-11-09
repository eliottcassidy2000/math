"""AUDIT G1: creator depth of the 193 order-12 escaper classes (own code: aud_core subs + census spectra).
depth 2 <=> some 10-subtournament U has H(U) not in S(9); depth 3 <=> not depth 2 and some 9-subtournament has
H not in S(8)."""
import os, subprocess
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
S = {n: set(int(x) for x in open(os.path.join(HERE, 'S_%d.txt' % n))) for n in range(1, 12)}
rows = [l.split() for l in open(os.path.join(HERE, 'classes193_props.txt'))]
codes = [r[0] for r in rows]
Hs = [int(r[1][2:]) for r in rows]
def subs(K, cs):
    inp = ''.join('12 %s\n' % c for c in cs)
    r = subprocess.run([os.path.join(HERE, 'aud_core'), 'subs', str(K)], input=inp, capture_output=True, text=True, check=True).stdout.splitlines()
    return [list(map(int, l.split()[1:])) for l in r]
s2 = subs(2, codes)
assert all(len(x) == 66 for x in s2)
fc10 = [sum(1 for h in x if h not in S[9]) for x in s2]
print('classes with a first-creator 10-subtournament (depth 2):', sum(1 for f in fc10 if f > 0), ' without:', sum(1 for f in fc10 if f == 0))
print('FC 10-subtournaments per depth-2 class: min %d max %d' % (min(f for f in fc10 if f), max(fc10)))
rest = [i for i, f in enumerate(fc10) if f == 0]
s3 = subs(3, [codes[i] for i in rest])
for i, x in zip(rest, s3):
    assert len(x) == 220
    k = sum(1 for h in x if h not in S[8])
    print('  depth>2 class H=%d: FC 9-subtournaments %d of 220 -> depth %s; code %s' % (Hs[i], k, 3 if k else '>3', codes[i]))
