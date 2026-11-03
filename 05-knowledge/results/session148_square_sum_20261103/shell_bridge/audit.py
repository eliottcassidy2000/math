from itertools import combinations
from math import isqrt
A={10,6,3,1,8,17};B={12,4,5,11};C={2,7,9,18,16}
X={2,3,5,7,11,13};Y={17,19,23,31,47};Z={29,41,59,71}
P=A|B|C; Q=X|Y|Z; f={**dict(zip(sorted(A),sorted(X))),**dict(zip(sorted(B),sorted(Z))),**dict(zip(sorted(C),sorted(Y)))}
def edge(a,b): return a!=b and isqrt(a+b)**2==a+b
src={tuple(sorted((f[a],f[b]))) for a,b in combinations(sorted(P),2) if edge(a,b)}
tgt={(a,b) for a,b in combinations(sorted(Q),2) if edge(a,b)}
if (len(src),len(tgt),len(src&tgt),len(src-tgt),len(tgt-src))!=(12,14,5,7,9):raise RuntimeError('shell edge mismatch')
if edge(2,3):raise RuntimeError('hostile invalid')
print('PASS: source12 target14 preserved5 lost7 created9; no graph isomorphism by edge count')
