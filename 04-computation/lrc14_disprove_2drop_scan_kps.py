"""
DISPROVE-side: 2-drop scan for LRC(14). Hunt L < 1/1260 and growing-lcm L->0.

Drop e1,e2 from AP {1..13}, add w1,w2 (both NOT in remaining core, distinct).
Use a float prescreen (SAFE lower bound is not needed here; we want SMALL L, and
float L is accurate to ~1e-12 for these), then exact-confirm any candidate < 0.01.

Focus windows informed by single-drop data:
  - drop-6 minimizes at w=69 (large); drop-10 at 20; drop-12 at 36.
  - So scan added speeds up to a decent window, and specifically probe the
    "coordinated multiple" structure (w1=c*e1', w2=c*e2' resonances).
"""
from fractions import Fraction as Fr
from math import gcd, lcm
from functools import reduce
from itertools import combinations

def danger_arcs(v):
    w = Fr(1, 14*v); A = []
    for k in range(v+1):
        c = Fr(k, v); lo = c - w; hi = c + w
        if lo < 0:
            A += [(Fr(0), hi), (1+lo, Fr(1))]
        elif hi > 1:
            A += [(lo, Fr(1)), (Fr(0), hi-1)]
        else:
            A.append((lo, hi))
    return A

_cache = {}
def darcs(v):
    r = _cache.get(v)
    if r is None:
        r = danger_arcs(v); _cache[v] = r
    return r

def L_exact(S):
    A = []
    for v in S:
        A.extend(darcs(v))
    A.sort(key=lambda t: (t[0], t[1]))
    tot = Fr(0); cl = ch = None
    for a, b in A:
        if b <= a: continue
        if ch is None:
            cl, ch = a, b
        elif a <= ch:
            if b > ch: ch = b
        else:
            tot += ch - cl; cl, ch = a, b
    if ch is not None:
        tot += ch - cl
    return 1 - tot

# float danger arcs for fast prescreen
import bisect
def Lf(S):
    A = []
    for v in S:
        wv = 1.0/(14*v)
        for k in range(v+1):
            c = k/v; lo = c-wv; hi = c+wv
            if lo < 0:
                A.append((0.0,hi)); A.append((1+lo,1.0))
            elif hi > 1:
                A.append((lo,1.0)); A.append((0.0,hi-1))
            else:
                A.append((lo,hi))
    A.sort()
    tot=0.0; cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch:
            if b>ch: ch=b
        else:
            tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot

THRESH = Fr(1,1260)
THRESHf = 1.0/1260
AP = list(range(1,14))

print("=== (B) 2-drop scan, added speeds in [14, Wmax] ===")
Wmax = 80
# all 2-drops
results = []
champ = None  # (L, S)
count = 0
for e1, e2 in combinations(range(1,14), 2):
    core = [x for x in AP if x not in (e1,e2)]
    cand = [w for w in range(2, Wmax+1) if w not in core]
    for i in range(len(cand)):
        w1 = cand[i]
        for j in range(i+1, len(cand)):
            w2 = cand[j]
            S = sorted(core + [w1, w2])
            count += 1
            lf = Lf(S)
            if lf < THRESHf + 1e-9:  # promising
                L = L_exact(S)
                if L == 0: continue
                if L < THRESH:
                    print(f"  *** BELOW 1/1260: drop({e1},{e2}) add({w1},{w2}) L={L}={float(L):.8f} lcm={lcm(*S)}")
                if champ is None or L < champ[0]:
                    champ = (L, S, (e1,e2), (w1,w2))
print(f"  scanned {count} configs")
if champ:
    L,S,ee,ww = champ
    print(f"  2-drop champion: L={L}={float(L):.8f} below 1/1260? {L<THRESH}")
    print(f"    S={S} drop{ee} add{ww} lcm={lcm(*S)}")
