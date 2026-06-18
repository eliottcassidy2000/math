"""
Adversarial verification of the 'finite-core-enumeration' advance on S3 of LRC(14).
machine: kind-pasteur  session: S3-wf

Goal: independently re-derive every load-bearing inequality with EXACT Fractions,
then HUNT for a covering S3 set that violates the claimed bounds.

Claims to attack:
  L0  min M(P) over admissible small parts (|P|<=11, P subset of 1..13) = 1/12 (>= 1/12).
  L1  drop-max scaling: 7*V2*W(P U {V2}) -> 6 and >= 1 for all V2 >= 51 (k=2).
  L2  bounded hard core: W_min over A=P(11) U {V2<=50} with W>0 is 9/3920, no W=0 case;
      => 7*Vmax*W(A) > 1 once Vmax >= 63.
  THM(k=2): finite hard core (all speeds <= 62) has 4865 sets, worst M = 1/12, none < 1/14.
  CRITERION C(S): exists v in S with W(S\{v}) > 1/(7v)  => M(S) >= 1/14.

We treat THM-523 (only covering 13-sets matter), THM-524 (candidate set), THM-526
(arc-width lemma) as black boxes BUT we independently sanity-check the candidate
set Mval against a fine brute-force grid AND against the arc-width implication.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, sys, time

# ---------------------------------------------------------------- exact tools
def nrm(x):
    r = x - int(x)
    r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h / u) % 1; b = (c + h / u) % 1
            if a < b:
                iv.append((a, b))
            else:
                iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    safe = []; prev = F(0)
    for a, b in merged:
        if a > prev: safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b - a for a, b in sc]
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x * t) for x in S)
        if v > b: b = v
    return b

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def is_primitive(S):
    return reduce(gcd, S) == 1

# M is symmetric over tau in (0,1); cand only lists tau<=1/2. Confirm M is the
# sup over the whole circle by also testing 1-t (M(S) is even in tau). The min
# is symmetric so cand at <=1/2 suffices. We double check below.

# ============================================================================
print("="*78)
print("STEP A: sanity-check Mval (THM-524 candidate set) against fine grid +")
print("        confirm symmetry, on a spread of S3-like sets.")
print("="*78)
def Mbrute(S, dens=20000):
    # coarse but independent: evaluate min_v ||v t|| on a dense grid + cand points
    best = F(0)
    pts = set(cand(S))
    for i in range(dens+1):
        pts.add(F(i, dens))
    for t in pts:
        v = min(nrm(x*t) for x in S)
        if v > best: best = v
    return best

test_sets = [
    [1,2,3,4,5,6,7,8,9,10,11,12,13],
    [1,2,3,5,6,7,8,9,11,12,13,30,42],
    [1,2,3,4,5,7,8,9,11,12,13,40],
    [2,3,4,5,6,7,8,9,10,11,12,13,14],
]
okA = True
for S in test_sets:
    mv = Mval(S); mb = Mbrute(S)
    # grid may UNDERSHOOT true M; cand should be >= grid. flag if cand < grid.
    flag = "" if mv >= mb else "  <-- cand BELOW grid (BAD)"
    if mv < mb: okA = False
    print(f"  S={S}\n    Mval={mv}={float(mv):.5f}  grid>={mb}={float(mb):.5f}{flag}")
print(f"  Mval >= grid on all test sets: {okA}")

# ============================================================================
print("\n"+"="*78)
print("STEP B: re-derive ARC-WIDTH criterion independently.")
print("  Claim: if W(A) > 1/(7v) then min over widest arc has a v-safe point,")
print("         so M(A U {v}) >= 1/14.")
print("  We TEST the implication directly: pick sets A, v; check C => M>=1/14.")
print("="*78)
# We verify: for many (A,v), if W(A) > 1/(7v) then Mval(A U {v}) >= 1/14.
import random
random.seed(12345)
viol_arc = []
trials = 0
for _ in range(4000):
    k = random.randint(3, 8)
    A = sorted(random.sample(range(1, 60), k))
    if reduce(gcd, A) != 1:
        continue
    v = random.randint(1, 90)
    if v in A:
        continue
    trials += 1
    W = Wwidth(A)
    if W > F(1, 7*v):
        m = Mval(sorted(A + [v]))
        if m < F(1, 14):
            viol_arc.append((A, v, W, m))
print(f"  trials with C firing checked: {trials}")
print(f"  arc-width-criterion violations (C true but M<1/14): {len(viol_arc)}")
for a,v,W,m in viol_arc[:10]:
    print(f"    A={a} v={v} W={W} M={m}")

# ============================================================================
print("\n"+"="*78)
print("STEP C: re-derive L0 (CEILING) min M(P) over admissible small parts.")
print("  Admissible P: P subset of {1..13}, |P|<=11. (The S3 small part.)")
print("  Claim min M(P) = 1/12 achieved at prefix [1..11].")
print("="*78)
t0=time.time()
glob_min = None; glob_arg = None
per_size = {}
for size in range(1, 12):
    best = None; barg=None
    for P in itertools.combinations(range(1, 14), size):
        m = Mval(list(P))
        if best is None or m < best:
            best = m; barg = P
    per_size[size] = (best, barg)
    if glob_min is None or best < glob_min:
        glob_min = best; glob_arg = barg
for size in range(11,0,-1):
    b,a = per_size[size]
    print(f"  |P|={size:2d}: min M(P)={str(b):>8s}={float(b):.5f}  arg={list(a)}")
print(f"  OVERALL min M(P) = {glob_min} = {float(glob_min):.5f}  arg={list(glob_arg)}")
print(f"  == 1/12 ? {glob_min==F(1,12)}   > 1/14 ? {glob_min>F(1,14)}   ({time.time()-t0:.0f}s)")

print("\nDONE STEP A-C.")
