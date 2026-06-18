"""
Part 2 of adversarial verification: the LOAD-BEARING partition for k=2.

The k=2 proof partitions covering S3 sets by Vmax:
  - Vmax >= 63: claim C(S) fires via v=Vmax because A = S\{Vmax} = P(11) U {V2}
                has W(A) >= W_min (or the L1 asymptotic bound) so 7*Vmax*W(A) > 1.
  - Vmax <= 62: finite hard core, exhaustively M-checked.

LOAD-BEARING sub-claims to attack:
  (i)  L2: over A = P U {V2} with P = 11-subset of {1..13}, V2 in 14..50,
       and A primitive+'covering-completable', W(A) has NO zero case and W_min=9/3920.
       ATTACK: find any A with W(A)=0 (=> M(A)<=1/14 => dropping Vmax gives no arc).
  (ii) Even if W(A)>0, the partition needs: for Vmax<=62 it's in the core; for
       Vmax>=63, 7*Vmax*W(A) > 1. But W(A) depends on V2 (which is <=Vmax). The
       claim bounds 7*Vmax*W(A) >= 7*63*W_min. We must check that the MINIMUM of
       W(A) over ALL admissible A (any V2 in 14..62) is indeed >= 9/3920, because
       Vmax can be 63 while V2 ranges up to 62.
  (iii) Independently RE-ENUMERATE the k=2 finite hard core (all speeds<=62),
        count it, find worst M, confirm 0 below 1/14.

We must be careful about which A are 'admissible': A = S\{Vmax} where S is a
covering primitive S3 set with k=2. So A = P U {V2}, |P|=11, P subset {1..13},
V2 in 14..Vmax (so V2 <= 61 when Vmax<=62, but V2 can be up to Vmax-... actually
V2 < Vmax). For the Vmax>=63 branch, V2 ranges 14..Vmax-... but the bound used is
W over A; we test ALL V2 in 14..62 to be safe (superset of what's needed for the
>=63 branch where V2 can be anything < Vmax).

Covering condition on S = P U {V2, Vmax}: S must contain a multiple of every q in
2..14. Since the SMALL part lives in {1..13}, q=14 can ONLY be covered by a large
speed (a multiple of 14) OR... no small speed is a multiple of 14. So 14 | V2 or
14 | Vmax. Also note 11,13 multiples: 11|something, 13|something in 2..13 region.
We do NOT pre-impose covering on A alone (A might not be covering); we enumerate
the FULL S and check covering on S. For the W_min sub-claim we enumerate A directly.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, time

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def safe_components(A, h=F(1, 14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u); a = (c - h / u) % 1; b = (c + h / u) % 1
            if a < b: iv.append((a, b))
            else: iv.append((a, F(1))); iv.append((F(0), b))
    iv.sort(); merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else: merged.append((a, b))
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

print("="*78)
print("PART 2A: L2 -- W_min over A = P(11-subset of 1..13) U {V2}, V2 in 14..62.")
print("  ATTACK: hunt for W(A)=0 (zero case) and find true W_min.")
print("="*78)
t0=time.time()
# 11-subsets of {1..13} = choose 2 to remove = C(13,2)=78 parts
parts = [list(c) for c in itertools.combinations(range(1,14), 11)]
print(f"  #11-subsets of 1..13 = {len(parts)}")
zero_cases = []
wmin = None; wmin_arg = None
count = 0
for P in parts:
    for V2 in range(14, 63):
        A = sorted(P + [V2])
        if reduce(gcd, A) != 1:
            continue
        count += 1
        W = Wwidth(A)
        if W == 0:
            zero_cases.append((P, V2))
        if wmin is None or W < wmin:
            wmin = W; wmin_arg = A[:]
print(f"  admissible A tested (primitive): {count}")
print(f"  ZERO cases W(A)=0 found: {len(zero_cases)}")
for P,V2 in zero_cases[:20]:
    print(f"    P={P} V2={V2}")
print(f"  W_min (over V2<=62) = {wmin} = {float(wmin):.6f}  at A={wmin_arg}")
# the claim used V2<=50 giving 9/3920. report both ranges.
wmin50=None; wmin50_arg=None
for P in parts:
    for V2 in range(14,51):
        A=sorted(P+[V2])
        if reduce(gcd,A)!=1: continue
        W=Wwidth(A)
        if wmin50 is None or W<wmin50: wmin50=W; wmin50_arg=A[:]
print(f"  W_min (V2<=50, claim's range) = {wmin50} = {float(wmin50):.6f} at {wmin50_arg}")
print(f"  claim 9/3920 = {float(F(9,3920)):.6f}; matches V2<=50? {wmin50==F(9,3920)}")
print(f"  ({time.time()-t0:.0f}s)")

# CRITICAL: the >=63 branch. There Vmax>=63 and V2 can be up to Vmax-1 (>=62 even).
# The bound used is 7*Vmax*W(A) >= 7*63*W_min(over all V2). We need W_min over ALL
# V2 that can appear in a Vmax>=63 set, i.e. V2 in 14..Vmax-1 -> unbounded V2!
# BUT for V2 large, L1 (asymptotic 6) takes over. The claim says: for V2<=50 use L2,
# for V2>=51 use L1. So the partition is by V2, not just Vmax. Let me re-read:
# "every k=2 covering S3 set with Vmax>=63 satisfies C via Vmax" -- combining L1+L2.
# So we must verify: for EVERY admissible A=P U {V2} (any V2>=14) with Vmax>=63
#   (Vmax>V2), 7*Vmax*W(A) > 1. We check the WORST over the SECOND speed:
#   min over A of 7*63*W(A) for V2<=62, AND 7*Vmax*(L1 bound) for V2>=51.
print()
print("  BRANCH CHECK: for Vmax=63 (smallest >=63), need 7*63*W(A) > 1 for all")
print("    admissible A=P U {V2}, V2 in 14..62. i.e. W(A) > 1/441.")
need = F(1,441)
bad=[]
for P in parts:
    for V2 in range(14,63):
        A=sorted(P+[V2])
        if reduce(gcd,A)!=1: continue
        W=Wwidth(A)
        if W <= need:
            bad.append((P,V2,W))
print(f"    1/441 = {float(need):.6f};  W_min found = {float(wmin):.6f}")
print(f"    A with W(A) <= 1/441 (would BREAK the Vmax=63 branch): {len(bad)}")
for P,V2,W in bad[:20]:
    print(f"      P={P} V2={V2} W={W}={float(W):.6f}  7*63*W={float(7*63*W):.4f}")
print(f"  ({time.time()-t0:.0f}s)")
print("DONE PART 2A.")
