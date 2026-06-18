"""
Independent verification of the 4-attack synthesis for LRC(14) S3.
kind-pasteur-S4-wf verifier pass.
"""
from fractions import Fraction as F
from math import gcd
import itertools

def gcd_list(L):
    g = 0
    for x in L: g = gcd(g, x)
    return g

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
        if merged and a <= merged[-1][1]: merged[-1] = (merged[-1][0], max(merged[-1][1], b))
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
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1: ws.append(sc[0][1] + (1 - sc[-1][0]))
    return max(ws)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2): C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x * t) for x in S)
        if v > b: b = v
    return b

def is_cov(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

R = []
def out(s): R.append(s); print(s)

out("="*72)
out("(A) ALL-MULT7-LARGE: does the critic's window-block actually arise IN SCOPE?")
out("="*72)
# Window-collapse reduces, on |u|<=1/(2 V*), to {w_i}-safety where 7 w_i are the
# multiples of 7 in S. ALL-MULT7-LARGE := every multiple of 7 in S exceeds V*.
# Covering REQUIRES a multiple of 14 (=2*7), which is even & a multiple of 7.
# So in ALL-MULT7-LARGE the multiple-of-14 runner is itself > V*, i.e. one of the
# large mult-of-7 runners is EVEN.  The critic's V*=19 example used mults {21,273}
# (w_i={3,39}), both odd => NO multiple of 14 => NOT covering.  Test that block:
ws = [3, 39]; Vstar = 19; bound = F(1, 2*Vstar)
sc = safe_components(ws)
hits = [(max(a,F(0)),min(b,bound)) for (a,b) in sc if min(b,bound)>max(a,F(0))]
out(f"  {{3,39}}-safe in (0,{bound}]: {hits}  (window blocked={not hits})")
out("  But {21,273} both ODD => no multiple of 14 => any S with ONLY these mult-of-7 runners is NOT covering.")

out("")
out("  Brute force: REAL covering primitive S3 ALL-MULT7-LARGE sets with V*=19.")
small_pool = [v for v in range(1,20) if v % 7 != 0]   # non-mult-of-7, <=19
found = []
# large mult-of-7 runners: need an EVEN one for q=14; include various combos
even7 = [28,42,56,70,84,98,112,126,140,154,168,182,196]
odd7  = [21,35,49,63,77,91,105,119,133,147,161,175,189,273]
# Build large = subset of mult-of-7, all >19, must contain >=1 even (for 14) and
# a multiple of 13 somewhere in S (could be 273=3*7*13, or small 13? 13 is in pool).
# Keep it tractable: large size 2 or 3.
large_candidates = []
for e in even7:
    large_candidates.append([e])           # single even mult of 7
    for o in odd7:
        if o != e: large_candidates.append(sorted([e,o]))
seen = set()
for large in large_candidates:
    if any(x <= 19 for x in large): continue
    lt = tuple(large)
    if lt in seen: continue
    seen.add(lt)
    need = 13 - len(large)
    if need < 1: continue
    pool = [v for v in small_pool]  # 19 must appear to pin V*=19
    if 19 not in pool: continue
    rest = [v for v in pool if v != 19]
    if need-1 > len(rest): continue
    for combo in itertools.combinations(rest, need-1):
        S = sorted([19] + list(combo) + list(large))
        if len(S) != 13: continue
        if gcd_list(S) != 1: continue
        if not is_cov(S): continue
        nonmult = [v for v in S if v % 7 != 0]
        if max(nonmult) != 19: continue
        if any(v % 7 == 0 and v <= 19 for v in S): continue  # ALL-MULT7-LARGE
        kbig = sum(1 for v in S if v > 13)
        if kbig < 2: continue  # S3
        m = Mval(S)
        found.append((m, tuple(S)))
        if len(found) > 4000: break
    if len(found) > 4000: break

found = sorted(set(found))
out(f"  found {len(found)} in-scope ALL-MULT7-LARGE V*=19 covering primitive S3 sets")
if found:
    mn = min(found)
    out(f"  MIN M over them = {mn[0]} = {float(mn[0]):.5f} at S={list(mn[1])}")
    out(f"  min M >= 1/14 ? {mn[0] >= F(1,14)}  (1/14={float(F(1,14)):.5f})")
    # Does the window-collapse witness exist for the min set? reduced w_i:
    Smin = list(mn[1]); m7 = [v//7 for v in Smin if v % 7 == 0]
    Vs = max(v for v in Smin if v % 7 != 0)
    bnd = F(1,2*Vs); scw = safe_components(m7)
    wh = [(max(a,F(0)),min(b,bnd)) for (a,b) in scw if min(b,bnd)>max(a,F(0))]
    out(f"  window-collapse for min set: reduced w_i={sorted(set(m7))}, u-window (0,{bnd}], in-window safe pts: {wh}")
    out(f"    => window-collapse {'GIVES' if wh else 'FAILS to give'} a global witness here")

out("")
out("="*72)
out("(B) Window-collapse sharpness equality: 1/7 - V*/(14 V*) = 1/14")
out("="*72)
out(f"  1/7 - 1/14 = {F(1,7)-F(1,14)} ; equals 1/14 ? {F(1,7)-F(1,14)==F(1,14)}")

out("")
out("="*72)
out("(C) 7-adic residue rule: ||v*k/7|| values for gcd(k,7)=1")
out("="*72)
bad = False
for k in [1,2,3,4,5,6]:
    for v in range(1,15):
        val = nrm(F(v*k,7))
        if v % 7 == 0:
            if val != 0: bad = True
        else:
            if val not in (F(1,7),F(2,7),F(3,7)): bad = True
out(f"  rule holds (0 iff 7|v, else in 1/7,2/7,3/7): {not bad}; all nonzero >= 1/7 > 1/14: {F(1,7)>F(1,14)}")

out("")
out("="*72)
out("(D) AP family two witnesses for S={1..12,m}, m=182k")
out("="*72)
def witness1(m):  # tau=2/27
    t = F(2,27)
    return min(nrm(F(v)*t) for v in list(range(1,13))+[m])
def witness2(m):  # tau = (m/13)/(m+1)
    assert m % 13 == 0
    t = F(m//13, m+1)
    return min(nrm(F(v)*t) for v in list(range(1,13))+[m])
w1_fail = []
for m in range(182, 182*60+1, 182):
    lv = witness1(m)
    if lv < F(1,14):
        w1_fail.append((m, lv, witness2(m)))
out(f"  witness1 (2/27) FAILS (<1/14) for m=182k: {[(m, str(lv), str(lv2)) for m,lv,lv2 in w1_fail][:10]} ...")
out(f"  for those m, witness2 level >= 1/14 ? {all(lv2>=F(1,14) for _,_,lv2 in w1_fail)}")
out(f"  witness2 min over the witness1-failures: {min((lv2 for _,_,lv2 in w1_fail), default='n/a')}")
# combined floor across all k<60:
floor = min(max(witness1(m), witness2(m)) for m in range(182,182*60+1,182))
out(f"  combined best-witness floor over k=1..59: {floor} = {float(floor):.5f} ; >=1/14: {floor>=F(1,14)}")
# cross-check with true Mval at small k
out(f"  true M({182})={Mval(list(range(1,13))+[182])} ; M({364})={Mval(list(range(1,13))+[364])}")
out(f"  covering test {{1..12,m}}: 182|m iff covering? check m in 14..400:")
mismatch = [m for m in range(14,401) if (is_cov(list(range(1,13))+[m])) != (m % 182 == 0)]
out(f"    mismatches: {mismatch}  (empty = confirmed)")

out("")
out("="*72)
out("(E) Hard case S* and realized floor descent (exact M)")
out("="*72)
Sstar = [1,2,3,5,7,8,9,10,11,12,13,38,42]
out(f"  S*={Sstar}: covering={is_cov(Sstar)}, primitive={gcd_list(Sstar)==1}, M={Mval(Sstar)} = {float(Mval(Sstar)):.5f}")
# floor descent claimed values
descent = {
 "{1,84} pair set":  [1,2,3,4,5,6,7,9,10,11,12,13,84],
 "{1,182} pair set": [1,2,3,4,5,6,7,8,9,10,11,12,182],
}
for name,S in descent.items():
    if len(set(S))==13:
        out(f"  {name}: S={S} cov={is_cov(S)} prim={gcd_list(S)==1} M={Mval(S)} = {float(Mval(S)):.5f}")
out(f"  14/183 = {float(F(14,183)):.5f}; > 1/14 ? {F(14,183)>F(1,14)} ; margin={F(14,183)-F(1,14)}")
out(f"  2/23   = {float(F(2,23)):.5f}; > 1/14 ? {F(2,23)>F(1,14)}")

out("")
out("="*72)
out("(F) k=3 board worst set")
out("="*72)
Sk3 = [1,2,3,4,5,7,8,11,12,13,18,20,28]
out(f"  S={Sk3}: len={len(set(Sk3))} cov={is_cov(Sk3)} prim={gcd_list(Sk3)==1} M={Mval(Sk3)} = {float(Mval(Sk3)):.5f}")
out(f"  k=#{{v>13}}={sum(1 for v in Sk3 if v>13)} ; Vmax={max(Sk3)} >= 13*Vmin={13*min(Sk3)} ? {max(Sk3)>=13*min(Sk3)}")

out("")
out("="*72)
out("(G) fixed-small-part in-scope exact M (P={1,2,3}, offsets 0..9, single cluster)")
out("="*72)
for V0 in [14,20,30,50,81]:
    L = [V0+d for d in range(10)]  # offsets 0..9 => 10 large runners
    P = [1,2,3]
    S = sorted(set(P)+L) if False else sorted(set(P+L))
    if len(S)!=13:
        out(f"  V0={V0}: |S|={len(S)} (skip, need 13)"); continue
    out(f"  V0={V0}: S={S} cov={is_cov(S)} prim={gcd_list(S)==1} M={Mval(S)} = {float(Mval(S)):.5f} >=1/14:{Mval(S)>=F(1,14)}")

with open(r"C:\Users\Eliott\Documents\GitHub\math\05-knowledge\results\lrc14_verify_synthesis_kps-S4-wf.out","w") as f:
    f.write("\n".join(R))
out("\n[written results file]")
