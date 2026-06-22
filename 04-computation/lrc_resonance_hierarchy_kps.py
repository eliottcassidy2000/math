"""
The resonance hierarchy & the proof<->disproof dialectic (kps-S31p).
f(t)=min_s ||st||. Its local maxima (WITNESSES) sit at resonances t=a/b; the lonely VALUE
there is ~1/b. M(S) = max witness = 1/(smallest UNKILLED resonance b). Killing resonance b
needs a speed ==0 mod b (a runner parked at origin at t=a/b). The 1/12-core {1..11,13}: b=12
witness (M=1/12); add v=24 (kills b=12) => drops to b=14 (M=1/14). ZETA(-1)=-1/12=-B2/2 is the
lonely-value avatar; ZETA(2) is the floor (resonance density). Functional-equation duals.
"""
from fractions import Fraction as F
def nf(x):
    r=x%1; return min(r,1-r)
def witnesses(S, topn=6):
    S=sorted(set(abs(s) for s in S if s!=0)); C=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for comb in {S[i]+S[j], abs(S[i]-S[j])}:
                if comb:
                    for k in range(1,comb): C.add(F(k,comb))
    vals=sorted(((min(nf(s*t) for s in S), t) for t in C if 0<t<F(1,2)), reverse=True)
    # dedupe by value, keep distinct local peaks
    out=[]; seen=set()
    for v,t in vals:
        if v in seen: continue
        seen.add(v); out.append((v,t))
        if len(out)>=topn: break
    return out
def killed_by(b, S):  # which speeds are ==0 mod b (park at origin at t=a/b)
    return [s for s in S if s%b==0]
print("RESONANCE HIERARCHY -- top witnesses (lonely value, t=a/b), b=denominator:")
for name,S in [("AP {1..13}", list(range(1,14))),
               ("1/12-core {1..11,13}", list(range(1,12))+[13]),
               ("GW {1..11,13,24}", list(range(1,12))+[13,24])]:
    print(f"\n {name}:  M = {witnesses(S,1)[0][0]}")
    for v,t in witnesses(S,5):
        b=t.denominator
        print(f"   witness t={str(t):>6s} (b={b:2d}): lonely value {str(v):>5s}={float(v):.4f};  speeds==0 mod {b}: {killed_by(b,S)}")
print("\n=> M = 1/(smallest b with NO speed ==0 mod b AND a speed at +-1 mod b). The 1/12-core's b=12")
print("   peak is killed by v=24 (24==0 mod 12) => M drops to the b=14 peak (1/14). To go BELOW 1/14")
print("   one must kill b=14 too (a speed ==0 mod 14) AND every smaller resonance -- the COVERING case.")
