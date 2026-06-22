"""
RIGOROUS re-verification (kps-S31n): is {1..11,13,12j} really tight (M=1/14) for j=2,3,4,5?
M = max_t min_s ||st||. Exact critical points: t = k/(s_i +- s_j) (where ||s_i t||=||s_j t||)
for all pairs, plus k/(2 s_i) (peaks). The TRUE max-min is among these rationals.
"""
from fractions import Fraction as F
def nf(x):
    r=x%1; return min(r,1-r)
def M_rigorous(S):
    S=sorted(set(abs(s) for s in S if s!=0))
    cands=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for comb in ({S[i]+S[j], abs(S[i]-S[j])}):
                if comb==0: continue
                for k in range(1, comb): cands.add(F(k, comb))
        for k in range(1, 2*S[i]): cands.add(F(k, 2*S[i]))   # peaks of ||s t||
    best=F(0); arg=None
    for t in cands:
        if not (0<t<1): continue
        mn=min(nf(s*t) for s in S)
        if mn>best: best=mn; arg=t
    return best, arg
tests = {
  "AP {1..13}": list(range(1,14)),
  "GW {1..11,13,24}": list(range(1,12))+[13,24],
  "{1..11,13,36}": list(range(1,12))+[13,36],
  "{1..11,13,48}": list(range(1,12))+[13,48],
  "{1..11,13,60}": list(range(1,12))+[13,60],
}
for name,S in tests.items():
    M,arg=M_rigorous(S)
    tag = "TIGHT (=1/14)" if M==F(1,14) else ("LOOSE (>1/14, NOT a tiler)" if M>F(1,14) else "M<1/14 !!")
    print(f"{name:20s}: M = {M} = {float(M):.6f}  at t={arg}  -> {tag}")
print("\n(1/14 = %.6f)" % (1/14))
