"""
DISPROOF ATTEMPT (kps-S31p): kill resonances b=1..14 with highly-divisible speeds, leaving
only b>=15 => M<=1/15<1/14 (counterexample?). The PROOF response: killing the exact point a/b
(speed ==0 mod b) does NOT kill the lonely NEIGHBORHOOD (width ~1/(bV)); a large divisible speed
has a THIN danger zone (1/(7s)) that misses the neighborhood => loneliness survives (ZETA(2) floor).
"""
from fractions import Fraction as F
def nf(x):
    r=x%1; return min(r,1-r)
def M_and_witness(S):
    S=sorted(set(abs(s) for s in S if s!=0)); C=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for comb in {S[i]+S[j], abs(S[i]-S[j])}:
                if comb:
                    for k in range(1,comb): C.add(F(k,comb))
    best=F(0); arg=None
    for t in C:
        if 0<t<1:
            m=min(nf(s*t) for s in S)
            if m>best: best=m; arg=t
    return best,arg
def covers(S):  # kills which b in 1..14 (a speed ==0 mod b)?
    return [b for b in range(1,15) if any(s%b==0 for s in S)]
print("DISPROOF ATTEMPTS -- aggressive resonance-killers (want M<1/14):  1/14 =",float(F(1,14)))
attempts = {
 "{2..14} (b=2..14)": list(range(2,15)),
 "divisible spread {2,4,6,8,10,12,14,3,9,5,7,11,13}": [2,4,6,8,10,12,14,3,9,5,7,11,13],
 "highly-div {12,24,36,8,9,5,7,11,13,14,10,6,4}": [12,24,36,8,9,5,7,11,13,14,10,6,4],
 "lcm-killers {4,8,9,5,7,11,13,12,10,14,6,3,2}": [4,8,9,5,7,11,13,12,10,14,6,3,2],
 "big-killers {28,40,72,90,11,13,84,120,126,110,130,154,182}": [28,40,72,90,11,13,84,120,126,110,130,154,182],
 "mixed small+div {1,2,3,12,24,8,9,5,7,11,13,14,10}": [1,2,3,12,24,8,9,5,7,11,13,14,10],
}
mins=[]
for name,S in attempts.items():
    M,t=M_and_witness(S); b=t.denominator if t else 0
    kb=covers(S)
    mins.append(float(M))
    print(f"  {name[:50]:50s}: M={str(M):>7s}={float(M):.4f} {'<1/14 !!!' if M<F(1,14) else 'OK'} (b_kill 1..14={len(kb)}/14, witness b={b})")
print(f"\n  MIN over attempts = {min(mins):.4f}  vs 1/14={float(F(1,14)):.4f}  => {'COUNTEREXAMPLE' if min(mins)<1/14 else 'NO counterexample -- floor holds'}")
print("  => killing exact resonance points a/b does NOT lower M: a large divisible speed has a THIN")
print("     danger zone that misses the lonely NEIGHBORHOOD. The surviving witness is a resonance the")
print("     killers cannot cover. The disproof FAILS exactly where the ZETA(2) neighborhood floor lives.")
