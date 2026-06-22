import cmath, math
from itertools import combinations, chain
def subsets(s): return chain.from_iterable(combinations(s,r) for r in range(len(s)+1))
w = cmath.exp(-2j*math.pi/7)
def gam(n): return (1-cmath.exp(-2j*math.pi*n/7))/(2j*math.pi*n)
sectors=[1,2,3,4,5,6]
def S_A(A): return sum(w**j for j in A)
def Gamma_s(k, s, freq=1):
    g=abs(gam(freq))
    tot=0.0
    for A in subsets(sectors):
        a=len(A)
        SA=sum(cmath.exp(-2j*math.pi*freq*j/7) for j in A)
        hatf=-gam(freq)*SA  # but use |.|^{2s}
        tot += ((-1)**a)*(abs(hatf)**(2*s))*((1-a/7.0)**(k-2*s))
    return tot
print("Frequency-1 additive-moment coefficients Gamma_k^(s) (s=2..6); sign => term-by-term direction:")
for k in [9,12]:
    row=[]
    for s in range(2,7):
        if k-2*s < 0: row.append((s,None)); continue
        row.append((s, Gamma_s(k,s,1)))
    print(f"  k={k}: " + "  ".join(f"s={s}:{(f'{v:+.3e}' if v is not None else 'NA')}" for s,v in row))
print("\nFrequency-2 leading (s=2) coefficient (the tail beyond freq 1), for size comparison:")
for k in [9,12]:
    print(f"  k={k}: Gamma_k^(2),freq2 = {Gamma_s(k,2,2):+.3e}   vs freq1 {Gamma_s(k,2,1):+.3e}  "
          f"(ratio {abs(Gamma_s(k,2,2)/Gamma_s(k,2,1)):.3f})")
