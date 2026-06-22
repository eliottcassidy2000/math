"""
Attacking the majorization p0(E) <= p0(AP) directly (kind-pasteur S31k).

Fourier expansion of the cover atom over the resonance lattice:
  p0(E) = sum_{A subset {1..6}} (-1)^|A| J(A,E),   J(A,E) = int_0^1 prod_e f_A(e x) dx,
  f_A(t) = 1 - sum_{j in A} 1_{I_j}(t),  I_j=[j/7,(j+1)/7),  hat f_{A,0}=1-|A|/7,
  hat f_{A,n} = -gamma_n * sum_{j in A} e^{-2pi i n j/7},  gamma_n=(1-e^{-2pi i n/7})/(2pi i n).
Expanding prod_e f_A(e x) and integrating keeps only Sum_e n_e e = 0 (the resonances):
  J(A,E) = sum_{(n_e): sum n_e e=0} prod_e hat f_{A,n_e}.
Diagonal (all n_e=0): (1-|A|/7)^k -> decorrelated p0.  Leading correction = the
ADDITIVE-ENERGY resonances (n_a=n_b=1, n_c=n_d=-1, e_a+e_b=e_c+e_d), Fourier frequency 1,
each contributing |hat f_{A,1}|^4 (1-|A|/7)^{k-4}.  So
  p0(E) - p0_decorr ~ A*(E) * Gamma_k,  Gamma_k := sum_A (-1)^|A| |hat f_{A,1}|^4 (1-|A|/7)^{k-4},
where A*(E) = (# nontrivial additive quadruples).  If Gamma_k > 0, p0 INCREASES with the
additive energy => maximized by the interval (sharp classical A(E)<=A(AP)).  Computed here.
"""
import cmath, math
from itertools import combinations, chain

def subsets(s):
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

w = cmath.exp(-2j*math.pi/7)   # e^{-2pi i/7}
gamma1 = (1 - cmath.exp(-2j*math.pi/7)) / (2j*math.pi*1)
sectors = [1,2,3,4,5,6]

def S_A(A):                     # sum_{j in A} e^{-2pi i j/7}
    return sum(w**j for j in A)

def Gamma(k):
    tot = 0.0
    for A in subsets(sectors):
        a = len(A)
        hatf1 = -gamma1 * S_A(A)          # hat f_{A,1}
        tot += ((-1)**a) * (abs(hatf1)**4) * ((1 - a/7.0)**(k-4))
    return tot

def p0_decorr(k):
    return sum(((-1)**len(A)) * ((1-len(A)/7.0)**k) for A in subsets(sectors))

if __name__=="__main__":
    print("Gamma_k = additive-energy coefficient in p0 (sign => direction of extremality):")
    for k in range(8,13):
        G = Gamma(k); pd = p0_decorr(k)
        print(f"  k={k}: Gamma_k = {G:+.6e}   (>0 => p0 increases with A(E) => INTERVAL maximizes)   p0_decorr={pd:.5f}")
    print("\nValidation: actual p0(E)-p0_decorr vs predicted (additive-energy * Gamma_k):")
    # reuse exact-ish p0
    def p0(E):
        E=sorted(set(e for e in E if e!=0)); bset={0.0,1.0}
        for e in E:
            for j in range(8):
                b=j/7.0; m=0
                while True:
                    xv=(b+m)/e
                    if xv>=1: break
                    if xv>=0: bset.add(xv)
                    m+=1
        B=sorted(bset); tot=0.0
        for lo,hi in zip(B,B[1:]):
            if hi<=lo: continue
            mid=(lo+hi)/2
            if len(set(int((e*mid)%1*7) for e in E)&set(range(1,7)))==6: tot+=hi-lo
        return tot
    def Astar(E):   # nontrivial additive quadruples (ordered): a+b=c+d, {a,b}!={c,d}
        from collections import Counter
        E=list(E); c=Counter(a+b for a in E for b in E)
        total=sum(v*v for v in c.values())
        trivial=2*len(E)**2 - len(E)   # (a,b,a,b) and (a,b,b,a) minus double-count diagonal
        return total - trivial
    k=9; G=Gamma(9); pd=p0_decorr(9)
    for E in [list(range(9)), [0,1,2,3,4,5,6,7,9], [0,2,4,6,8,10,12,14,16], [0,1,3,7,12,20,30,44,60]]:
        actual = p0(E)-pd; pred = Astar(E)*G
        print(f"  E={E}: A*={Astar(E):5d}  actual dev={actual:+.4f}  pred(lead)={pred:+.4f}")
