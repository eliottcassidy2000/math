r"""
The gamma-trick for the 14-covering residual (kps-S31ad, creative).
A multiple of 14 is 1/14-PERIODIC in t: ||14m(t+1/14)|| = ||14m t||. So fix gamma*=frac(14 t) making
ALL multiples of 14 safe (via LRC on {m_i}); then the 14 points {(gamma*+j)/14 : j=0..13} keep the
multiples safe, and we need R safe at ONE. PIGEONHOLE: each R-runner r_k coprime to 14 marks <=2 of the
14 points as bad (the 14 points hit the 14th-roots, only 2 within 1/14); so |R|<=6 => <=12 bad => a GOOD
point => M(S)>=1/14. Closes r>=7 (|R|<=6, R coprime-to-14). Residual: R has a multiple of 7 (gcd=7 => up
to 7 bad), the apex-7 inner self-reference of 14=2*7.
"""
from fractions import Fraction as F
from math import gcd
def nf(x):
    r=x%1; return min(r,1-r)
def M_set(S):  # max-min over rational crit points
    S=sorted(set(abs(s) for s in S if s!=0)); C=set()
    for i in range(len(S)):
        for j in range(i,len(S)):
            for comb in {S[i]+S[j], abs(S[i]-S[j])}:
                if comb:
                    for k in range(1,comb): C.add(F(k,comb))
    return max((min(nf(s*t) for s in S) for t in C if 0<t<1), default=F(0))
def gamma_star(mults):  # frac(14t)=gamma making multiples safe = witness of {m_i}={mult/14}
    ms=[m//14 for m in mults]
    # find gamma in (0,1) maximizing min ||m_i gamma||
    C=set()
    for i in range(len(ms)):
        for j in range(i,len(ms)):
            for comb in {ms[i]+ms[j],abs(ms[i]-ms[j])}:
                if comb:
                    for k in range(1,comb): C.add(F(k,comb))
    best=F(0); g=None
    for t in C:
        if 0<t<1:
            m=min(nf(mm*t) for mm in ms)
            if m>best: best=m; g=t
    return g,best
print("gamma-trick + pigeonhole on 14-covering sets:")
tests=[("r=7, R coprime-14",[14,28,42,56,70,84,98]+[1,3,5,9,11,13]),       # 7 mult of 14 + 6 coprime
       ("r=8, R coprime-14",[14,28,42,56,70,84,98,112]+[1,3,5,9,11]),
       ("r=7, R has 7 (mult of 7!)",[14,28,42,56,70,84,98]+[7,1,3,5,9,11]),
       ("r=2 (union bound regime)",[14,28]+[1,2,3,4,5,6,8,9,10,11,13])]
for name,S in tests:
    S=sorted(set(S))
    if len(S)!=13: 
        print(f"  {name}: |S|={len(S)} skip"); continue
    mults=[s for s in S if s%14==0]; R=[s for s in S if s%14!=0]
    g,mm=gamma_star(mults)
    # the 14 points (gamma*+j)/14; count R-bad (some r_k within 1/14)
    bad=0; good_pts=[]
    for j in range(14):
        p=(g+j)/14
        if all(nf(r*p)>=F(1,14) for r in R): good_pts.append(j)
        else: bad+=1
    Rcop=all(gcd(r,14)==1 for r in R)
    M=M_set(S)
    print(f"  {name}: |R|={len(R)} R-coprime-14? {Rcop}; multiples safe at gamma*={g} (margin {mm}); "
          f"bad pts={bad}/14, GOOD pts={good_pts}; M(S)={M}={float(M):.4f} {'>1/14 OK' if M>F(1,14) else ''}")
print("\n=> R coprime-to-14, |R|<=6 (r>=7): <=2|R|<=12 bad => a GOOD point exists => M>=1/14 (pigeonhole).")
print("   Union bound closes r<=6. RESIDUAL: R has a multiple of 7 (gcd=7, up to 7 bad) -- the apex-7")
print("   self-reference (14=2*7, the 7-content recurses). Narrows the 14-covering crux to double-7 sets.")
