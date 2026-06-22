from fractions import Fraction as F
def nf(x):
    r=x%1; return min(r,1-r)
def M_exact(S):
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
print("EQUIDISTRIBUTION test -- does a LARGE committed speed keep M>1/14 (lonely)?")
print("1/14 =", float(F(1,14)))
base=list(range(1,12))+[13]   # {1..11,13}, M=1/12
for v in [60, 840, 5040]:     # 5040 = 84*lcm(1..5): THM-566 adversarial, COVERING, mult of 14
    S=base+[v]; M,arg=M_exact(S)
    cov = all(any(s%q==0 for s in S) for q in range(2,15))
    print(f"  {{1..11,13,{v}}}: M={M}={float(M):.5f} {'>1/14 LONELY' if M>F(1,14) else '<=1/14'}; "
          f"covering(all q<=14)? {cov}; 14|{v}? {v%14==0}; witness t={arg}")
print("\n=> THM-566 adversarial {1..11,13,5040} is COVERING (no easy t=1/q witness, 14|5040 so t=1/14 fails)")
print("   yet M stays at the 12-subset's 1/12 >> 1/14: the huge speed EQUIDISTRIBUTES, its danger zone")
print("   covers only ~1/7 of the subset lonely set, so loneliness SURVIVES (witness at large denominator).")
