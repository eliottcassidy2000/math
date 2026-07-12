from math import gcd
from fractions import Fraction as F
import random
def lcm(a,b): return a//gcd(a,b)*b
def posbunch(E):
    nz=[abs(e) for e in E if e]; has0=len(nz)<len(E)
    L=1
    for e in nz: L=lcm(L,e)
    D=7*L; pts=set([0,D])
    for e in nz:
        st=L//e; pts.update(range(0,D+1,st))
    pts=sorted(pts); pn=[0]*8; b0=1 if has0 else 0
    for t1,t2 in zip(pts,pts[1:]):
        s=t1+t2; hit=b0
        for e in nz: hit|=1<<((7*e*s//(2*D))%7)
        pn[7-bin(hit).count("1")]+=t2-t1
    p=[F(x,D) for x in pn]; T=[sum(p[j:]) for j in range(8)]
    return 6*T[1]+4*T[2]+2*T[3], 2*T[5]+4*T[6]
def norm(E):
    g=0
    for e in E: g=gcd(g,e)
    return tuple(sorted(e//g for e in E)) if g>1 else tuple(sorted(E))
def hillclimb(objective, sign, nseeds=40, moves=300, box=90):
    # sign=+1 maximize, -1 minimize
    best=None;best_s=None
    for _ in range(nseeds):
        E=list(norm(tuple(sorted(random.sample(range(1,box),9)))))
        cur=objective(tuple(E))
        for _ in range(moves):
            i=random.randrange(9); old=E[i]
            new=random.randint(1,box)
            if new in E: continue
            E2=sorted(E[:i]+[new]+E[i+1:])
            if len(set(E2))<9: continue
            v=objective(tuple(E2))
            if (v-cur)*sign>0:
                E=E2; cur=v
        En=norm(tuple(E))
        vv=objective(En)
        if best is None or (vv-best)*sign>0: best=vv;best_s=En
    return best,best_s
random.seed(11)
POSf=lambda E: posbunch(E)[0]; BUNf=lambda E: posbunch(E)[1]
mp,mps=hillclimb(POSf,-1)  # minimize POS
mb,mbs=hillclimb(BUNf,+1)  # maximize BUNCH
consecPOS=posbunch(tuple(range(1,10)))[0]; mod7B=F(6,19)
print("ADVERSARIAL HILL-CLIMB (guard vs another box artifact):")
print(f"  min POS found = {float(mp):.5f} at {mps}")
print(f"    consec POS = {float(consecPOS):.5f}  => {'consec still min' if mp>=consecPOS else '*** BEATS consec: '+str(mps)}")
print(f"  max BUNCH found = {mb} ~ {float(mb):.5f} at {mbs}")
print(f"    mod-7 pole BUNCH = 6/19 = {float(mod7B):.5f}  => {'mod-7 still max' if mb<=mod7B else '*** BEATS mod-7: '+str(mbs)}")
FL=F(432,91)
print(f"  corrected separation minPOS-maxBUNCH = {float(min(mp,consecPOS)-max(mb,mod7B)):.5f} vs floor {float(FL):.5f}: "
      f"margin {float(min(mp,consecPOS)-max(mb,mod7B)-FL):+.5f}")
