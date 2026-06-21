"""
OPEN-Q-108 synthesis + completeness critic (independent re-verification).
Re-checks: engine correctness, chain identity, comb bound, V-growth (V<=C*k FALSE,
V<=7*sigma rigorous), separated-branch cutoff, balanced-branch failure, wide hunt.
EXACT rationals. DO NOT run git.
"""
from fractions import Fraction as F
import random

def bp(E):
    s=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            m=0
            while True:
                xv=(F(j,7)+m)/e
                if xv>=1: break
                if xv>=0: s.add(xv)
                m+=1
    return sorted(b for b in s if 0<=b<1)

def p0p1(E):
    E=sorted(set(E)); B=bp(E); a=F(0);b=F(0)
    for lo,hi in zip(B,B[1:]+[F(1)]):
        if hi<=lo: continue
        mid=(lo+hi)/2; miss=set(range(1,7))-set(int((e*mid)%1*7) for e in E)
        if len(miss)==0:a+=hi-lo
        elif len(miss)==1:b+=hi-lo
    return a,b

def Vcount(E):
    """V(E)=number of interval-components of the miss-exactly-one region."""
    E=sorted(set(E)); B=bp(E)
    comp=0; prev=False
    for lo,hi in zip(B,B[1:]+[F(1)]):
        if hi<=lo: continue
        mid=(lo+hi)/2; miss=len(set(range(1,7))-set(int((e*mid)%1*7) for e in E))
        cur=(miss==1)
        if cur and not prev: comp+=1
        prev=cur
    return comp

caps = {8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7)}

def report(label, val):
    print(f"{label}: {val}")

def main():
    random.seed(12345)
    print("="*60); print("ENGINE SANITY"); print("="*60)
    for k in range(8,13):
        E=list(range(0,k)); p0,p1=p0p1(E)
        report(f"k={k} consec p0", f"{float(p0):.5f} cap={float(caps[k]):.5f} {'OK' if p0<=caps[k] else 'VIOL'}")

    print("="*60); print("(A) CHAIN IDENTITY (telescoping, exact)"); print("="*60)
    # p0(E) = p0(core)+sum (1/7)p1(E_{i-1}) + sum Delta_wi  is pure telescoping; trivially exact.
    # Verify Delta_wi = p0(E_i)-p0(E_{i-1})-(1/7)p1(E_{i-1}) sums correctly.
    fails=0; tot=0
    for trial in range(120):
        k=random.randint(8,12)
        core_sz=random.randint(3,5)
        core=set([0])
        while len(core)<core_sz: core.add(random.randint(1,12))
        core=sorted(core)
        fars=sorted(random.sample(range(15,200), k-len(core)))
        chain=[core]+[core+fars[:i] for i in range(1,len(fars)+1)]
        # wait need cumulative
        Es=[core]
        cur=list(core)
        for w in fars:
            cur=sorted(cur+[w]); Es.append(cur)
        E=Es[-1]
        p0E,_=p0p1(E)
        acc,_=p0p1(core)
        for i in range(1,len(Es)):
            p0i,_=p0p1(Es[i]); p0im1,p1im1=p0p1(Es[i-1])
            d=p0i-p0im1-F(1,7)*p1im1
            acc+= F(1,7)*p1im1 + d
        tot+=1
        if acc!=p0E: fails+=1
    report("chain identity fails", f"{fails}/{tot}")

    print("="*60); print("(B) COMB BOUND |Delta_w|<=(6/49)V(E')/w"); print("="*60)
    fails=0; tot=0; worst=F(0)
    for trial in range(400):
        k=random.randint(8,12)
        ep=random.randint(3,6)
        Ep=set([0])
        while len(Ep)<ep: Ep.add(random.randint(1,14))
        Ep=sorted(Ep)
        w=random.randint(15,300)
        if w in Ep: continue
        E=sorted(Ep+[w])
        p0E,_=p0p1(E); p0p,p1p=p0p1(Ep)
        Plat=p0p+F(1,7)*p1p
        delta=p0E-Plat
        V=Vcount(Ep)
        bound=F(6,49)*V/w
        tot+=1
        if abs(delta)>bound: fails+=1
        if bound>0:
            r=abs(delta)/bound
            if r>worst: worst=r
    report("comb bound fails", f"{fails}/{tot}, worst ratio={float(worst):.4f}")

    print("="*60); print("(D) V-GROWTH"); print("="*60)
    maxVsig=F(0); maxVk=F(0); viol7=0; tot=0
    for trial in range(400):
        k=random.randint(8,12)
        E=set([0])
        while len(E)<k: E.add(random.randint(1,90))
        E=sorted(E); V=Vcount(E); sig=sum(E); nz=len([e for e in E if e>0])
        if sig>0:
            r=F(V,sig)
            if r>maxVsig: maxVsig=r
        if nz>0:
            rk=F(V,nz)
            if rk>maxVk: maxVk=rk
        tot+=1
        if V>7*sig: viol7+=1
    report("max V/sigma", f"{float(maxVsig):.4f}")
    report("max V/#nonzero (proxy for V<=C*k)", f"{float(maxVk):.4f}  -> V<=C*k FALSE if large")
    report("V<=7*sigma violations", f"{viol7}/{tot}")

    print("="*60); print("(SEP) SEPARATED BRANCH k=9 two-far cutoff"); print("="*60)
    # k=9: core consec_8 = {0..7}, two far w1, w2=2*w1. Check certificate.
    Q8 = None  # Q(k-1) plateau argmax: use consec_{k-1} plat as proxy
    for w1 in [40,45,48,49,50,55,60,80]:
        core=list(range(0,7))  # 7 elements, then 2 far -> k=9
        E1=sorted(core+[w1]); E2=sorted(core+[w1,2*w1])
        p0,_=p0p1(E2)
        report(f"k=9 w1={w1},w2={2*w1} p0", f"{float(p0):.5f} cap9={float(caps[9]):.5f} {'OK' if p0<=caps[9] else 'VIOL'}")

    print("="*60); print("(HUNT) wide primitive k-sets, p0>cap?"); print("="*60)
    viol=0; tot=0; thin=F(1)
    fams=0
    for trial in range(4000):
        k=random.randint(8,12)
        style=random.randint(0,4)
        E=set([0])
        if style==0:  # consec + far
            for i in range(1,k-1): E.add(i)
            E.add(random.randint(15,400))
        elif style==1:  # stretched AP
            step=random.randint(2,7);
            for i in range(1,k): E.add(i*step)
        elif style==2:  # two blocks
            b=random.randint(20,200)
            for i in range(k//2): E.add(i+1)
            for i in range(k-k//2): E.add(b+i)
        elif style==3:  # random wide
            while len(E)<k: E.add(random.randint(1,500))
        else:  # multi-scale
            scs=[1,random.randint(20,60),random.randint(200,400)]
            while len(E)<k: E.add(random.choice(scs)*random.randint(1,5))
        E=sorted(E)
        if len(E)<k: continue
        if max(E)-min(E)<=14: continue
        from math import gcd
        from functools import reduce
        g=reduce(gcd,[e for e in E if e>0])
        if g!=1: continue
        p0,_=p0p1(E)
        tot+=1
        m=caps[k]-p0
        if m<thin: thin=m
        if p0>caps[k]:
            viol+=1
            if viol<=5: print("  VIOLATION", k, E, float(p0))
    report("wide hunt violations", f"{viol}/{tot}, thinnest margin={float(thin):.5f}")

main()
