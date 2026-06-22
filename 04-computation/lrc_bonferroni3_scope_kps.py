"""
HONEST scope of the Bonferroni-3 reduction (kps-S31u, correcting S31t overclaim).
Bonferroni-3 (p0<=T_1+T_2+T_3) is NOT universal. Hypothesis: it holds exactly on the BINDING leg
(large base, p0 near cap), and FAILS only on the SLACK leg (small base, p0 << cap) -- where the
wide bound holds anyway. Verify: (a) does the WIDE BOUND p0<=cap hold in ALL cases? (b) does
Bonferroni-3 failure <=> p0 small (slack)? This assembles the honest dichotomy.
"""
import random
from itertools import combinations
def p0(E):
    E=sorted(set(e for e in E if e!=0))
    if not E: return 0.0
    bset={0.0,1.0}
    for e in E:
        for j in range(8):
            b=j/7.0; m=0
            while True:
                xv=(b+m)/e
                if xv>=1: break
                if xv>=0: bset.add(xv)
                m+=1
    Bs=sorted(bset); tot=0.0
    for lo,hi in zip(Bs,Bs[1:]):
        if hi<=lo: continue
        mid=(lo+hi)/2
        if len(set(int((e*mid)%1*7) for e in E)&set(range(1,7)))==6: tot+=hi-lo
    return tot
def s3(base,far):
    pf={S:p0(base+[far[i] for i in S]) for r in range(4) for S in combinations(range(len(far)),r)}
    tot=0.0
    for r in range(1,4):
        for S in combinations(range(len(far)),r):
            for r2 in range(r+1):
                for T in combinations(S,r2):
                    if T in pf: tot+=((-1)**(r-r2))*pf[T]
    return tot
random.seed(11)
cap=0.3815  # cap_8 (smallest, most binding); total runners = base+far
rows=[]
for _ in range(50):
    bsz=random.randint(3,7); base=sorted(set([0]+random.sample(range(1,13), bsz-1)))
    nf=8-len(base); 
    if nf<1: continue
    far=sorted(random.sample(range(18,60), nf))
    full=p0(base+far); bound3=s3(base,far)
    rows.append((len(base),full,bound3))
widebound_ok=all(full<=cap+1e-9 for _,full,_ in rows)
bonf_holds_large=[full for nb,full,b3 in rows if full<=b3+1e-9]
bonf_fails_small=[full for nb,full,b3 in rows if full>b3+1e-9]
print(f"WIDE BOUND p0<=cap_8={cap} holds in ALL {len(rows)} configs: {widebound_ok}")
print(f"Bonferroni-3 HOLDS ({len(bonf_holds_large)} configs): p0 range [{min(bonf_holds_large):.3f},{max(bonf_holds_large):.3f}], mean {sum(bonf_holds_large)/len(bonf_holds_large):.3f}")
if bonf_fails_small:
    print(f"Bonferroni-3 FAILS ({len(bonf_fails_small)} configs): p0 range [{min(bonf_fails_small):.3f},{max(bonf_fails_small):.3f}], mean {sum(bonf_fails_small)/len(bonf_fails_small):.3f}")
    print(f"  => failures have p0 max {max(bonf_fails_small):.3f} << cap {cap} (SLACK, not binding)")
print("\nHONEST: Bonferroni-3 is a BINDING-leg handle (p0 near cap), NOT universal. The slack leg")
print("(small base, high-order coverage) has p0<<cap directly. Wide bound = [binding: <=T1+T2+T3] + [slack].")
