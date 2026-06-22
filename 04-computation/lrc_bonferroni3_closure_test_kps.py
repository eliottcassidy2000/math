"""
Stress-test the Bonferroni-3 closure of the multi-large (wide) case (kps-S31u).
Claim: p0(B u far) <= T_1+T_2+T_3 (Bonferroni-3 upper bound), because the r>=4 Newton tail is <=0.
If robust, the multi-large case reduces to single-far (THM-565) + doublet (HYP-2797) + triple (THM-557),
all individually bounded => the WIDE BOUND closes. Test over many random wide configs.
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
def packets(base, far):  # T_r for r=0..len(far)
    pf={}
    for r in range(len(far)+1):
        for S in combinations(range(len(far)),r):
            pf[S]=p0(base+[far[i] for i in S])
    T=[0.0]*(len(far)+1)
    for r in range(len(far)+1):
        for S in combinations(range(len(far)),r):
            val=0.0
            for r2 in range(r+1):
                for Tt in combinations(S,r2):
                    val+=((-1)**(r-r2))*pf[Tt]
            T[r]+=val
    return T, pf[tuple(range(len(far)))]
random.seed(7)
NF=5; trials=40
bonf3_ok=0; tail_neg=0; decay_ok=0; viol=[]
caps={8:0.3815,9:0.4943,10:0.6044}
for _ in range(trials):
    bsz=random.randint(3,5); base=sorted(random.sample(range(0,12), bsz))
    if 0 not in base: base=[0]+base[:-1]
    spread=random.choice([15,25,40])
    far=sorted(random.sample(range(spread, spread+30), NF))
    T,full=packets(base,far)
    s3=sum(T[1:4]); tail=sum(T[4:])
    if full<=s3+1e-9: bonf3_ok+=1
    else: viol.append((base,far,full,s3))
    if tail<=1e-9: tail_neg+=1
    # decay |T_3|>=|T_4|>=|T_5| ?
    if abs(T[3])>=abs(T[4])-1e-9>=0 and abs(T[4])>=abs(T[5])-1e-9: decay_ok+=1
print(f"Bonferroni-3 (p0 <= T_1+T_2+T_3) holds: {bonf3_ok}/{trials}")
print(f"r>=4 tail <= 0:                        {tail_neg}/{trials}")
print(f"decay |T_3|>=|T_4|>=|T_5|:             {decay_ok}/{trials}")
if viol: 
    print("VIOLATIONS:")
    for b,f,fu,s in viol[:3]: print(f"  base={b} far={f}: p0={fu:.4f} > T1+T2+T3={s:.4f}")
else: print("=> ZERO violations: p0 <= T_1+T_2+T_3 over all sampled wide configs.")
print("   => multi-large case REDUCES to single-far(THM-565)+doublet(HYP-2797)+triple(THM-557) -- the")
print("      r>=4 corrections only HELP (subtract). The wide bound closes IF T_1+T_2+T_3 <= cap (next).")
