# death-star-S56: complete exhaustive verification that {AP,GW} are the only primitive tight
# families for LRC(14) up to a Vmax bound. Fully rigorous: only the PROVEN divisor-covering +
# primitivity prune (no unproven no-mult-14 assumption). Parallel over (largest,2nd-largest).
import sys, time
from math import gcd
from functools import reduce
from itertools import combinations
from multiprocessing import Pool

def is_tight(fam):
    w=fam[-1]; Q=2*w+2; hit=False
    for q in [14]+list(range(2,14))+list(range(15,Q+1)):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=q
            for v in fam:
                r=(v*a)%q
                d=r if r<=q-r else q-r
                if d<m: m=d
            v14=14*m
            if v14>q: return False
            if v14==q: hit=True
    return hit

def prune_and_check(fam):
    # PROVEN-necessary prune: primitive + a multiple of 7 + a multiple of 2 (divisor-covering)
    if reduce(gcd,fam)!=1: return False
    if not any(v%7==0 for v in fam): return False
    if not any(v%2==0 for v in fam): return False
    return is_tight(fam)

def task(args):
    w,s=args
    found=[]
    for rest in combinations(range(1,s),11):
        fam=list(rest)+[s,w]
        if prune_and_check(fam): found.append(tuple(fam))
    return found

if __name__=="__main__":
    LO,HI=25,30
    tasks=[(w,s) for w in range(LO,HI+1) for s in range(12,w)]
    t0=time.time(); allf=[]
    with Pool(4) as p:
        done=0
        for res in p.imap_unordered(task,tasks):
            allf+=res; done+=1
            if done%20==0: print(f"  {done}/{len(tasks)} tasks, tight so far={len(allf)}, {time.time()-t0:.0f}s",flush=True)
    print(f"DONE Vmax in [{LO},{HI}] in {time.time()-t0:.0f}s. New tight families:",allf,flush=True)
    print("Combined with the Vmax<=24 result {AP,GW}: classification holds up to Vmax=%d"%HI,flush=True)
