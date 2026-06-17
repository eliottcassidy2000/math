"""
tight_locus_largescan_kps.py  (kind-pasteur, 2026-06-16)

Aggressive hunt for tight configs that THREATEN finiteness:
 - large-entry / large-lcm primitive tight 13-sets
 - random search over wide ranges
 - "AP + far outlier" families (does {1..12, W} stay tight as W -> inf?)
 - dilations of the AP and of the sporadic core (these are NON-primitive at d>1,
   but a *perturbed* dilation can be primitive AND tight -> these are the real
   threat to bounded-lcm).

If the tight locus is finite & bounded-lcm, NONE of these should produce a tight
config with lcm/max above the known small ones (beyond the AP & sporadic family).
"""
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations
import random, time, json

def danger_arcs(v):
    w = Fr(1, 14*v); A = []
    for k in range(v+1):
        c = Fr(k, v); lo = c - w; hi = c + w
        if lo < 0: A += [(Fr(0), hi), (1+lo, Fr(1))]
        elif hi > 1: A += [(lo, Fr(1)), (Fr(0), hi-1)]
        else: A.append((lo, hi))
    return A

def L_exact(S):
    A = []
    for v in S: A.extend(danger_arcs(v))
    A = sorted((a, b) for a, b in A if b > a)
    tot = Fr(0); cl = ch = None
    for a, b in A:
        if ch is None: cl, ch = a, b
        elif a <= ch: ch = max(ch, b)
        else: tot += ch - cl; cl, ch = a, b
    if ch is not None: tot += ch - cl
    return 1 - tot

def gcd_list(S): return reduce(gcd, S)
def lcm2(a,b): return a*b//gcd(a,b)
def lcm_list(S): return reduce(lcm2, S)
def primitive(S): return gcd_list(S)==1
def is_tight(S): return L_exact(tuple(S))==0

def main():
    t0=time.time()
    tight=[]; seen=set()
    def consider(S, src):
        S=tuple(sorted(S))
        if S in seen or len(set(S))!=13 or any(x<=0 for x in S): return
        seen.add(S)
        if not primitive(S): return
        if is_tight(S):
            tight.append({"S":S,"lcm":lcm_list(S),"max":max(S),"source":src})
            print(f"  TIGHT: lcm={lcm_list(S)} max={max(S)} S={list(S)} [{src}]", flush=True)

    # 1) {1..12, W} far-outlier family: does it ever go tight as W grows?
    print("=== {1..12, W} far-outlier, W=14..2000 ===", flush=True)
    base12=list(range(1,13))
    for W in range(14,2001):
        consider(base12+[W], f"{{1..12,{W}}}")
    # 2) {1..11,13,W}: the sporadic family (W=24 known tight). scan W to 2000
    print("=== {1..11,13,W}, W=14..2000 ===", flush=True)
    b=list(range(1,12))+[13]
    for W in range(14,2001):
        if W==13: continue
        consider(b+[W], f"{{1..11,13,{W}}}")
    # 3) drop-2-add-2 from AP but allow ONE very large entry up to 2000
    print("=== drop-2-add-2 AP, one entry small (<=30) one large (<=2000) ===", flush=True)
    AP=list(range(1,14))
    cnt=0
    for di,dj in combinations(range(13),2):
        kept=[AP[i] for i in range(13) if i not in (di,dj)]
        ks=set(kept)
        smalls=[x for x in range(1,31) if x not in ks]
        for a in smalls:
            ks2=ks|{a}
            for big in range(31,2001,1):
                if big in ks2: continue
                consider(kept+[a,big], "d2a2-onelarge")
                cnt+=1
    print(f"  scanned {cnt}", flush=True)
    # 4) perturbed dilations: d*(AP \ {j}) ∪ {w}, d up to 12, w arbitrary up to d*60
    print("=== perturbed dilations d*(AP\\{j}) ∪ {w}, d<=12 ===", flush=True)
    for d in range(1,13):
        for j in range(1,14):
            core=[d*x for x in range(1,14) if x!=j]
            for w in range(1, d*60+1):
                consider(core+[w], f"pdil d={d} j={j} w={w}")
    print(f"  tight so far {len(tight)}", flush=True)
    # 5) perturbed dilations of SPORADIC core {1..11,13,24}
    print("=== perturbed dilations of sporadic core, d<=12 ===", flush=True)
    spor=list(range(1,12))+[13,24]
    for d in range(1,13):
        for j in range(13):
            core=[d*spor[i] for i in range(13) if i!=j]
            for w in range(1, d*60+1):
                consider(core+[w], f"pdil-spor d={d} j={j} w={w}")
    print(f"  tight so far {len(tight)}", flush=True)
    # 6) random wide search
    print("=== random wide search (200k trials, entries 1..300) ===", flush=True)
    random.seed(14)
    for _ in range(200000):
        S=random.sample(range(1,301),13)
        if primitive(S) and is_tight(S):
            consider(S,"random-wide")
    print(f"  tight so far {len(tight)}", flush=True)
    # 7) random near-AP: AP + small jitter
    print("=== random near-AP jitter (200k trials) ===", flush=True)
    for _ in range(200000):
        S=set()
        for x in range(1,14):
            S.add(x+random.randint(-3,30))
        S={s for s in S if s>=1}
        if len(S)!=13: continue
        S=sorted(S)
        if primitive(S) and is_tight(S):
            consider(S,"random-nearAP")
    print(f"  tight so far {len(tight)}", flush=True)

    tight_sorted=sorted(tight,key=lambda d:(d["lcm"],d["max"]))
    print("\n===== ALL TIGHT (largescan) =====", flush=True)
    for d in tight_sorted:
        print(f"  lcm={d['lcm']:>10} max={d['max']:>5} S={list(d['S'])} [{d['source']}]", flush=True)
    print(f"\nTotal: {len(tight)}", flush=True)
    if tight:
        print(f"Max lcm: {max(d['lcm'] for d in tight)}  Max entry: {max(d['max'] for d in tight)}", flush=True)
    out={"tight":[{"S":list(d["S"]),"lcm":d["lcm"],"max":d["max"],"source":d["source"]} for d in tight_sorted],
         "n_tight":len(tight),"elapsed_sec":time.time()-t0}
    with open("05-knowledge/results/tight_locus_largescan_kps.json","w") as f:
        json.dump(out,f,indent=2)
    print(f"Elapsed {time.time()-t0:.1f}s", flush=True)

if __name__=="__main__":
    main()
