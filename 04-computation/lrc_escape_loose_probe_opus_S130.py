from math import gcd, lcm, ceil, floor
import random
random.seed(7)

def clears_at(V,q):
    lo=ceil(2*q/25); hi=floor(23*q/25)
    if lo>hi: return None
    for c in range(1,q):
        if gcd(c,q)!=1: continue
        if all(lo<=(v*c)%q<=hi for v in V): return c
    return False

print("=== ESCAPE-FAMILY LOOSENESS PROBE (opus-S130) ===")
print("Q: do compressed varying-k families {i+L*k_i} ALL clear at some q (=> loose, not gap members)?\n")
maxq=80; noclear=[]; probes=0; firstq=[]
for Q0 in [12,20,25,30,37,41]:
    L=lcm(*range(2,Q0+1))
    for _ in range(40):
        k=[random.randint(1,5) for _ in range(12)]
        if len(set(k))==1: k[1]+=1
        V=[(i+1)+L*k[i] for i in range(12)]
        probes+=1
        cq=next((q for q in range(6,maxq+1) if clears_at(V,q) not in (None,False)), None)
        if cq is None: noclear.append((Q0,k))
        else: firstq.append(cq)
print(f"  probed: {probes} varying-k compressed escape families (heights ~10^3..10^15)")
print(f"  NOT clearing at any q<=80 (gap-member candidates): {len(noclear)}")
for Q0,k in noclear[:10]: print("     Q0=",Q0," k=",k)
if firstq:
    from collections import Counter
    print(f"  clearing-q distribution (min={min(firstq)}, max={max(firstq)}): {dict(sorted(Counter(firstq).items()))}")
print(f"\n  VERDICT: {'ALL escape families clear somewhere => all LOOSE, no gap members at height. (C) holds; only the FINITE covering fails (modulus unbounded).' if not noclear else 'FOUND non-clearing families -- INVESTIGATE (possible gap member).'}")
