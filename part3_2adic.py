"""Part 3: the 2-adic tree and n even. n=4 (m=3), n=6 (m=5).
Tabulate speeds by 2-adic valuation; lonely structure at q=2,4,8,16.
Does 2-adic tree alone give lonely time?"""
from lrc import lonely_count_all_a, is_safe

def v2(x):
    c=0
    while x%2==0: x//=2; c+=1
    return c

def lonely_times_at_q(speeds, q):
    n=len(speeds)+1
    out=[]
    for a in range(1,q):
        if all(is_safe((v*a)%q,q,n) for v in speeds):
            out.append(a)
    return out

# n=4 sets (m=3) and n=6 sets (m=5)
SETS = {
    "n=4 {1,2,3}":[1,2,3],
    "n=4 {1,3,4}":[1,3,4],
    "n=4 {2,3,7}":[2,3,7],
    "n=4 {1,3,5}":[1,3,5],
    "n=4 {3,5,7}":[3,5,7],
    "n=6 {1,2,3,4,5}":[1,2,3,4,5],
    "n=6 {1,3,5,7,9}":[1,3,5,7,9],
    "n=6 {1,2,3,4,6}":[1,2,3,4,6],
}
for name,speeds in SETS.items():
    n=len(speeds)+1
    val=[(s,v2(s)) for s in speeds]
    print(f"\n{name}  n={n} threshold=1/{n}")
    print("  2-adic valuations:", {s:v for s,v in val})
    for k in range(1,5):
        q=2**k
        t=lonely_times_at_q(speeds,q)
        print(f"   q=2^{k}={q}: lonely a in [1,{q-1}] -> {t if t else 'NONE'}")
