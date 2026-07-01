from itertools import combinations
from math import factorial
from collections import Counter
exec(open("sccount.py").read().split("A={")[0])  # reuse fix_comp, partitions, rep_perm, class_size
def contributing(n):
    out=[]
    for p in partitions(n):
        if fix_comp(rep_perm(p,n),n)>0: out.append(tuple(p))
    return out
print("Contributing anti-automorphism cycle types vs 'parts=2 mod4 (+ one fixed pt if odd)':")
for n in range(4,15):
    ct=set(contributing(n))
    # predicted: partitions of (n if even else n-1) into parts in {2,6,10,14,...}, plus a '1' if odd
    target=n if n%2==0 else n-1
    parts2mod4=[k for k in range(2,n+1) if k%4==2]
    def parts_into(t,allowed,mx=None):
        if t==0: yield (); return
        for k in [a for a in allowed if a<=t and (mx is None or a<=mx)]:
            for rest in parts_into(t-k,allowed,k): yield (k,)+rest
    pred=set()
    for pp in parts_into(target,parts2mod4):
        q=tuple(sorted(pp,reverse=True))+((1,) if n%2==1 else ())
        pred.add(tuple(sorted(q,reverse=True)))
    print(f"  n={n:>2}: #contrib={len(ct)}  matches 'parts=2mod4 (+1 if odd)'? {ct==pred}   types={sorted(ct,reverse=True)[:4]}")
print("\n=> the VALID anti-automorphisms = permutations whose cycles are ALL length =2 mod 4 (2,6,10,...),")
print("   plus exactly ONE fixed point iff n is odd. #SC(n)=Burnside over just these => the unified (even&odd) rule.")
print("   EVEN/ODD MIRROR: odd-n rule = even-n rule with one extra fixed vertex (the 1-cycle); NOT #SC(2k)=A000568(2k-1).")
