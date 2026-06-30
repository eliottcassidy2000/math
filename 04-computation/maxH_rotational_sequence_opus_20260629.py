"""
(1) explicit MAXIMIZING circulant connection sets (look for structure).
(2) a fast TECHNIQUE: do structured circulant families admit a linear recurrence in n?
    Test the 'rotational' R_n (i beats i+1..i+(n-1)/2) and 'consecutive' families.
(3) more new sequences: #distinct H over circulants; the maximizer's H growth ratio.
"""
from itertools import combinations

def H_circ(n, outset):
    A=[[((j-i)%n) in outset for j in range(n)] for i in range(n)]
    dp=[[0]*n for _ in range(1<<n)]; dp[1][0]=1
    for mask in range(1<<n):
        if not mask&1: continue
        row=dp[mask]
        for last in range(n):
            c=row[last]
            if not c or not (mask>>last)&1: continue
            al=A[last]
            for nxt in range(n):
                if (mask>>nxt)&1 or not al[nxt]: continue
                dp[mask|1<<nxt][nxt]+=c
    return n*sum(dp[(1<<n)-1])

def best_circ(n):
    half=(n-1)//2; best=-1; bestout=None; allH=set()
    for r in range(half+1):
        for S in combinations(range(1,half+1),r):
            out=set(k if k in S else n-k for k in range(1,half+1))
            h=H_circ(n,out); allH.add(h)
            if h>best: best=h; bestout=set(out)
    return best,bestout,allH

print("(1) maximizing circulant connection sets (the 'out'-differences):")
for n in [3,5,7,9,11,13,15]:
    b,o,allH=best_circ(n)
    qr=set((x*x)%n for x in range(1,n)) if True else None
    print(f"  n={n:>2}: maxH={b:>11}  out-set={sorted(o)}  (QR mod n={sorted(qr & set(range(1,n)))[:6]}...) #distinctH={len(allH)}")

print("\n(2) rotational R_n (i beats next (n-1)/2): H sequence + recurrence search")
Rseq=[]
for n in range(3,22,2):
    out=set(range(1,(n-1)//2+1))
    Rseq.append((n,H_circ(n,out)))
for n,h in Rseq: print(f"  R_{n}: H={h}")
vals=[h for _,h in Rseq]
# simple ratio test (geometric-ish?) and check OEIS-friendly differences
print("  ratios H(R_{n+2})/H(R_n):", [round(vals[i+1]/vals[i],3) for i in range(len(vals)-1)])

# (3) growth of the MAX-H maximizer
print("\n(3) max-H growth and the Alon ratio H_max(n)*2^{n-1}/n! (Alon: bounded by ~c n^{3/2}):")
import math
A={3:3,5:15,7:189,9:3357,11:95095,13:3711175,15:198464295}
for n,h in A.items():
    ratio=h*2**(n-1)/math.factorial(n)
    print(f"  n={n:>2}: H_max={h:>11}  H_max*2^(n-1)/n! = {ratio:.4f}  (/n^1.5={ratio/n**1.5:.5f})")
