"""
Improve A038375 (max Hamiltonian paths per tournament size) computation, and find new
sequences, using the H-atom recursion (max is a single strongly-connected prime, conj.
circulant). Search only the 2^{(n-1)/2} CIRCULANT tournaments instead of 2^{C(n,2)}.

Circulant tournament on Z_n (n odd): pick S subset of {1..(n-1)/2}, one of each {k,n-k};
i->j iff (j-i mod n) in connectset, where connectset = S union {n-k : k in {1..(n-1)/2}\S}.
"""
from itertools import combinations

def H_circulant(n, S):
    # connectset: for each k in 1..(n-1)/2, k in S => k is 'out', else n-k is 'out'
    out=set()
    for k in range(1,(n-1)//2+1):
        out.add(k if k in S else n-k)
    A=[[ ((j-i)%n) in out for j in range(n)] for i in range(n)]
    # H via DP, using vertex-transitivity: H = n * (#paths from 0)  -- but compute full to be safe
    full=(1<<n)-1
    # paths from fixed start 0, times n (vertex transitive)
    dp=[[0]*n for _ in range(1<<n)]
    dp[1][0]=1  # start at 0
    for mask in range(1<<n):
        if not (mask&1): continue
        row=dp[mask]
        for last in range(n):
            c=row[last]
            if not c or not (mask>>last)&1: continue
            al=A[last]
            for nxt in range(n):
                if (mask>>nxt)&1 or not al[nxt]: continue
                dp[mask|1<<nxt][nxt]+=c
    return n*sum(dp[full])

def QR_set(n):
    return set((x*x)%n for x in range(1,n) if (x*x)%n!=0) & set(range(1,(n-1)//2+1)) or \
           {k for k in range(1,(n-1)//2+1) if k in set((x*x)%n for x in range(1,n))}

A038375={1:1,2:1,3:3,4:5,5:15,6:45,7:189,8:661,9:3357}

print("Extend A038375 via circulant search (odd n). max H over 2^{(n-1)/2} circulants:")
print(f"{'n':>3} {'#circ':>7} {'maxH(circ)':>12} {'A038375':>9} {'argmax S':>20} {'QR achieves?':>12}")
new_minSC={}
for n in [3,5,7,9,11,13,15]:
    half=(n-1)//2
    best=-1; bestS=None; allH=set()
    for r in range(half+1):
        for S in combinations(range(1,half+1), r):
            h=H_circulant(n,set(S)); allH.add(h)
            if h>best: best=h; bestS=set(S)
    known=A038375.get(n,'?')
    # is QR maximal? compute QR connection set (residues in 1..(n-1)/2)
    qr={k for k in range(1,half+1) if k in set((x*x)%n for x in range(1,n))}
    hqr=H_circulant(n,qr)
    print(f"{n:>3} {2**half:>7} {best:>12} {str(known):>9} {str(sorted(bestS)):>20} {('YES' if hqr==best else f'no({hqr})'):>12}")
    new_minSC[n]=min(allH)

print("\nNEW sequence candidate -- MIN H over circulant tournaments (smallest 'rotational core'):")
for n in [3,5,7,9,11,13,15]: print(f"  n={n}: minH(circ)={new_minSC[n]}")
