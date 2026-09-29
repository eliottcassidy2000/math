"""Census of j-fold consecutive fair splits of cyclic binary necklaces.
A necklace w of length K=jm with X=jx ones admits a j-fold consecutive fair split
iff some rotation r has every window [r+im, r+(i+1)m) containing exactly x ones.
j=2: always (discrete IVT). j>=3: not always. We count, per shape (K,X) and j | gcd(K,X),
the primitive necklaces admitting such a split, and print the failing fraction."""
import sys, math
from itertools import combinations
from math import gcd

def necklaces(K, X):
    """canonical (lexicographically minimal rotation) primitive necklaces with X ones as tuples"""
    seen=set(); out=[]
    for ones in combinations(range(K), X):
        w=[0]*K
        for p in ones: w[p]=1
        rots=[tuple(w[i:]+w[:i]) for i in range(K)]
        m=min(rots)
        if m in seen: continue
        seen.add(m)
        # primitive iff all rotations distinct
        if len(set(rots))==K: out.append(m)
    return out

def has_fair_split(w, j):
    K=len(w); m=K//j; X=sum(w); x=X//j
    ww=w+w
    for r in range(m):   # a rotation by m maps the split to itself; r in [0,m) suffices
        if all(sum(ww[r+i*m:r+(i+1)*m])==x for i in range(j)):
            return r
    return None

if __name__=='__main__':
    Kmax=int(sys.argv[1]) if len(sys.argv)>1 else 20
    print("K X j #prim_necklaces #with_j_fair_split fraction_failing")
    rows=[]
    for K in range(2,Kmax+1):
        for X in range(1,K):
            g=gcd(K,X)
            if g==1: continue
            neck=necklaces(K,X)
            for j in range(2,g+1):
                if g%j: continue
                ok=sum(1 for w in neck if has_fair_split(w,j) is not None)
                print(K,X,j,len(neck),ok,"%.4f"%(1-ok/len(neck)) if neck else "-")
                rows.append((K,X,j,len(neck),ok))
                if j==2: assert ok==len(neck), ("IVT violated?!",K,X)
    print("IVT check passed: every primitive necklace with K,X even has a 2-fold fair split.")
    # sequence: for j=3, shapes (3m, 3x): counts of necklaces with a 3-fold fair split
    print("\nj=3 table: (K,X) prim  with3split")
    for (K,X,j,n,ok) in rows:
        if j==3: print((K,X),n,ok)

# ---- criteria cross-checks (run as: python3 fairsplit_census.py K criteria)
def gaps_of(w):
    K=len(w); ones=[i for i in range(K) if w[i]]
    return [ (ones[(t+1)%len(ones)]-ones[t])%K or K for t in range(len(ones)) ]
def gap_criterion(w,j):
    """X=j beads: split exists iff every run of t consecutive cyclic gaps (1<=t<=j-1) has sum in ((t-1)m,(t+1)m)"""
    K=len(w); m=K//j; g=gaps_of(w); assert len(g)==j
    gg=g+g
    for t in range(1,j):
        for s in range(j):
            S=sum(gg[s:s+t])
            if not ((t-1)*m < S < (t+1)*m): return False
    return True
def level_criterion(w,j):
    """general X: split exists iff the discrepancy G(t)=F(t)-tX/K takes a common value on some coset r+mZ"""
    K=len(w); m=K//j; X=sum(w)
    F=[0]
    for b in w: F.append(F[-1]+b)
    # G(t) scaled by K to stay integral: K*F(t)-t*X
    for r in range(m):
        vals={K*F[(r+i*m)]-(r+i*m)*X for i in range(j)} if r+ (j-1)*m<=K else None
        # positions r+im for i=0..j-1 are all < K when r<m; use F on [0,K]
        if len(vals)==1: return True
    return False
if __name__=='__main__' and len(sys.argv)>2 and sys.argv[2]=='criteria':
    Kmax=int(sys.argv[1]); bad=0; tested=0
    for K in range(2,Kmax+1):
        for X in range(1,K):
            g=gcd(K,X)
            for j in range(2,g+1):
                if g%j: continue
                for w in necklaces(K,X):
                    bf=has_fair_split(list(w),j) is not None
                    lv=level_criterion(list(w),j)
                    if bf!=lv: bad+=1; print("LEVEL MISMATCH",K,X,j,w)
                    if X==j:
                        gc=gap_criterion(list(w),j)
                        if bf!=gc: bad+=1; print("GAP MISMATCH",K,X,j,w,gaps_of(list(w)))
                    tested+=1
    print(f"criteria cross-check: {tested} (necklace,j) pairs, mismatches={bad}")
    # closed count for X=j=3: primitive necklaces with a 3-fold fair split = m(m-1)
    for m in range(2,8):
        K=3*m; cnt=sum(1 for w in necklaces(K,3) if has_fair_split(list(w),3) is not None)
        print(f"  (3m,3) m={m}: with-split={cnt} = m(m-1)={m*(m-1)}; primitive total={len(necklaces(K,3))} = 3C(m,2)={3*m*(m-1)//2}")
