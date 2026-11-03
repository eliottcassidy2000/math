def require(x):
 if not x: raise AssertionError("state bank check failed")
from itertools import product, permutations
from collections import defaultdict

def automaton(word, ab):
    k=2; total=0
    def edge_a(a,b): return bool(ab) if (a,b)==(0,1) else not bool(ab)
    def xbeats(x,a): return bool(x & (1<<a))
    for order in permutations(range(k)):
        d={(0,)*(k+1):1}
        for x in word:
            z=defaultdict(int)
            for state,n in d.items():
                for j in range(k+1):
                    left=order[j-1] if j else None
                    right=order[j] if j<k else None
                    if state[j]==0 and left is not None and xbeats(x,left): continue
                    passed=right is None or xbeats(x,right)
                    s=list(state); s[j]=2 if passed else 1
                    z[tuple(s)]+=n
            d=z
        for s,n in d.items():
            ok=True
            for j,t in enumerate(s):
                if t==1: ok=False
                if t==0 and j>0 and j<k and not edge_a(order[j-1],order[j]): ok=False
            if ok: total+=n
    return total

def rows(word, ab):
    m=len(word); n=m+2; A=[0]*n
    for i,x in enumerate(word):
        for j in range(i+1,m): A[i]|=1<<j
        for a in range(2):
            if x>>a&1: A[i]|=1<<(m+a)
            else: A[m+a]|=1<<i
    A[m+(0 if ab else 1)] |= 1 << (m+(1 if ab else 0))
    return A

def strong(A):
    n=len(A)
    for a in range(n):
        reach=1<<a
        while True:
            old=reach
            for i in range(n):
                if reach>>i&1: reach|=A[i]
            if old==reach: break
        if reach!=(1<<n)-1:return False
    return True

def dp(A):
    n=len(A); d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            if d[s][i]:
                for j in range(n):
                    if not (s>>j&1) and A[i]>>j&1:d[s|1<<j][j]+=d[s][i]
    return sum(d[-1])

if __name__=='__main__':
    require(automaton([1,1,3,0,3],0)==35)
    require(automaton([1,1,1,2,3],0)==49)
    for m in range(0,5):
        for word in product(range(4),repeat=m):
            for ab in range(2):
                A=rows(word,ab); a=automaton(word,ab); b=dp(A)
                if a!=b: raise AssertionError((word,ab,a,b))
    print('independent DP comparison: all 682 banks, m=0..4')
    for m in range(2,8):
        vals=set(); count=0
        for word in product(range(4),repeat=m):
            for ab in range(2):
                A=rows(word,ab)
                if strong(A):
                    count+=1; h=automaton(word,ab); vals.add(h)
                    t=max(word.count(i) for i in range(4))
                    if h<2**((t+1)//2)+1: raise AssertionError((word,ab,h,t))
        print('m=%d strong_codes=%d distinct_counts=%d max=%d small_counts=%s'%(m,count,len(vals),max(vals) if vals else 0, sorted(v for v in vals if v<=75)))

