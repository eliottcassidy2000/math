"""Finite controls for the independent quotient and actual marked starts."""
from functools import lru_cache
from itertools import permutations

IR=(4,1,10,67,15,31,55)
WR=(8,1,3,6)
IC=81
WC=13

def ck(ok, receipt):
    if not ok: raise AssertionError(receipt)

def inflate(word):
    m=sum(word)
    u=set(range(7,7+m))
    w=set(range(7+m,11+m))
    v=11+m
    high=(1<<(v+1))-1
    rows=[0]*(v+1)
    for i,r in enumerate(IR):
        rows[i]=r|sum(1<<j for j in u|w)|(((IC>>i&1)^1)<<v)
    for i,r in enumerate(WR):
        rows[7+m+i]=(r<<(7+m))|(((WC>>i&1)^1)<<v)
    rows[v]=IC|(WC<<(7+m))
    p=7
    for size in word:
        later=high ^ ((1<<(p+size))-1)
        for j in range(size):
            rows[p+j]=later
            if size==3:
                rows[p+j]|=1<<(p+(j+1)%3)
        p+=size
    return tuple(rows)

def root_modules(A,x):
    V=(1<<len(A))-1
    return [m for m in range(1,V) if m>>x&1 and m.bit_count()>1 and
            all((m>>i&1) or (A[i]&m) in (0,m) for i in range(len(A)))]

def Hstart(A,r):
    @lru_cache(None)
    def count(u,R):
        if not R:return 1
        allowed=R&A[u]
        total=0
        while allowed:
            b=allowed&-allowed
            allowed^=b
            total+=count(b.bit_length()-1,R^b)
        return total
    return count(r,((1<<len(A))-1)^(1<<r))

def alpha0():
    return sum(1 for p in permutations(range(7)) if IC>>p[0]&1 and
          all(IR[x]>>y&1 for x,y in zip(p,p[1:])))

def beta():
    ans=0
    for p in permutations(range(4)):
        for j in range(5):
            if j and WC>>p[j-1]&1: continue
            if all(WR[x]>>y&1 for q in (p[:j],p[j:]) for x,y in zip(q,q[1:])):
                ans+=1
    return ans

def run():
    ck(not root_modules(inflate((1,)),7),"root module certificate")
    ck((alpha0(),beta())==(12,14),"outer port counts")
    total=0
    for word in ((1,),(1,1),(1,3),(3,1),(1,3,1),(3,1,1)):
        p=7
        for i,size in enumerate(word):
            if size==1:
                a=word.count(3)
                j=word[i+1:]
                e=2*j.count(3)+j.count(1)
                ck(Hstart(inflate(word),p)==168*3**a*2**e,(word,i,p))
                total+=1
            p+=size
    pairs=[(a,s,7*3**a*(2**s+1)) for a in range(35) for s in range(2,85)]
    ck(len({z for a,s,z in pairs})==len(pairs),"finite arithmetic control")
    print("root module candidates exhausted; outer port constants 12 and 14")
    print("literal marked starts:",total)
    print("finite lattice controls:",len(pairs))

if __name__=='__main__':run()
