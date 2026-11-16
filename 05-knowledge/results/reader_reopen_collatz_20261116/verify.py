from functools import lru_cache
from pathlib import Path
import json

BANK=json.loads(Path(__file__).with_name('readers.json').read_text())

def verify_tournament(a):
    for i,r in enumerate(a):
        assert 0<=r<1<<len(a) and not r>>i&1
        for j in range(i):assert (r>>j&1)+(a[j]>>i&1)==1

def profile(a,x):
    n=len(a); F=(1<<n)-1
    @lru_cache(None)
    def f(S,i):
        R=S^(1<<i)
        return sum(f(R,j) for j in range(n) if R>>j&1 and a[i]>>j&1) if R else 1
    @lru_cache(None)
    def H(S):return sum(f(S,j) for j in range(n) if S>>j&1) if S else 1
    r=s=0
    for S in range(1,F+1):
        if not S>>x&1:continue
        z=f(S,x)*H(F^S)
        if S.bit_count()%2:s+=z
        else:r+=z
    h,g=H(F),H(F^(1<<x))
    return h,g,r,s

def H(a):
    n=len(a);N=1<<n
    dp=[[0]*n for _ in range(N)]
    for j in range(n):dp[1<<j][j]=1
    for S in range(1,N-1):
        for i,k in enumerate(dp[S]):
            if k:
                O=a[i]&~S
                while O:
                    j=O&-O;O^=j
                    dp[S|j][j.bit_length()-1]+=k
    return sum(dp[-1])

def mat(labels,edge):
    return [sum(1<<j for j,v in enumerate(labels) if u!=v and edge(u,v)) for u in labels]

def motif(s,reserve):
    # Labels L0..L3=0..3,R0..R2=4..6,sentinel=7,head=8,
    # long reset spacer=9..9+s-1, newcomers=9+s,10+s.
    T=list(range(9,9+s))
    old=list(range(9+s))
    order=[0,1,2,3]+T+[4,5,6,7]
    pos={x:i for i,x in enumerate(order)}
    marks={1,3,4,5,7}
    def e0(u,v):
        return u in marks if v==8 else v not in marks if u==8 else pos[u]<pos[v]
    seed=(6,12,40,33,15,83,159,63)
    roles={u:j for u,j in enumerate((6,7,4,3,0,8,9,10,5))}
    roles[9+s],roles[10+s]=1,2
    if reserve:
        lead=T.pop(0)
        T.append(5)
        del roles[5]
    else:lead=-1
    vpos={x:i for i,x in enumerate(T)}
    def e1(u,v):
        if u==lead:return True
        if v==lead:return False
        if u in vpos and v in vpos:return vpos[u]<vpos[v]
        if u in vpos:return roles[v]==0 or roles[v]>=8
        if v in vpos:return 0<roles[u]<8
        i,j=roles[u],roles[v]
        return bool(seed[i]>>j&1) if i<8 and j<8 else i<j
    a=mat(old,e0)
    b=mat(old+[9+s,10+s],e1)
    diff=[(i,j) for i in old for j in old if i<j and (a[i]>>j&1)!=(b[i]>>j&1)]
    assert len(diff)==(9 if reserve else 3)
    if s==6:assert H(a)==121*2**s+1 and H(b)==121*2**s+3
    return diff

if __name__=='__main__':
    assert [row[0] for row in BANK]==list(range(-5,70,2))
    for B,A,x,a in BANK:
        verify_tournament(a)
        h,g,r,s=profile(a,x)
        assert (h-2*s+2*g,r+s-g)==(B,A)
        assert A//(A&-A)==121 and A<=121*16
    print('38 consecutive offsets: all seed profiles independently verified.')
    for reserve in (False,True):
        print('boundary old arc reversals:',motif(6,reserve))

def shared(s,row):
    B,A,x,seed=row
    j=(A//121).bit_length()-1
    seats=[9+s,10+s,3,2,8,0,1,9]
    roles=dict(zip(seats,[z for z in range(len(seed)) if z!=x]))
    roles[4]=x
    neutral=set(seats)-roles.keys()
    full=list(range(s+11))
    common=list(range(10,s+6))
    surplus=list(range(s+6,s+9))+[5]
    live=common+surplus[:4-j]
    terminal=surplus[4-j:]+[6,7]
    p={u:i for i,u in enumerate(neutral)}
    l={u:i for i,u in enumerate(live)}
    e={u:i for i,u in enumerate(terminal)}
    def rule(u,v):
        if u in p or v in p:
            return u in p and (v not in p or p[u]<p[v])
        if u in e or v in e:
            return v in e and (u not in e or e[u]<e[v])
        if u in l and v in l:return l[u]<l[v]
        if u in l:return roles[v]==x
        if v in l:return roles[u]!=x
        return bool(seed[roles[u]]>>roles[v]&1)
    a=mat(full,rule)
    verify_tournament(a)
    return a,common

if __name__=='__main__':
    for s in (6,7,12):
        prev=None
        for row in BANK:
            a,T=shared(s,row)
            if s==6 and row[0] in (-5,1,3,5,7,69):
                assert H(a)==121*2**s+row[0]
            if prev:
                for i in T+[6,7]:assert a[i]==prev[i]
                assert sum((x^y).bit_count() for x,y in zip(a,prev))<=156
            prev=a
    print('Fixed-tail 13-seat bank and sampled full literal counts verified.')
    A=[((1<<14)-1)^((1<<(i+1))-1) for i in range(14)]+[0]
    for i in range(14):
        if i in (1,3,10,11,13):A[i]|=1<<14
        else:A[14]|=1<<i
    assert H(A)==7745
    for i,j in ((3,11),(0,14),(1,14),(13,14)):
        A[i]^=1<<j; A[j]^=1<<i
    assert H(A)==7747
    print('First-height same-order four-flip alternative verified.')
