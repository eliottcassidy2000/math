"""Finite exact proofs of three-port perfect-number rails; stdlib only."""
from collections import Counter, defaultdict
from itertools import product
from math import prod

LOWER=(2,8,3,37,15,151,63,95)
LPORTS=(6,4,1)
UPPER=(2,16,11,19,69,31,47,127)
UPORTS=(6,2,0)


def strong(rows):
    allv=(1<<len(rows))-1
    for i in range(len(rows)):
        z=1<<i
        while True:
            a=z
            for j in range(len(rows)):
                if z>>j&1:a|=rows[j]
            if a==z:break
            z=a
        if z!=allv:return False
    return True


def cap_two(rows,marks):
    for a in marks:
        E=[(o,i) for o in range(len(rows)) if rows[a]>>o&1
                 for i in range(len(rows)) if (rows[i]>>a&1) and (rows[o]>>i&1)]
        assert all(not ((u!=v or u in marks) and (x!=y or x in marks))
                   for u,x in E for v,y in E)


def signature(rows,marks):
    cap_two(rows,marks)
    n=len(rows)
    U=((1<<n)-1)^sum(1<<j for j in marks)
    dp={(0,(0,)*len(marks),-1):1}
    data=Counter()
    for _ in range(n+len(marks)+1):
        m=defaultdict(int)
        for (state,counts,end),ways in dp.items():
            if state==U and all(counts):
                data[counts]+=ways
            for j in range(n):
                if end>=0 and not rows[end]>>j&1:continue
                if j in marks:
                    i=marks.index(j)
                    if counts[i]==2:continue
                    v=counts[:i]+(counts[i]+1,)+counts[i+1:]
                    m[state,v,j]+=ways
                elif not state>>j&1:
                    m[state|(1<<j),counts,j]+=ways
        dp=m
    return dict(data)


def power_coeffs(sig):
    r=[0]*8
    for s,v in sig.items():
        y=sum(1<<j for j in range(3) if s[j]==2)
        q=y
        while True:
            r[q]+=v*(-2)**((q^y).bit_count())
            if not q:break
            q=(q-1)&y
    return tuple(r)


def H(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    for j in range(n):dp[1<<j][j]=1
    for k in range(1,1<<n):
        for j,v in enumerate(dp[k]):
            if not v:continue
            z=rows[j]&~k
            while z:
                b=z&-z;t=b.bit_length()-1
                dp[k|b][t]+=v
                z^=b
    return sum(dp[-1])


def inflate(rows,ports,widths):
    values=dict(zip(ports,widths))
    seq=[i for i in range(len(rows)) for _ in range(values.get(i,1))]
    return tuple(sum(1<<j for j,w in enumerate(seq)
                  if (w!=v and rows[v]>>w&1) or (w==v and i<j))
                 for i,v in enumerate(seq))


def cube(sig,widths):
    return sum(c*prod((1<<w)-2 if k==2 else 1 for w,k in zip(widths,r))
               for r,c in sig.items())


def Fold(a):
    out=1
    while a:
        a,b=divmod(a,2)
        if not b:out+=a+(a&1)
    return out


def lower_read(u,w):return 2*u*w+u+(u&1)+Fold(u)+Fold(w)-1

def check():
    sig=signature(LOWER,LPORTS)
    assert strong(LOWER)
    print('odd lower skeleton:',sorted(sig.items()),flush=True)
    assert power_coeffs(sig)==(-1,-1,-3,-3,3,3,4,4)
    T=[(p,(p,(p-3)//2,(p-3)//2)) for p in (5,7,9,31,127)]
    for p,m in T:
        E=(1<<(p-1))*((1<<p)-1)
        assert cube(sig,m)==E-1
        if p<=7: assert H(inflate(LOWER,LPORTS,m))==E-1
    print('verified odd exponents 5,7,9,31,127',flush=True)
    for e in range(3,33):
        T=1<<(e-1)
        n=2*T*T-T-1
        actual=[lower_read(u,w) for u,w in ((T-2,T-1),(T-2,T),(T-1,T+1),(T-1,T+2))]
        assert actual==[2*T*T-9*T//2+3,2*T*T-3*T//2-1,2*T*T+3*T//2-1,2*T*T+15*T//4-1-(T==4)]
        assert n not in actual
    assert (Fold(10),Fold(12))==(9,11)
    print('one-apex obstruction rails checked through p=32',flush=True)


if __name__=='__main__':check()
