#!/usr/bin/env python3
from itertools import combinations
from collections import defaultdict

def rows(n, code):
    a=[set() for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:a[i].add(j)
        else:a[j].add(i)
    return a

def ends(a):
    n=len(a);f=(1<<n)-1
    masks=[sum(1<<j for j in r) for r in a]
    d=[defaultdict(int)for _ in range(f+1)]
    for j in range(n):d[1<<j][j]=1
    for s in range(1,f+1):
        for i,w in d[s].items():
            u=masks[i]&~s
            while u:
                z=u&-u;u-=z
                d[s|z][z.bit_length()-1]+=w
    return d

def subsets(s):
    orig=s
    while True:
        yield s
        if not s:break
        s=(s-1)&orig

def H_Q(a):
    d=ends(a)
    h=[1]+[sum(v.values())for v in d[1:]]
    f=len(h)-1
    return h[f],sum(h[s]*h[f^s]for s in range(f+1))

def markstats(a,O):
    e=ends(a)
    s=ends([set(range(len(a)))-{i}-r for i,r in enumerate(a)])
    h=[1]+[sum(v.values())for v in e[1:]]
    mf=[sum(v[i]for i in O)for v in s]
    me=[sum(v[i]for i in O)for v in e]
    r=[];m=[]
    for mask in range(len(h)):
        r.append(h[mask]-sum(me[x]*mf[mask^x]for x in subsets(mask)))
        m.append(sum(h[x]*mf[mask^x]for x in subsets(mask)))
    f=len(h)-1
    return m[f],sum(h[x]*r[f^x]for x in subsets(f)),sum(h[x]*m[f^x]for x in subsets(f)),r[f],sum(h[x]*me[f^x]for x in subsets(f))

def tt(n):return [set(range(i+1,n)) for i in range(n)]

def close(a,O,b):
    p,q=len(a),len(b)
    return [u|set(range(p,p+q))|({p+q} if i not in O else set())for i,u in enumerate(a)]+[{j+p for j in u}|{p+q}for u in b]+[set(O)]

def outer(a,alpha=6):
    # 3 Q−H is the closer with initial pair TT2, both marked.
    # Optional terminal sink in module multiplies only Q by 2.
    if alpha==6:a=[r|{len(a)} for r in a]+[set()]
    assert alpha in (3,6)
    return close(tt(2),{0,1},a)

def expect_test(n, w, expected, alpha=6):
    code,mask=w
    a=rows(n,code);O={i for i in range(n)if mask>>i&1}
    assert markstats(a,O)==expected
    M,R,U,r,F=expected
    for t in range(1,4):
        b=close(a,O,tt(t));h,q=H_Q(b)
        assert (h,q)==(r+M*2**t,2*R*2**t+2*U*3**t-2*F)
        N=H_Q(outer(b,alpha))[0]
        assert N==2*alpha*U*3**t+(2*alpha*R-M)*2**t-r-2*alpha*F
    return a,O

if __name__=='__main__':
    expect_test(7,(1382658,38),(397,464,2705,-17,433))
    expect_test(7,(67680,76),(397,464,2705,-11,433))
    expect_test(8,(71836161,45),(727,-1558,7750,-313,1310))
    expect_test(8,(143663362,45),(751,-1556,7750,-311,1322))
    expect_test(8,(8388864,65),(81,96,1674,-27,306),alpha=3)
    expect_test(7,(1025,73),(177,112,1674,-25,306),alpha=3)
    print('Literal replay verified two same-response order-7 hosts and two would-be order-8 hosts (whose finite difference is 146, not 2) and alpha=3 adjacent pair at three heights.')
