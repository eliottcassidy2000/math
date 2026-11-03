#!/usr/bin/env python3
"""Exact controls for complementary ears; no asserts, no external packages."""
from functools import lru_cache

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def add(a,mask):
    n=len(a)
    return tuple(a[i]|((0 if mask>>i&1 else 1)<<n) for i in range(n))+(mask,)
def left(a,q): return add(a,((1<<len(a))-1)^(1<<q))
def right(a,s): return add(a,1<<s)
def delete(a,q):
    verts=[i for i in range(len(a)) if i!=q]
    return tuple(sum(1<<j for j,v in enumerate(verts) if a[u]>>v&1) for u in verts)
def reverse(a): return tuple((((1<<len(a))-1)^(1<<i)^a[i]) for i in range(len(a)))
@lru_cache(None)
def endings(a):
    n=len(a)
    if n==0:return ()
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for i in range(n):
            z=dp[mask][i]
            if not z:continue
            avail=a[i]&~mask
            while avail:
                b=avail&-avail;avail-=b;j=b.bit_length()-1
                dp[mask|b][j]+=z
    return tuple(dp[-1])
def H(a):return sum(endings(a)) if a else 1
def starts(a):return endings(reverse(a))
def strong(a):
    n=len(a)
    for i in range(n):
        seen=1<<i;old=0
        while old!=seen:
            old=seen
            for v in range(n):
                if seen>>v&1:seen|=a[v]
        if seen!=(1<<n)-1:return False
    return True
def profile(a,q):
    ins=[i for i in range(len(a)) if a[i]>>q&1]
    require(len(ins)==1,'not a guard')
    r=ins[0];b=delete(a,q);r-=r>q
    return H(a),H(a)-starts(a)[q],starts(b)[r]
def observer(a,p):
    h=H(a);b=add(a,1<<p);h1=H(b);c=add(b,1<<len(a));h2=H(c)
    return h,(h1-h)//2,(h2-h1)//2

def C(v):h,u,d=v;return h+u+d,h+d,h-u
def D(v):h,u,d=v;return h+u+d,2*u,0

def run(a,q,word):
    for ch in word:
        a=left(a,q)
        if ch=='C':q=len(a)-1
        elif ch!='D':raise RuntimeError('bad word')
    return a,q

def labelled(n,code):
    a=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
            k+=1
    return tuple(a)

def protected(a,q,s):
    n=len(a);require(q!=s,'ports equal')
    require(sum(a[i]>>q&1 for i in range(n))==1,'left guard invalid')
    require(a[s].bit_count()==1,'right guard invalid')
    L=[b for b in range(n) if b!=q and all(b==i or i==q or a[b]>>i&1 for i in range(n))]
    R=[c for c in range(n) if c!=s and all(c==i or i==s or a[i]>>c&1 for i in range(n))]
    require(any(len({q,s,b,c})==4 for b in L for c in R),'disjoint protected pairs absent')

def permute(a,order):
    return tuple(sum(1<<j for j,v in enumerate(order) if a[u]>>v&1) for u in order)

def main():
    cases=0
    for n in [3,4]:
        for code in range(1<<(n*(n-1)//2)):
            a=labelled(n,code)
            if not strong(a):continue
            for q in range(n):
                if sum(a[i]>>q&1 for i in range(n))!=1:continue
                v=profile(a,q);b=left(a,q)
                require(strong(b),'left attachment loses strongness')
                require(profile(b,n)==C(v),'C formula')
                require(profile(b,q)==D(v),'D formula')
                for w in ['CDDC','CDDD']:
                    c,p=run(a,q,w);h,u,d=v
                    expect=(9*h+8*d,5*h+4*d,h) if w=='CDDC' else (9*h+8*d,8*h+8*d,0)
                    require(profile(c,p)==expect,'complementary double-offset profile')
                    require(strong(c),'double offset loses strongness')
                cases+=1
    print('All labelled strong guard cases, orders 3 and 4:',cases)
    examples=[((8,1,3,6),2,3,1),((20,1,2,7,14),3,4,2)]
    for a,q,p,j in examples:
        v=profile(a,q);require(v[2]==j,'named initial defect')
        b,r=run(a,q,'CD');c,t=run(b,r,'DC')
        require([H(a),H(b),H(c)]==[H(a),3*H(a)+2*j,9*H(a)+8*j],'offset chain')
        require(observer(c,p)==(9*observer(a,p)[0]+8*j,9*observer(a,p)[1],9*observer(a,p)[2]),'named fixed observer')
        print('Offset chain:',H(a),H(b),H(c),'final guard profile',profile(c,t),'fixed observer',observer(c,p))
        print('Rows:',c,'guard',t,'old observer',p)
    seeds=[((8,1,19,6,43,15),4,0,(17,12,12,8)),((20,1,2,39,14,87,31),5,1,(31,22,22,16))]
    tested=0
    for a,q,s,coef in seeds:
        protected(a,q,s);h,u,v,w=coef
        require(H(a)==h,'corner00')
        require(H(left(a,q))==h+u,'corner10')
        require(H(right(a,s))==h+v,'corner01')
        require(H(right(left(a,q),s))==h+u+v+w,'corner11')
        # The two insertion orders differ only by the two new labels.
        lr=right(left(a,q),s);rl=left(right(a,s),q)
        require(lr==permute(rl,list(range(len(a)))+[len(a)+1,len(a)]),'left/right commute')
        for k in range(3):
            for ell in range(3):
                b=a
                for _ in range(k):b=left(b,q)
                for _ in range(ell):b=right(b,s)
                protected(b,q,s)
                require(strong(b),'two-sided strongness')
                A=(1<<k)-1;B=(1<<ell)-1
                require(H(b)==h+A*u+B*v+A*B*w,'bivariate interpolation')
                require(profile(b,q)[2]==0,'left zero defect')
                require(profile(reverse(b),s)[2]==0,'right zero defect')
                tested+=1
        print('Two-sided corners:',coef,'h(2,0)=',h+3*u,'h(2,2)=',h+3*u+3*v+9*w)
    print('Two-sided finite controls:',tested)
    print('PASS')
if __name__=='__main__':main()
