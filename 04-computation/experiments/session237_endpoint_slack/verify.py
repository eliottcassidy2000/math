#!/usr/bin/env python3
"""Independent subset DP controls for the endpoint slack bridge."""
import argparse, json, itertools

def require(x,msg):
    if not x: raise RuntimeError(msg)
def tournament(n,code):
    r=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            a,b=(i,j) if code>>k&1 else (j,i);r[a]|=1<<b;k+=1
    return r
def ends(r):
    n=len(r);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,1<<n):
        for i,x in enumerate(d[s]):
            if x:
                a=r[i]&~s
                while a:
                    b=a&-a;a-=b;j=b.bit_length()-1;d[s|b][j]+=x
    return d[-1]
def sub(r,vs):return [sum(1<<j for j,b in enumerate(vs) if r[a]>>b&1) for a in vs]
def transitive(r):return sorted(x.bit_count() for x in r)==list(range(len(r)))
def equality_shape(r,a):
    # Intrinsic recognition, without factorial isomorphism tests.
    n=len(r);reach=[r[i]|1<<i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    C=[v for v in range(n) if reach[a]>>v&1 and reach[v]>>a&1]
    P=[v for v in range(n) if v not in C]
    if len(C)<3 or not transitive(sub(r,P)):return False
    if any(not (r[p]>>v&1) for p in P for v in C):return False
    O=[v for v in C if r[a]>>v&1]
    if len(O)!=1:return False
    b=O[0];chain=[v for v in C if v!=b]
    if not transitive(sub(r,chain)):return False
    # The only b-outneighbor in this chain must be its source.
    order=sorted(chain,key=lambda v:-(r[v]&sum(1<<u for u in chain)).bit_count())
    return [v for v in chain if r[b]>>v&1]==[order[0]] and order[-1]==a

def apex(l,mask):
    a=l;b=l+1;r=[0]*(l+2)
    for i in range(l+1):
        for j in range(i+1,l+1):r[i]|=1<<j
    for i in range(l):
        if mask>>i&1:r[b]|=1<<i
        else:r[i]|=1<<b
    r[a]|=1<<b
    return r,a

def energy(pos):
    return sum(2 if j-i==1 else 1<<(j-i-1) for i,j in itertools.combinations(pos,2))
def main():
    p=argparse.ArgumentParser();p.add_argument('--fail',action='store_true');a=p.parse_args()
    require(not a.fail,'deliberate checker failure')
    counts=[];eq=0;roots=0
    for n in range(1,7):
        z=0
        for code in range(1<<(n*(n-1)//2)):
            r=tournament(n,code);e=ends(r);h=sum(e)
            for v in range(n):
                if not r[v] or not e[v]:continue
                roots+=1; require(2*e[v]+1<=h,'half bound')
                same=(h==2*e[v]+1)
                require(same==equality_shape(r,v),f'equality shape {n},{code},{v}')
                if same:require(e[v]&(e[v]-1)==0,'power two');z+=1
        counts.append([n,z]);eq+=z
    apex_count=0
    for l in range(0,10):
        for mask in range(1<<l):
            r,v=apex(l,mask);e=ends(r);h=sum(e);pos=[i for i in range(l) if mask>>i&1]
            q=energy(pos);m=sum(1<<(l-1-i) for i in pos)
            require(h-2*e[v]==1+q,'slack energy')
            require(e[v]==m-q,'endpoint energy')
            require(h==1+2*m-q,'root sequence')
            require(q==0 if len(pos)<=1 else q>0,'Q zero')
            apex_count+=1
    print(json.dumps({'status':'PASS','rooted_nonsink_positive_cases':roots,'equality_by_order':counts,'equality_total':eq,'all_apex_words_through_length9':apex_count,'universe':'all labelled tournaments n<=6; all apex masks l<=9'},sort_keys=True,indent=2))
if __name__=='__main__':main()
