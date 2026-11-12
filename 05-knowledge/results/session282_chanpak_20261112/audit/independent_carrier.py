#!/usr/bin/env python3
from itertools import combinations
import json

def check(p,msg):
    if not p:raise RuntimeError(msg)
def path_tables(T):
    n=len(T);N=1<<n
    ends=[[0]*n for _ in range(N)];starts=[[0]*n for _ in range(N)]
    H=[0]*N;H[0]=1
    for S in range(1,N):
        for i in range(n):
            if not(S>>i&1):continue
            R=S^(1<<i)
            if not R:ends[S][i]=starts[S][i]=1
            else:
                ends[S][i]=sum(ends[R][j]for j in range(n)if R>>j&1 and T[j]>>i&1)
                starts[S][i]=sum(starts[R][j]for j in range(n)if R>>j&1 and T[i]>>j&1)
        H[S]=sum(ends[S]);check(H[S]==sum(starts[S]),'dual counts')
    return H,starts,ends

def packet(T):
    D=tuple((row>>1) for row in T[1:]);H,S,E=path_tables(D);n=len(D);full=(1<<n)-1
    out=[i for i in range(n)if T[0]>>(i+1)&1];inc=[i for i in range(n)if T[i+1]&1]
    U=V=d=e=H[full];Q=k=0
    total=0
    for mask in range(1<<n):
        inv=full^mask;sign=(-1)**mask.bit_count()
        Q+=H[mask]*H[inv];k+=sign*H[mask]*H[inv]
        if mask:
            x=sum(S[mask][i]for i in out)*H[inv]
            y=sum(E[mask][i]for i in inc)*H[inv]
            U+=x;V+=y;d+=sign*x;e+=sign*y
        before=1 if not mask else sum(E[mask][i]for i in inc)
        after=1 if not inv else sum(S[inv][i]for i in out)
        total+=before*after
    direct=path_tables(T)[0][-1]
    check(total==direct,'root splitting vs direct')
    return dict(H=total,h=H[full],Q=Q,U=U,V=V,d=d,e=e,k=k)

if __name__ == '__main__':
    X=(594,2020,2001,1991,1994,2013,1152,513,193,320,897)
    Y=(1542,2020,2000,1991,1995,2013,1153,513,193,320,896)
    px,py=packet(X),packet(Y)
    check(px['H']==2747 and py['H']==2749,'new scalar unit gap')
    expected=dict(h=169,Q=10816,U=4030,V=6955,d=104,e=65,k=0)
    check(all(px[k]==py[k]==v for k,v in expected.items()),'new full central packet')
    A=(56,13,25,48,34,6);B=(52,13,24,49,34,6)
    pa,pb=packet(A),packet(B)
    check(pb['H']-pa['H']==2 and all(pa[k]==pb[k]for k in pa if k!='H'),'old central control')
    print(json.dumps(dict(status='PASS',new_packets=[px,py],old_packets=[pa,pb],method='independent deletion start/end tables and root splitting; direct full path table control'),indent=2,sort_keys=True))
