#!/usr/bin/env python3
from itertools import product
from test import sig

def rows(n,code):
    a=[0]*n;jj=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>jj&1:a[i]|=1<<j
            else:a[j]|=1<<i
            jj+=1
    return a

def trans(n):
    return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]

def graph(L,l,C,R,r):
    a=[];n=len(l)+len(C)+len(r)+1
    marks=l+C+r
    for B in (L,trans(len(C)),R):
        start=len(a);end=start+len(B)
        later=((1<<(n-1))-1)^((1<<end)-1)
        a.extend((v<<start)|later for v in B)
    for i,b in enumerate(marks):
        if not b:a[i]|=1<<(n-1)
    a.append(sum(int(b)<<i for i,b in enumerate(marks)))
    return a

def paths(a):
    n=len(a);f=(1<<n)-1
    d=[[0]*n for i in range(f+1)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,f):
        for i,h in enumerate(d[s]):
            if h==0:continue
            t=a[i]&~s
            while t:
                z=t&-t;t-=z;d[s|z][z.bit_length()-1]+=h
    return sum(d[-1])

def w(s):return [int(x) for x in s]
centre=list(map(w,("10100101","10011001","10000001")))
up=[]
for m in range(3,15):
    steps=set()
    for i in range(1<<(m-2)):
        a=i<<1;b=a|1|(1<<(m-1))
        s,e,I=sig(a,m);ss,ee,J=sig(b,m)
        D=J-I
        assert (s,e)==(ss,ee)
        if m>3:
            L=m-2;z=format(i,f'0{L}b')
            assert 2*D==(1<<L)-i-int(z[::-1],2)-int(z[0])-int(z[-1])
        if D==1:up.append(format(a,f'0{m}b'))
        if D>0:steps.add(D)
    assert min(steps)==(1 if m in (3,5) else 2)
assert up==["000","00010","01000"]
for i,s in enumerate(centre):
    assert sig(int(''.join(map(str,s)),2),8)==(0,0,30+i)

for L,marks in [(trans(1),[1]),(trans(2),[1,0]),(rows(3,2),[1,0,1]),(rows(4,17),[1,1,0,0])]:
    k=len(L);g=paths(L);R=[0]*k
    for i,a in enumerate(L):
        for j in range(k):
            if a>>j&1:R[k-1-j]|=1<<(k-1-i)
    rm=[1-x for x in marks[::-1]]
    N=[paths(graph(L,marks,C,R,rm)) for C in centre]
    assert N[2]-N[1]==N[1]-N[0]==2*g*g
    if k==1:assert N==[573,575,577]
    print('matched',g,N)
# Independent exterior: mismatched transitive masks yield the predicted shear.
L,R=trans(2),trans(2)
N=[paths(graph(L,[1,1],C,R,[0,1])) for C in centre]
# S_L=3, reverse-complement end address T_R=1.
assert (N[1]-N[0],N[2]-N[1])==(2+12*2,2+24*2)
print('mismatch',N)
print('unit states',up)
