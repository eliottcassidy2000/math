#!/usr/bin/env python3
"""Standalone exact Session282 template/Jordan controls; assertions not used."""
from itertools import product,permutations,combinations
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction as F
from math import comb
import json,sys

def ck(p,s):
    if not p:raise RuntimeError(s)
def tournaments(n):
    es=list(combinations(range(n),2))
    for mask in range(1<<len(es)):
        T=[0]*n
        for k,(i,j)in enumerate(es):
            if mask>>k&1:T[i]|=1<<j
            else:T[j]|=1<<i
        yield tuple(T)
@lru_cache(None)
def hpall(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)];H=[1]+[0]*((1<<n)-1)
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for j in range(n):
            if mask>>j&1:
                old=mask^(1<<j)
                if old:dp[mask][j]=sum(dp[old][i]for i in range(n)if old>>i&1 and T[i]>>j&1)
        H[mask]=sum(dp[mask])
    return tuple(H)
def hp(T):return hpall(T)[-1]
def blowup(Q,ms):
    st=[0]
    for m in ms:st.append(st[-1]+m)
    T=[0]*st[-1]
    for i,m in enumerate(ms):
        for u in range(st[i],st[i+1]):
            for v in range(u+1,st[i+1]):T[u]|=1<<v
            for j in range(len(Q)):
                if Q[i]>>j&1:
                    for v in range(st[j],st[j+1]):T[u]|=1<<v
    return tuple(T)
@lru_cache(None)
def onto(m,r):return sum((-1)**(r-j)*comb(r,j)*j**m for j in range(r+1))
def compile_runs(Q,ms):
    n=len(Q);layer={(tuple(int(j==i)for j in range(n)),i):1 for i in range(n)};coef=defaultdict(int)
    for length in range(1,sum(ms)+1):
        nxt=defaultdict(int)
        for (r,i),v in layer.items():
            if all(r):coef[r]+=v
            for j in range(n):
                if Q[i]>>j&1 and r[j]<ms[j]:
                    rr=list(r);rr[j]+=1;nxt[tuple(rr),j]+=v
        layer=nxt
    H=0
    for r,v in coef.items():
        for m,k in zip(ms,r):v*=onto(m,k)
        H+=v
    return H,dict(coef)
def cycles(Q):
    if len(Q)<3:return 0
    return sum(all(Q[w[i]]>>w[(i+1)%len(w)]&1 for i in range(len(w)))for p in permutations(range(1,len(Q)))for w in[(0,)+p])
def tape(w):
    T=[0]*(len(w)+1)
    for i,b in enumerate(w,1):
        if b:T[0]|=1<<i
        else:T[i]|=1
        for j in range(i+1,len(w)+1):T[i]|=1<<j
    return tuple(T)
def slots(w):
    if not w:return 1
    return int(w[-1]==0)+int(w[0]==1)+sum(w[i:i+2]==(0,1)for i in range(len(w)-1))+sum(2**(i-j-1)for j in range(len(w))for i in range(j+1,len(w))if w[j]==1 and w[i]==0)
def join(*Ts):
    out=[];off=0;total=sum(map(len,Ts))
    for T in Ts:
        future=((1<<total)-1)^((1<<(off+len(T)))-1)
        out.extend((row<<off)|future for row in T);off+=len(T)
    return tuple(out)
def tt(n):return tuple(sum(1<<j for j in range(i+1,n))for i in range(n))
def converse(T):return tuple(sum(1<<j for j in range(len(T))if T[j]>>i&1)for i in range(len(T)))
def rooted(T):
    n=len(T);H=hpall(T);full=(1<<n)-1;dp=[[0]*n for _ in range(1<<n)];dp[1][0]=1;P=[0]*(n+1)
    for mask in range(1,1<<n):
        P[mask.bit_count()]+=sum(dp[mask])*H[full^mask]
        for v in range(n):
            if dp[mask][v]:
                rem=T[v]&~mask
                while rem:
                    bit=rem&-rem;rem-=bit;dp[mask|bit][bit.bit_length()-1]+=dp[mask][v]
    return P
def profile(T):
    H=hpall(T);R=((1<<len(T))-1)^1;U=rooted(T);V=rooted(converse(T));Q=k=0;S=R
    while True:
        v=H[S]*H[R^S];Q+=v;k+=(-1)**S.bit_count()*v
        if not S:break
        S=(S-1)&R
    return dict(H=H[-1],h=H[R],Q=Q,k=k,U=sum(U),V=sum(V),d=sum((-1)**(j+1)*v for j,v in enumerate(U)),e=sum((-1)**(j+1)*v for j,v in enumerate(V)))

templates=floors=0
for n in(2,3,4):
    for Q in tournaments(n):
        h=hp(Q);c=cycles(Q)
        for ms in product(range(1,4 if n<4 else 3),repeat=n):
            if sum(ms)>8:continue
            actual=hp(blowup(Q,ms));compiled,_=compile_runs(Q,ms)
            ck(actual==compiled,'run compiler literal')
            if c:
                ck(actual>=h+c*sum(2**m-2 for m in ms),'cycle floor')
                for m in ms:ck(actual>=2**m+1,'module height')
                floors+=1
            templates+=1
C3=(2,4,1)
for m in range(1,11):ck(hp(blowup(C3,(m,1,1)))==2**m+1,'sharp C3 floor')
ck(compile_runs(C3,(3,1,1))[1]=={(1,1,1):3,(2,1,1):1},'fixed run tensor')
for m in range(1,30):ck(3*onto(m,1)+onto(m,2)-2**m==1,'allowed constant cancellation')
tapes=0
for n in range(8):
    for w in product((0,1),repeat=n):ck(slots(w)==hp(tape(w)),'tape literal');tapes+=1
for a,b,m in product(range(1,6),range(1,6),range(1,10)):
    w=(0,)*a+(1,)*b;R=2**(a+b);D=(2**a-1)*(2**b-1)
    value=F(D*(R**m-1),(R-1)**2)+(1-F(D,R-1))*m
    ck(value.denominator==1 and value==slots(w*m),'block tape formula')
for s in range(1,9):
    for m in range(1,15):ck((2**s+1)**2*slots(((0,)*s+(1,)*s)*m)==4**(s*m)-1+2*(2**s+1)*m,'balanced formula')
for m in range(1,20):
    ck(9*slots((0,1)*2*m)-25*slots((0,0,1,1)*m)==2*m,'linear extractor')
    ck(slots((0,1)*3*m)==9*slots((0,0,0,1,1,1)*m),'scalar twins')
ck(9*slots(())-25*slots(())==-16,'empty hostile')
profiles=[]
for s,m in((1,1),(1,3),(2,1),(2,2),(3,1)):
    w=((0,)*s+(1,)*s)*m;R=4**(s*m);uv=1+(R-1)//(2**s+1);got=profile(tape(w))
    ck(got==dict(H=slots(w),h=1,Q=R,k=0,U=uv,V=uv,d=0,e=0),'balanced carrier');profiles.append(got)
P=tape((0,1)*3);Q=tape((0,0,0,1,1,1));L=join((0,),tt(6),P);R=join((0,),C3,C3,Q)
left=profile(L);right=profile(R)
ck(len(L)==len(R)==14,'same order')
for k in('H','h','V','d','e'):ck(left[k]==right[k],'receiver equality')
ck(left['H']==9 and left['Q']==left['U']==21504 and right['Q']==right['U']==18432,'sending contrast')
ck(left['k']==right['k']==0,'nontrivial boundary')
def eye(n):return [[int(i==j)for j in range(n)]for i in range(n)]
def mul(A,B):return [[sum(A[i][k]*B[k][j]for k in range(len(B)))for j in range(len(B[0]))]for i in range(len(A))]
def power(A,n):
    out=eye(len(A))
    while n:
        if n&1:out=mul(out,A)
        A=mul(A,A);n//=2
    return out
states=[(b,a)for b in(2,0,1)for a in(0,1)];idx={v:k for k,v in enumerate(states)}
def gap(bit):
    A=eye(6)
    for (b,a),j in idx.items():
        A[idx[bit,a]][j]+=1
        if a or bit:A[idx[b,1]][j]+=1
    return A
v0=[[1],[0],[0],[0],[0],[0]];acc=[[2*int(b!=1)for b,a in states]]
def fullF(s,m):
    if s==1:return F(8,25)*9**m+(F(4,15)*m+F(178,225))*4**m+F(8,9)
    return F(128,1225)*729**m+(F(148,315)*m+F(31282,99225))*64**m+F(128,81)
initials={}
for s in(1,3):
    M=mul(power(gap(1),s),power(gap(0),s));ann=eye(6)
    for root,powr in((9**s,1),(4**s,3),(1,2)):
        factor=[[M[i][j]-root*int(i==j)for j in range(6)]for i in range(6)]
        ann=mul(ann,power(factor,powr))
    ck(all(v==0 for row in ann for v in row),'six state annihilator')
    vals=[]
    for m in range(12):
        val=mul(mul(acc,power(M,m)),v0)[0][0];vals.append(val);ck(val==fullF(s,m),'full two cover recurrence')
    initials[s]=vals[:6]
for m in range(1,21):
    delta=144*fullF(3,m)-64*fullF(1,3*m)
    expected=-F(6656,1225)*729**m+(F(576,35)*m-F(19232,3675))*64**m+F(512,3)
    ck(delta==expected and delta<0,'all-mode sender contrast')
if '--fail' in sys.argv:ck(False,'requested failure')
print(json.dumps(dict(status='PASS',literal_template_controls=templates,cycle_floor_controls=floors,literal_binary_tapes=tapes,full_tape_profiles=profiles,full_two_cover_initials=initials,same_order_receiver_twins=[left,right],m0_difference=-16,chan_pak='OPEN'),indent=2,sort_keys=True))
