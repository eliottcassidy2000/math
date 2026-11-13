"""Literal controls for a same-target one-arc repair and dyadic compensation."""
from itertools import permutations
from math import gcd
import numpy as np
from inflation_probe import zprofile,end_table,align
from fibre_probe import HERE,save,need

K9=(4,33,26,35,43,5)
K11=(4,37,24,35,43,5)
K3=(2,32,11,19,39,13)

def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def near(n):
    T=list(transitive(n));T[0]^=1<<(n-1);T[-1]|=1;return tuple(T)

def join(A,B):
    n=len(A);high=((1<<len(B))-1)<<n
    return tuple(t|high for t in A)+tuple(t<<n for t in B)

def birth(T,r):
    n=len(T);U=list(T)+[((1<<n)-1)^(1<<r)];U[r]|=1<<n;return tuple(U)

def profile(T,r):return tuple(map(int,zprofile(np.array(T,dtype=np.int64),r)))

def hp2(T):
    D=end_table(np.array(T,dtype=np.int64));H=list(map(int,D.sum(axis=1)));H[0]=1;full=len(H)-1
    return H[-1],sum(H[S]*H[full^S] for S in range(full+1))

def strong(T):
    n=len(T);R=[T[i]|1<<i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i]>>k&1:R[i]|=R[k]
    return all(v==(1<<n)-1 for v in R)

def literal_profile(T,r):
    def path(p):return all(T[a]>>b&1 for a,b in zip(p,p[1:]))
    h=s=F=0
    for p in permutations(range(len(T))):
        ok=path(p);h+=ok;s+=ok and p[0]==r;cut=p.index(r)+1;F+=path(p[:cut]) and path(p[cut:])
    q=sum(path(p) for p in permutations(i for i in range(len(T)) if i!=r))
    return h,s,q,F

def chain(T,r,k):
    need(k>=2,'chain guard');z=profile(T,r)
    for _ in range(k):T=birth(T,r)
    w=len(T)-1;T=birth(T,r);y=len(T)-1;T=birth(T,y)
    h,s,q,F=z;v=1<<(k-1);d=F-q;c=h-s+4*q-F
    S=v*d;Q=3*S+c;expected=(Q+3*S,S,Q,4*Q+2*S)
    observed=profile(T,w)
    need(observed==expected and strong(T),'literal chart and strong entry')
    return T,w,dict(k=k,rows=T,root=w,ports=[y,w],profile=observed,content=gcd(Q,S))

def run():
    base=[]
    for T,expect in ((K9,(23,2,9,60)),(K11,(23,2,11,60)),(K3,(23,2,3,30))):
        need(profile(T,5)==expect==literal_profile(T,5) and strong(T),'independent base control')
        base.append(dict(rows=T,root=5,profile=expect))
    need([i for i in range(6) if K9[i]!=K11[i]]==[1,2],'one off-root flip')
    need(K9[1]^(1<<2)==K11[1] and K9[2]^(1<<1)==K11[2],'literal flip')
    B,flips,permutation=align(K9,K3,5)
    controls=[]
    factors=[transitive(a) for a in range(1,6)]+[(2,4,1),near(4),near(6)]
    for L in factors:
        h,f=hp2(L);states=[]
        for q,K,extra in ((9,K9,False),(11,K11,False),(3,B,True)):
            left=join(L,(0,)) if extra else L;T=join(left,K);r=len(left)+5;z=profile(T,r)
            need(z==(23*h,0,q*h,60*f),'prefix law and compensated amplitude')
            U=birth(T,r);u=profile(U,r);M=23*h+60*f
            need(u==(M,q*h,q*h,120*f-q*h) and strong(U),'exact next target')
            expected_g=gcd((23+3*q)*h,60*f-q*h)
            # One ordinary R has already occurred. Every later L_k has this odd content.
            V,w,entered=chain(U,r,2)
            need(entered['content']==expected_g,'exact chain entry gcd')
            need(T[0]>>r&1 and U[0]>>r&1,'retained prefix incoming port')
            states.append(dict(q=q,extra_singleton=extra,source_rows=T,root=r,source_profile=z,next_rows=U,next_profile=u,entry=entered,formula_content=expected_g))
        if gcd(h,30*f)==1:
            need(states[0]['formula_content']==states[2]['formula_content']==1,'repaired universal guard')
            need(states[1]['formula_content']==(7 if (f-h)%7==0 else 1),'exact factor7 phase')
        controls.append(dict(prefix_rows=L,h=h,f=f,guard=gcd(h,30*f)==1,states=states))
    # A full second physical return, not just an arithmetic orbit.
    T=join(transitive(1),K9);r=len(T)-1;T=birth(T,r)
    T,r,first=chain(T,r,2);T,r,second=chain(T,r,2)
    need(first['content']==second['content']==1,'literal renewable entry')
    clocks=[]
    for a in range(1,121):
        x=1<<a;M=23+60*x;g9=gcd(50,60*x-9);g11=gcd(56,60*x-11);g3=gcd(32,60*x-3)
        need(g9==g3==1 and g11==(7 if a%3==0 else 1),'all-phase arithmetic control')
        clocks.append([a,M,g9,g11,g3])
    # R6 has H17 and P2=224, giving an unbounded nontrivial prefix grammar.
    need(hp2(near(6))==(17,224),'nontrivial prefix generator')
    grid=[]
    for b in range(9):
        for a in range(9):
            if a==b==0:continue
            h=17**b;f=224**b*2**a
            need(gcd(h,30*f)==1,'prefix grammar reentry')
            grid.append([a,b,23*h+60*f])
    save('repair23',dict(scope='Literal3 base cores,8 prefix controls and all three next-target variants, physical chain entries through18 vertices, a full second chain return; universal formulas and guards are proved in REPAIR23.md',base=base,compensated_core_rows=B,compensated_core_flips=flips,compensated_core_permutation=permutation,prefix_controls=controls,literal_renewals=[first,second],arithmetic_clock_controls=clocks,prefix_grid_controls=grid))
    print('PASS',len(controls),'prefix controls;3 preparations each; two literal chain entries')
    print('COMPENSATION',B,'flips',flips)
    print('FIRST_RENEWALS',first['profile'],second['profile'])
    print('FIRST_TARGETS',[v[1] for v in clocks[:7]])

if __name__=='__main__':run()
