#!/usr/bin/env python3
"""Standalone exact controls; paired-end generator adapted from S323 with attribution."""
from itertools import product, permutations
from math import isqrt
import json

R13=(24,9,35,4,14,27)
R15=(108,97,2,6,79,92,12)

def need(ok,msg='replay mismatch'):
    if not ok:raise RuntimeError(msg)

def bits(mask):
    while mask:
        b=mask&-mask;yield b.bit_length()-1;mask-=b

def hp_all(r):
    n=len(r);H=[0]*(1<<n);H[0]=1;E=[[0]*n for _ in H]
    for i in range(n):E[1<<i][i]=1
    for S in range(1,1<<n):
        for i in bits(S):
            v=E[S][i];H[S]+=v
            if v:
                for j in bits(r[i]&~S):E[S|(1<<j)][j]+=v
    return H

def hp(r):return hp_all(r)[-1]

def validate(r):
    n=len(r)
    for i,x in enumerate(r):
        need(0<=x<1<<n and not x>>i&1,'invalid row')
        for j in range(i):need((x>>j&1)+(r[j]>>i&1)==1,'invalid tournament')

def strong(r):
    n=len(r);F=(1<<n)-1
    for rows in (r,[F^x^(1<<i)for i,x in enumerate(r)]):
        S=1
        while True:
            T=S
            for i in bits(S):T|=rows[i]
            if T==S:break
            S=T
        if S!=F:return False
    return True

def profile(r):
    H=hp_all(r);F=len(H)-1;Q2=[0]*len(H)
    for S in range(len(H)):
        A=S
        while True:
            Q2[S]+=H[A]*H[S^A]
            if not A:break
            A=(A-1)&S
    return H[F],Q2[F],sum(H[A]*Q2[F^A]for A in range(F+1))

def ordinal(blocks):
    n=sum(map(len,blocks));out=[];off=0
    for r in blocks:
        k=len(r);future=((1<<n)-1)^((1<<(off+k))-1)
        out.extend((x<<off)|future for x in r);off+=k
    return tuple(out)

def inflate(quotient,U,V):
    blocks=[U]+[(0,)]*(len(quotient)-2)+[V]
    offsets=[];m=0
    for B in blocks:offsets.append(m);m+=len(B)
    out=[]
    for i,B in enumerate(blocks):
        cross=sum(((1<<len(blocks[j]))-1)<<offsets[j]for j in bits(quotient[i]))
        out.extend((row<<offsets[i])|cross for row in B)
    validate(out);return tuple(out)

def word_rows(w):
    n=len(w);r=[sum(1<<j for j in range(i+1,n))|((not w[i])<<n)for i in range(n)]
    r.append(sum(a<<i for i,a in enumerate(w)));return tuple(r)

def word_q(w,z):
    n=len(w);out=z**(n+1)
    for i in range(n):
        if w[i]:
            for j in range(i+1,n):
                if not w[j]:
                    d=j-i-1
                    out+=z**(n-d-1)*((z+1)**d+(z-1)**d)
    return out

def candidates(n):
    # S323 paired-end pruning, retaining the final exact dyadic test below.
    if n==3:
        yield (1,0,0);yield (1,1,0);return
    if n<5 or n%4!=1:return
    J=(n-1)//4
    def extend(j,a,b,carry):
        if j==J:
            for v in (0,1):yield tuple(a+[v]+[1-i for i in b[::-1]])
            return
        alpha=(3**(n-2-2*j)+1)//4;beta=(3**(n-3-2*j)+1)//2
        good=0
        for e,f,g,h in product((0,1),repeat=4):
            aa=a+[e,g];bb=b+[f,h]
            C=sum(aa[i]*bb[2*j-i]for i in range(2*j+1))
            D=sum(aa[i]*bb[2*j+1-i]for i in range(2*j+2))
            v=carry+alpha*C+beta*D
            if v%4==0:
                good+=1;yield from extend(j+1,aa,bb,v//4)
        need(good<=(6 if j%2==0 else 4),'carry branching bound')
    alpha=(3**(n-2)+1)//4;beta=(3**(n-3)+1)//2
    for v in (0,1):yield from extend(1,[1,v],[1,1-v],(alpha+beta)//4)

def dyadic_letters(m):
    # (cost,H,word), one least-cost actual word per H.
    best={};n=2
    while 2*3**(n-2)<=1<<m:
        for w in candidates(n):
            q=word_q(w,2)
            if q&(q-1)==0 and q<=1<<m:
                h=word_q(w,1);c=q.bit_length()-1
                if h not in best or c<best[h][0]:best[h]=(c,h,w)
        n+=1
    return sorted(best.values())

def decode(N,A,h):
    if h<=0 or A<=h or N<=A:return None
    m=1
    while A*(1<<m)<=N:m+=1
    if not A*(1<<(m-1))<N<A*(1<<m):return None
    t=A*(1<<m)-N
    if t%h:return None
    D=t//h
    if D<=0 or D%2==0 or D>1<<(m-1):return None
    letters=dyadic_letters(m);visited=0
    def dfs(start,cost,H,parts):
        nonlocal visited;visited+=1
        if H==D:return (m-cost,parts)
        for i in range(start,len(letters)):
            c,d,w=letters[i]
            if c+cost>m:break
            if H*d>D:continue
            ans=dfs(i,cost+c,H*d,parts+(w,))
            if ans is not None:return ans
        return None
    ans=dfs(0,0,1,())
    return None if ans is None else dict(m=m,D=D,padding=ans[0],words=ans[1],visited=visited)

def make_reservoir(ans):
    return ordinal([word_rows(w)for w in ans['words']]+[(0,)]*ans['padding'])

def source_letters(Y):
    ans=[(3,(0,),(1,2,3))];n=2
    while 3*4**(n-2)<=Y:
        for mid in product((0,1),repeat=n-2):
            w=(1,)+mid+(0,);p=tuple(word_q(w,z)for z in (1,2,3))
            if p[2]<=Y:ans.append((p[2],word_rows(w),p))
        n+=1
    return ans

def sourcefree(N):
    cap=N//2
    # Tagged alphabets; no score, H, or graph-order shortcut is used.
    letters=[(weight,0,rows,p)for weight,rows,p in source_letters(cap)]
    letters += [(2,1,(0,),(1,2))]
    for cost,h,w in dyadic_letters(cap.bit_length()-1):
        q=1<<cost
        if q<=cap:letters.append((q,1,word_rows(w),(h,q)))
    letters.sort(key=lambda z:(z[0],z[1],z[2]))
    visited=0
    def dfs(start,weight,up,vp,Ub,Vb):
        nonlocal visited;visited+=1
        if Ub and Vb:
            h,b,c=up;d,q=vp
            for gamma in (2,6):
                if (-3*h+6*b+gamma*c)*q-h*d==N:
                    return gamma,ordinal(Ub),ordinal(Vb)
        for i in range(start,len(letters)):
            w,role,r,p=letters[i]
            if weight*w>cap:break
            if role==0:
                ans=dfs(i,weight*w,tuple(a*b for a,b in zip(up,p)),vp,Ub+(r,),Vb)
            else:ans=dfs(i,weight*w,up,tuple(a*b for a,b in zip(vp,p)),Ub,Vb+(r,))
            if ans:return ans
        return None
    ans=dfs(0,1,(1,1,1),(1,1),(),())
    return ans,visited

def isprime(n):
    return n>=2 and all(n%p for p in range(2,isqrt(n)+1))

def main():
    # Independent bounded validation of the inherited carry compiler.
    cutchecks=0
    for n in range(2,15):
        brute=set()
        for mid in product((0,1),repeat=n-2):
            w=(1,)+mid+(0,);q=word_q(w,2)
            if not q&(q-1):brute.add(w)
            cutchecks+=1
        pruned={w for w in candidates(n)if not word_q(w,2)&(word_q(w,2)-1)}
        need(brute==pruned,'complete carry pruning')
    # Actual graphs check Q(1),Q(2),Q(3) rather than the compiler formula alone.
    graphprofiles=0
    for n in range(2,6):
        for mid in product((0,1),repeat=n-2):
            w=(1,)+mid+(0,);r=word_rows(w)
            need(strong(r),'active core strong')
            need(profile(r)==tuple(word_q(w,z)for z in (1,2,3)),'literal slot profile')
            graphprofiles+=1
    examples=[]
    for N,A,h,gamma in ((3821,15,1,2),(859,27,1,6),(27629,27,1,6)):
        z=decode(N,A,h);need(z is not None,'positive factor-free inverse')
        V=make_reservoir(z);d,q,_=profile(V)
        need(A*q-h*d==N and q==1<<z['m'] and d==z['D'],'reservoir exact target')
        T=inflate(R13 if gamma==2 else R15,(0,),V)
        need(strong(T)and hp(T)==N,'actual strong quotient target')
        examples.append(dict(N=N,A=A,h=h,gamma=gamma,m=z['m'],D=z['D'],padding=z['padding'],words=z['words'],order=len(T),visited=z['visited']))
    need(isprime(859)and isprime(3947)and 27629==7*3947,'two atom lanes')
    # Budget-complete inversion tests against direct two-letter exponent products.
    # Through m=18 the complete generator supplies exactly the verified 5 and 19 letters.
    matches=0;rejects=0
    for m in range(1,19):
        letters=dyadic_letters(m)
        need({(c,h)for c,h,w in letters}<={(5,5),(8,19)},'bounded alphabet changed')
        possible={5**a*19**b for a in range(m//5+1)for b in range(m//8+1)if 5*a+8*b<=m}
        tests=set(possible)|set(range(1,min(101,1<<m),2))
        for D in tests:
            N=15*(1<<m)-D
            if not 15*(1<<(m-1))<N<15*(1<<m):continue
            ans=decode(N,15,1)
            need((ans is not None)==(D in possible),'multiset inverse completeness on bounded budget')
            if ans is None:rejects+=1
            else:matches+=1
    free=[]
    for N in (7,21,29,53,59,119,859,3821,27629):
        ans,vis=sourcefree(N)
        if N in (7,21):need(ans is None,'small exclusions below chart minimum')
        else:
            need(ans is not None,'source-free positive target')
            gamma,U,V=ans;T=inflate(R13 if gamma==2 else R15,U,V)
            need(strong(T)and hp(T)==N,'source-free actual graph')
        free.append(dict(N=N,found=ans is not None,visited=vis))
    # Removing either kind of component is checked in actual inflated graphs.
    tri=word_rows((1,0));five=word_rows((1,0,0));one=(0,)
    peel=0;commute=0
    for gamma,Q in ((2,R13),(6,R15)):
        for Ub,Vb in (((one,tri),(one,one)),((five,one),(tri,one)),((tri,one),(five,))):
            U=ordinal(Ub);V=ordinal(Vb);T=inflate(Q,U,V);N=hp(T)
            p=profile(U);v=profile(V)
            need(N==(-3*p[0]+6*p[1]+gamma*p[2])*v[1]-p[0]*v[0],'selective chart formula')
            for j in range(len(Ub)):
                if len(Ub)>1:
                    child=inflate(Q,ordinal(Ub[:j]+Ub[j+1:]),V);M=hp(child)
                    need(strong(child)and 0<2*M<N,'physical source halving');peel+=1
            for j in range(len(Vb)):
                if len(Vb)>1:
                    child=inflate(Q,U,ordinal(Vb[:j]+Vb[j+1:]));M=hp(child)
                    need(strong(child)and 0<2*M<N,'physical reservoir halving');peel+=1
            for uu in permutations(Ub):
                for vv in permutations(Vb):
                    need(profile(ordinal(uu))==p and profile(ordinal(vv))==v,'commuting slot profiles')
                    commute+=1
    print(json.dumps(dict(status='PASS: finite controls; all-order results are proved in COMMUTATIVE_INVERSE.md',carry_words_checked=cutchecks,one_apex_graph_profiles=graphprofiles,factor_free_examples=examples,bounded_inverse_accepts=matches,bounded_inverse_rejects=rejects,source_free_controls=free,physical_halving_checks=peel,commutative_profile_controls=commute),indent=2,sort_keys=True))

if __name__=='__main__':main()
