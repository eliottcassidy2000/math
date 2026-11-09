#!/usr/bin/env python3
"""Self-contained exact checks: no third-party packages or repository imports."""
import itertools, json, random
from fractions import Fraction
from pathlib import Path

def tournament(n, bits):
    a=[0]*n
    for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
        if bits>>k&1: a[i]|=1<<j
        else: a[j]|=1<<i
    return a

def converse(a):
    return [sum(1<<j for j in range(len(a)) if a[j]>>i&1) for i in range(len(a))]

def counts(a):
    n=len(a); N=1<<n; dp=[[0]*n for _ in range(N)]; h=[0]*N; h[0]=1
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,N):
        for v in range(n):
            ways=dp[mask][v]
            if not ways: continue
            avail=a[v]&~mask
            while avail:
                bit=avail&-avail; avail-=bit
                dp[mask|bit][bit.bit_length()-1]+=ways
        h[mask]=sum(dp[mask])
    return h

def state(a,u):
    h=counts(a); full=(1<<len(a))-1; eo=[0,0]
    def dfs(v,mask):
        eo[mask.bit_count()%2]+=h[full^mask]
        avail=a[v]&~mask
        while avail:
            bit=avail&-avail; avail-=bit; dfs(bit.bit_length()-1,mask|bit)
    dfs(u,1<<u)
    return (h[-1],*eo)

def pad(a,m):
    n=len(a)
    return [row|sum(1<<j for j in range(n,n+m)) for row in a]+[sum(1<<j for j in range(i+1,n+m)) for i in range(n,n+m)]

def seam(a,u,b,w):
    n=len(a); k=len(b)
    c=[row|(((1<<k)-1)<<n) for row in a]+[row<<n for row in b]
    c[u]^=1<<(n+w); c[n+w]|=1<<u
    return c

def require(ok,msg):
    if not ok: raise RuntimeError(msg)

def apex_h(w):
    if w==0: return 1
    old=w//2
    return apex_h(old)+(old+old%2 if w%2==0 else 0)

def block_data(bits):
    alpha=beta=0;c=1;e=0;p=bits[-1]
    for bit in bits:
        if bit==0:alpha+=c;beta+=e+p
        c,e,p=2*c,2*e+bit,bit
    return c,e,alpha,beta

def apex_tournament(bits):
    n=len(bits);a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]+[0]
    for i,bit in enumerate(bits):
        if bit:a[n]|=1<<i
        else:a[i]|=1<<n
    return a

def check_periodic():
    checked=0;direct=0;examples=[]
    for k in range(1,11):
        for bits in itertools.product((0,1),repeat=k):
            q,d,alpha,beta=block_data(bits)
            D=Fraction(beta)-Fraction(alpha*d,q-1)
            mixed=0<sum(bits)<k
            require((D>0)==mixed,'mixed drift sign')
            require(mixed or D==0,'constant word drift')
            w0=2+bits[-1];w=w0;H0=apex_h(w0)
            A=Fraction(alpha,q-1)*(w0+Fraction(d,q-1))
            for m in range(5):
                require(apex_h(w)==A*q**m+D*m+H0-A,'periodic exact formula')
                w=q*w+d
            if k<=4 and mixed:examples.append(dict(word=''.join(map(str,bits)),alpha=alpha,beta=beta,drift=str(D)))
            checked+=1
    for k in range(1,8):
        for bits in itertools.product((0,1),repeat=k):
            w=int(''.join(map(str,bits)),2)
            require(counts(apex_tournament(bits))[-1]==apex_h(w),'one-apex word orientation')
            direct+=1
    # Omitted initial-parity sidecar produces the wrong first block.
    q,d,alpha,beta=block_data((0,1))
    require(apex_h(d)!=apex_h(0)+beta,'parity hostile')
    require(Fraction(2,3)==Fraction(beta)-Fraction(alpha*d,q-1),'smallest mixed example')
    # Negative zero-tail offset exists, but its leading coefficient is not a power of two.
    for m in range(12):require(apex_h(6*2**m)==6*2**m-1,'negative offset')
    for r in range(1,20):require(apex_h(2**r)==2**r+1,'pure clock offset')
    return dict(periodic_words=checked,max_period_length=10,direct_apex_tournaments=direct,max_direct_word_length=7,small_examples=examples)

def main():
    roots=[]
    for n in range(1,4):
        for bits in range(1<<(n*(n-1)//2)):
            a=tournament(n,bits)
            roots.extend((a,u) for u in range(n))
    cases=0; padding=0
    for a,u in roots:
        H,E,O=state(a,u)
        for m in range(1,4):
            require(state(pad(a,m),u)==(H,2**(m-1)*(E+O),2**(m-1)*(E+O)), 'padding')
            padding+=1
        for b,w in roots:
            K,Qe,Qo=state(converse(b),w)
            for m in (1,2):
                actual=counts(seam(pad(a,m),u,b,w))[-1]
                expected=H*K+2**m*(E+O)*(Qe+Qo)
                require(actual==expected,'bridge')
                cases+=1
    rng=random.Random(236)
    for _ in range(50):
        n=rng.randrange(3,6); k=rng.randrange(2,5); m=rng.randrange(1,4)
        a=tournament(n,rng.randrange(1<<(n*(n-1)//2))); b=tournament(k,rng.randrange(1<<(k*(k-1)//2)))
        u=rng.randrange(n);w=rng.randrange(k)
        H,E,O=state(a,u);K,Qe,Qo=state(converse(b),w)
        require(counts(seam(pad(a,m),u,b,w))[-1]==H*K+2**m*(E+O)*(Qe+Qo),'random')
        cases+=1
    examples=[]
    for bits,u,C,r in ((9,1,5,5),(8,2,9,4)):
        a=converse(tournament(5,bits));H,E,O=state(a,u)
        require((H,E+O)==(C,2**r),'ray seed')
        for m in range(1,5):
            value=counts(seam(pad(a,m),u,[0],0))[-1]
            require(value==2**(m+r)+C,'ray')
        examples.append(dict(adjacency=a,port=u,h=H,response=E+O,formula=f'2^(m+{r})+{C}',m_min=1))
    # Hostile: padding cannot suppress the unequal initial parity coordinates at m=0.
    require(state([0],0)==(1,0,1),'hostile seed')
    require(counts(seam([0],0,[0],0))[-1]==1,'m0 actual')
    require(1!=1+2**0,'m0 extrapolation must fail')
    # Hostile: equal h does not determine the dyadic coefficient.
    require(state(tournament(2,1),0)[0]==state(tournament(2,1),1)[0]==1,'same h')
    require(sum(state(tournament(2,1),0)[1:])!=sum(state(tournament(2,1),1)[1:]),'response distinction')
    periodic=check_periodic()
    out=dict(periodic=periodic,status='FINITE-EXACT verification of separately proved all-order identities',exhaustive_rooted_tournaments_orders='1..3, labelled',padding_checks=padding,bridge_checks=cases,positive_offset_rays=examples,hostile_controls=5)
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
