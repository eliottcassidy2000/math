from pathlib import Path
from math import isqrt
import json,random
STATES=[tuple(map(int,l.split())) for l in Path(__file__).parents[1].joinpath('CORE_STATES.txt').read_text().splitlines()]
def f(n):
    ans=0
    while n:
        if n&1:ans+=(n>>1)+((n>>1)&1)
        n>>=1
    return ans
def rev(n,b):
    ans=0
    for _ in range(b):ans=(ans<<1)|(n&1);n>>=1
    return ans
def val(c,p,r,b):
    h,q,A,G,T,S,E=STATES[c]
    return (q*p+A)*((2<<b)-rev(r,b))+h*(2*r+1-f(r))-p*G-h*f(p)-T-(p&1)*S-E*(r>>(b-1) if b else 0)
def ceildiv(x,y):return -((-x)//y)
def solve(N):
    out=set()
    for c,(h,q,A,G,T,S,E) in enumerate(STATES):
      b=0
      while h*(1<<b)<=N-h:
        W=1<<b
        if W*W<N:
          for r in range(W):
            K=q*(2*W-rev(r,b))-G
            C=A*(2*W-rev(r,b))+h*(2*r+1-f(r))-T-E*(r>>(b-1) if b else 0)
            depth=0
            while q*q*W*W*4**depth<2*h*N:depth+=1
            tree=[(K,S,C,0)]
            for t in range(depth):
                tree=[z for K,s,C,l in tree for z in ((2*K,0,C,l),(2*K-h,h,C+K-s,l+(1<<t)))]
            for K,s,C,l in tree:
                x=N-C
                lo=max(1,ceildiv(x,K));hi=(x-h+s)//(K-h)
                if lo<hi:raise Exception(('prefix-width',N,c,b,r,lo,hi))
                candidates=[0] if x==0 else []
                if lo==hi:candidates.append(lo)
                for u in candidates:
                    p=l+(u<<depth)
                    if (p or A) and val(c,p,r,b)==N:out.add((c,p,r,b))
        else:
          for p in range(max(0,((N-h)//W-A)//q)+1):
            alpha=q*p+A
            if not alpha or W*alpha>N-h:continue
            tau=p*G+h*f(p)+T+(p&1)*S
            X=2*W*alpha+h-tau-N
            t=2
            while alpha*4**t<24*h*W:t+=1
            t=min(t,b)
            tree=[(alpha,E,0,0)]
            for i in range(t):
              tree=[z for K,s,C,l in tree for z in ((2*K,0,C,l),(2*K+h,h,C+K+s-2*h*(1<<(b-i-1)),l+(1<<i)))]
            m=b-t
            for K,s,C,l in tree:
              x=X-C
              lo=max(0,ceildiv(x-h,K+h));hi=min((1<<m)-1,(x+2*h*((1<<m)-1))//K)
              if lo<hi:raise Exception(('suffix-width',N,c,p,b,lo,hi))
              if lo==hi:
                z=l+(lo<<t);r=rev(z,b)
                if val(c,p,r,b)==N:out.add((c,p,r,b))
        b+=1
    return out

def brute_table(M):
    out={}
    # Independent exact all-length universe from alpha<=M-h.
    for c,(h,q,A,G,T,S,E) in enumerate(STATES):
      for p in range((M-h-A)//q+1):
        a0=q*p+A
        if not a0:continue
        b=0
        while (1<<b)*a0<=M-h:
          for r in range(min((1<<b)-1,(M-h-(1<<b)*a0)//h)+1):
            N=val(c,p,r,b)
            if N<=M:out.setdefault(N,set()).add((c,p,r,b))
          b+=1
    return out
M=1500
tab=brute_table(M)
for N in range(1,M+1):
    got=solve(N);expected=tab.get(N,set())
    if got!=expected:raise Exception(('mismatch',N,sorted(got-expected)[:5],sorted(expected-got)[:5]))
large=10000
large_tab=brute_table(large)
rng=random.Random(217)
extra=[265,287,1237,4095,4097,9999]+rng.sample(range(M+1,large+1),100)
for N in extra:
 if solve(N)!=large_tab.get(N,set()):raise Exception(('large-mismatch',N))
print(json.dumps({'status':'PASS','all_targets':[1,M],'compared':'entire recipe sets; exact all-length alpha-pruned brute','core_states':80,'brute_recipe_total':sum(map(len,tab.values())),'extra_targets':len(extra),'extra_bound':large,'includes_even_targets':True,'omissions_47_63':not tab.get(47) and not tab.get(63)},indent=2))
