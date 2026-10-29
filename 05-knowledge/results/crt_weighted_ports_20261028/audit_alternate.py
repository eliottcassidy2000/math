#!/usr/bin/env python3
"""Same-agent alternate methods: literal words, portraits, powers, LCA kernels."""
from fractions import Fraction as F
from itertools import product
import json
from crt_probe import make
from linear_repair import repair

CHECKS=0
def ck(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(msg)
def digits(x,k):return tuple(x//3**d%3 for d in range(k-1,-1,-1))
def number(w):
    x=0
    for d in w:x=3*x+d
    return x
def compose(a,b):return tuple(a[x] for x in b)
def power(p,e):
    out=tuple(range(len(p)))
    while e:
        if e%2:out=compose(out,p)
        p=compose(p,p);e//=2
    return out
def portrait_groups(k,centralizer):
    prefixes=[w for d in range(k) for w in product(range(3),repeat=d)]
    reps=[w for w in prefixes if any(w) and next(d for d in w if d)==1] if centralizer else prefixes
    for values in product(range(3),repeat=len(reps)):
        shifts={w:0 for w in prefixes}
        for w,a in zip(reps,values):
            shifts[w]=a
            if centralizer:shifts[tuple((-d)%3 for d in w)]=(-a)%3
        yield tuple(number([(d+shifts[w[:i]])%3 for i,d in enumerate(w)])
                    for w in product(range(3),repeat=k))
def arc(a,b):
    for u,v in zip(a,b):
        if u!=v:return (v-u)%3==1
    return False
def main():
    crt=[]
    for n in [9,15]:
        qs,strides,R,W=make(n)
        roots=[2 if q==4 else next(a for a in range(q) if (a*a-a-1)%q==0 and pow(a,n,q)==1) for q in qs]
        count=0
        def read(mask,reverse=False):
            coords=[]
            for q,a in zip(qs,roots):
                value=0
                for i in range(n):
                    if mask>>i&1:
                        j=(-i)%n if reverse else i
                        if q==4:value^=[1,2,3][j%3]
                        else:value=(value+pow(a,j,q))%q
                coords.append(value)
            return sum(s*x for s,x in zip(strides,coords))
        for mask in range(1<<n):
            if mask&(mask<<1) or ((mask&1) and (mask>>(n-1)&1)):continue
            x=read(mask);y=read(mask,True);count+=1
            ck(W[x]==mask and R[x]==y,'literal bit-mask reader matches recursive CRT')
        ck(count==len(R),'literal word universe')
        crt.append(dict(n=n,words=count))
    loops=[]
    for k in [2,3]:
        n=3**k;words=[digits(x,k) for x in range(n)]
        alpha=tuple(number([(-d)%3 for d in w]) for w in words)
        ident=tuple(range(n));groups=list(portrait_groups(k,True))
        ck(len(set(groups))==3**((n-2*k-1)//4),'independent complete centralizer portraits')
        for c in groups:
            actual=ident
            for t,u in repair(k,c):
                for mark in [t,u]:
                    ck(compose(mark,mark)==ident,'full returned mark involution')
                    for x in range(n):
                        for y in range(x):ck(arc(words[mark[x]],words[mark[y]])!=arc(words[x],words[y]),'full returned mark anti')
                # Aut(W_k) has exponent dividing 3^k; use powers, not cyclewise roots.
                e=(3**k+1)//2
                gs=[power(compose(a,b),e) for a,b in [(t,alpha),(u,t),(u,alpha)]]
                for g,(a,b) in zip(gs,[(t,alpha),(u,t),(u,alpha)]):ck(compose(g,g)==compose(a,b),'independent root by group exponent')
                last_inverse=power(gs[2],3**k-1)
                actual=compose(actual,compose(last_inverse,compose(gs[1],gs[0])))
            ck(actual==c,'portrait target equals synthesized full permutation')
        loops.append(dict(k=k,full_centralizer=len(groups)))
    # A nonproduct source law; literal finite group average versus LCA kernel.
    G=list(portrait_groups(2,False));n=9
    raw=[2,3,5,7,11,13,17,19,23];mu=[F(x,sum(raw)) for x in raw]
    prefix=[sum(mu[3*i:3*i+3]) for i in range(3)]
    expected=1+2*(1+3*sum(x*x for x in prefix))
    actual=sum((9*sum(mu[i] for i in range(n) if g[i]==i))**2 for g in G)/len(G)
    ck(actual==expected,'literal nonproduct weighted energy')
    kernel=0
    for x in range(n):
        for y in range(n):
            common=0
            for a,b in zip(digits(x,2),digits(y,2)):
                if a!=b:break
                common+=1
            value=3**min(2,common+1)
            ck(F(sum(g[x]==x and g[y]==y for g in G),len(G))*81==value,'literal pair law equals LCA kernel')
            kernel+=mu[x]*mu[y]*value
    ck(kernel==expected,'LCA source quadratic form')
    # Sibling signed null mode at depth1 becomes energy12 at depth2.
    for k in [1,2]:
        G=list(portrait_groups(k,False));signed=[F(1,3**(k-1))]*3**(k-1)+[F(-1,3**(k-1))]*3**(k-1)+[F(0)]*3**(k-1)
        second=sum((3**k*sum(signed[i] for i in range(3**k) if g[i]==i))**2 for g in G)/len(G)
        ck(second==[0,12][k-1],'one-generation activation of signed kernel')
    print(json.dumps(dict(crt=crt,loops=loops,nonproduct_second_moment=str(expected),active_checks=CHECKS),indent=2,sort_keys=True))
if __name__=='__main__':main()
