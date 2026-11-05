#!/usr/bin/env python3
"""Independent exact checks. Run normally and with python3 -O; compare JSON."""
from collections import Counter, deque
from hashlib import sha256
from itertools import product
from math import isqrt
from pathlib import Path
import json

T=(1,0,2,1)
B=(3,1,3,1)

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def move(s,c):
    D,E,U,d=s
    if c=='O': return D,0,2*U,d
    if c=='Z': return D-E+U,d,2*U-d,d
    if c=='J': return D-E+U,D-E,D+U,D
    raise ValueError(c)

def run(s,word):
    for c in word: s=move(s,c)
    return s

def value(s,closure):
    D,E,U,d=s
    return D-E+U if closure=='F' else D+2*U

def K(v,w):
    D,E,U,d=v; A,B,C,a=w
    return D*a+d*A+E*B+U*C-E*a-d*B-U*a-d*C+d*a

def graph(word):
    adj=[2,0]; root=0
    for c in word:
        n=len(adj); adj.append(0)
        for i in range(n):
            if c in 'ZJ' and i==root: adj[n]|=1<<i
            else: adj[i]|=1<<n
        if c=='J':root=n
    return adj,root

def counts(adj,root=0):
    n=len(adj); size=1<<n
    # First implementation: paths ending at each vertex, for every induced subset.
    dp=[[0]*n for _ in range(size)]
    start=[[0]*n for _ in range(size)]
    for j in range(n):dp[1<<j][j]=1
    start[1<<root][root]=1
    for mask in range(1,size):
        for end in range(n):
            if not (mask>>end)&1:continue
            rem=adj[end]&~mask
            while rem:
                bit=rem&-rem; rem-=bit; j=bit.bit_length()-1
                dp[mask|bit][j]+=dp[mask][end]
                start[mask|bit][j]+=start[mask][end]
    total=[sum(row) for row in dp];total[0]=1
    full=size-1
    H=total[full]; E=dp[full][root]; S=sum(start[full]); d=total[full^(1<<root)]
    R=0
    for a in range(size):
        if (a>>root)&1:continue
        R+=total[a]*sum(start[full^a])
    C=sum(start[full][j] for j in range(n) if (adj[j]>>root)&1)
    return (H,E,R,d), S, C

def critical_fibers(prefix,k,r,closure):
    s=run(T,prefix); mod=1<<r; hist=Counter()
    # Indexed words are intentionally not deduplicated by profile.
    states=[s]
    for _ in range(k):
        states=[move(v,c) for v in states for c in 'OZ']
    for v in states:hist[value(v,closure)%mod]+=1
    return hist

def polynomial(s,word,closure):
    D,E,U,d=s; k=len(word)
    if k==0:return value(s,closure)
    bits=[c=='Z' for c in word]
    aa=sum((1<<i)*b for i,b in enumerate(bits))
    rr=sum((1<<(k-1-i))*b for i,b in enumerate(bits))
    cc=sum(((1<<(j-i-1))+(j==i+1))*bits[i]*bits[j]
           for i in range(k) for j in range(i+1,k))
    out=(D+U*aa-E*bits[0]-d*cc,d*bits[-1],(1<<k)*U-d*rr,d)
    need(out==run(s,word),'quadratic profile formula')
    return value(out,closure)

def bank(cap):
    seen={T:''}; todo=deque([T])
    while todo:
        s=todo.popleft()
        for c in 'OZJ':
            v=move(s,c)
            if v[2]<=cap and v not in seen:
                seen[v]=seen[s]+c;todo.append(v)
    return seen

def main():
    checks=0; prefix_words=['','Z','ZJ','OZJ','OJZJ']
    parity_receipts=[]
    for r in range(1,5):
        critical=(1<<r)-2
        for prefix in prefix_words:
            for closure in 'FG':
                for k in (critical,critical+1,critical+2):
                    hist=critical_fibers(prefix,k,r,closure)
                    need(set(hist)==set(range(1,1<<r,2)),'missing/extra dyadic residue')
                    need(all((v%2)==(k==critical) for v in hist.values()),'sharp parity')
                    checks+=1
                    if prefix=='':
                        parity_receipts.append({'r':r,'k':k,'closure':closure,'fibers':dict(sorted(hist.items()))})
    # The general degree bound must also hold at noncritical lengths.
    for prefix in prefix_words[:3]:
        for k in range(0,13):
            s=run(T,prefix)
            for word_tuple in product('OZ',repeat=k):
                word=''.join(word_tuple)
                for closure in 'FG':polynomial(s,word,closure)
                checks+=1
    # Literal graph controls use no profile recursion for their counts.
    literal_graphs=0
    for prefix in ('','ZJ','OZ'):
        for k in range(5):
            for wt in product('OZ',repeat=k):
                word=prefix+''.join(wt)
                adj,root=graph(word)
                prof,starts,_=counts(adj,root)
                need(prof==run(T,word) and starts==1,'literal source response')
                for closure,tail in [('F','J'),('G','OJ')]:
                    a,r=graph(word+tail);p,_,cycles=counts(a,r)
                    need(p[0]==value(prof,closure),'literal closure count')
                    need(cycles==1,'literal closure cycles')
                    literal_graphs+=1
    # Critical graph support: parity of square-path perfect matchings.
    fib=[0,1]
    for i in range(2,513):fib.append(fib[-1]+fib[-2])
    for r in range(1,11):
        t=(1<<(r-1))-1
        need(fib[t+1]%2==1,'critical matching parity')
        # Match counts have independent two-branch recurrence a_t=F_(t+1).
        aa,bb=1,1
        if t==0:a=aa
        elif t==1:a=bb
        else:
            for _ in range(2,t+1):aa,bb=bb,aa+bb
            a=bb
        need(a==fib[t+1],'square-path matching recurrence')
        checks+=1
    sharp=[]
    for a in range(1,6):
        h=(1<<a)+1; v=run(T,'Z'*a);w=run(B,'Z'*(a-1))
        need(v==w==(h,1,h,1),'primitive balanced sources')
        need(K(v,w)==h*h,'sharp boundary pair')
        adj,root=graph('Z'*a+'O'*a+'J');prof,_,cyc=counts(adj,root)
        need(prof[0]==h*h and cyc==1,'sharp boundary literal graph')
        literal_graphs+=1
        sharp.append({'a':a,'height':h,'target':7,'modulus':h*h-7,'H':h*h})
    small=bank(20); certs=0
    for v in small:
        for w in small:
            h=K(v,w);p=v[2]*w[2]
            need(p//2+v[3]*w[3]<=h<=p,'physical K bounds')
            need(h%2==1,'physical K oddness')
            for m in range(1,41):
                q=m if m%2==0 else 2*m
                n=h%q
                if n>=3 and q>n and p<n+q:
                    need(h==n,'height certificate exactness');certs+=1
    cutoff=isqrt(4*(7-1)); tiny=bank(cutoff)
    support=sorted({K(v,w) for v in tiny for w in tiny if K(v,w)<=7})
    need(support==[3,5],'C1 finite hostile at 7')
    root=Path(__file__).resolve().parent
    print(json.dumps({'status':'PASS','checks':checks,'literal_graphs':literal_graphs,
        'height_certificates':certs,'small_source_bank':len(small),
        'dyadic_fibers':parity_receipts,'sharpness':sharp,
        'C1_hostile':{'target':7,'U_cutoff':cutoff,'bank_size':len(tiny),'support':support},
        'proof_sha256':sha256((root/'PROOF.md').read_bytes()).hexdigest(),
        'verifier_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))

if __name__=='__main__':main()
