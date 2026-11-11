#!/usr/bin/env python3
"""Exact new constructor tests; no inherited result census. Standard library only."""
import itertools, json
from collections import Counter
from functools import lru_cache
from pathlib import Path

def apex(w):
    n=len(w)+1; a=[0]*n
    for i in range(1,n):
        for j in range(i+1,n): a[i]|=1<<j
        if w[i-1]=='1': a[0]|=1<<i
        else: a[i]|=1
    return tuple(a)

def join(A,B,alpha=None,beta=None):
    n=len(A); m=len(B)
    out=[x|(((1<<m)-1)<<n) for x in A]+[x<<n for x in B]
    if alpha is not None:
        out[alpha]&=~(1<<(n+beta)); out[n+beta]|=1<<alpha
    return tuple(out)

def paths(A,start=None,end=None):
    out=[]
    def walk(p,used):
        if end is None or p[-1]==end: out.append(p)
        if end is not None and p[-1]==end: return
        for v in range(len(A)):
            if not used>>v&1 and A[p[-1]]>>v&1: walk(p+(v,),used|1<<v)
    for s in ([start] if start is not None else range(len(A))): walk((s,),1<<s)
    return out

def cycles(A):
    out=[]
    for s in range(len(A)):
        def walk(p,used):
            if len(p)>=3 and len(p)%2 and A[p[-1]]>>s&1: out.append(p)
            for v in range(s+1,len(A)):
                if not used>>v&1 and A[p[-1]]>>v&1: walk(p+(v,),used|1<<v)
        walk((s,),1<<s)
    return out

def cmask(c): return sum(1<<v for v in c)
def pack(cyc):
    ans=Counter()
    def rec(i,used,k):
        if i==len(cyc): ans[k,used.bit_count()]+=1; return
        rec(i+1,used,k)
        m=cmask(cyc[i])
        if not m&used: rec(i+1,used|m,k+1)
    rec(0,0,0)
    return ans

def mul(P,Q):
    R=Counter()
    for (k,s),v in P.items():
        for (l,t),w in Q.items(): R[k+l,s+t]+=v*w
    return R

def predicted(A,B,a,b):
    CA,CB=cycles(A),cycles(B); P=mul(pack(CA),pack(CB)); new=[]
    for p in paths(A,start=a):
        for q in paths(B,end=b):
            if (len(p)+len(q))%2:
                new.append(p+tuple(v+len(A) for v in q))
                R=mul(pack([c for c in CA if not cmask(c)&cmask(p)]),pack([c for c in CB if not cmask(c)&cmask(q)]))
                for (k,s),v in R.items(): P[k+1,s+len(p)+len(q)]+=v
    return P,new

def h_subsets(A):
    n=len(A); H=[0]*(1<<n); H[0]=1; dp={}
    for mask in range(1,1<<n):
        row=[0]*n
        for v in range(n):
            if mask>>v&1:
                prev=mask^(1<<v)
                row[v]=1 if not prev else sum(dp[prev][u] for u in range(n) if prev>>u&1 and A[u]>>v&1)
        dp[mask]=row; H[mask]=sum(row)
    return H

def q_direct(A,r):
    H=h_subsets(A); size=len(H); f=[0]*size; f[0]=1
    for _ in range(r):
        g=[0]*size
        for mask in range(size):
            sub=mask
            while True:
                g[mask]+=f[mask^sub]*H[sub]
                if sub==0: break
                sub=(sub-1)&mask
        f=g
    return f[-1]

def q_poly(P,n,r): return sum(v*2**k*r**(n-s+k) for (k,s),v in P.items())
def poly_json(P): return {f'{k},{s}':v for (k,s),v in sorted(P.items())}

def main():
    words=[''.join(w) for n in range(6) for w in itertools.product('01',repeat=n)]
    clique=[]
    for w in words:
        A=apex(w); C=cycles(A)
        p=sum((1 if j-i==1 else 2**(j-i-2)) for i in range(len(w)) for j in range(i+1,len(w)) if w[i]=='1' and w[j]=='0')
        assert p==len(C) and all(0 in c for c in C)
        assert h_subsets(A)[-1]==1+2*p
        if w in ['10','100','1000','1010','1100']: clique.append({'word':w,'addresses':p,'H':1+2*p})
    small=[w for w in words if len(w)<=3]; trials=0
    for u in small:
        for v in small:
            A,B=apex(u),apex(v); T=join(A,B,0,0)
            pred,new=predicted(A,B,0,0); direct=pack(cycles(T))
            assert pred==direct
            assert len(cycles(T))==len(cycles(A))+len(cycles(B))+len(new)
            for r in (1,2,3): assert q_direct(T,r)==q_poly(pred,len(T),r)
            trials+=1
    A=apex('10'); T=join(A,A,0,0); C=cycles(T); P=pack(C)
    assert P==Counter({(0,0):1,(1,3):4,(1,5):2,(2,6):1})
    nonedges=[(i,j) for i in range(len(C)) for j in range(i+1,len(C)) if not cmask(C[i])&cmask(C[j])]
    assert len(nonedges)==1
    # Literal signed-address decoder: each permitted support receives all 2^|S| signs.
    decoder_states=0
    for bits in itertools.product((0,1),repeat=len(C)):
        chosen=[C[i] for i,x in enumerate(bits) if x]
        if sum(len(c) for c in chosen)==cmask(tuple(set(v for c in chosen for v in c))).bit_count(): decoder_states+=2**sum(bits)
    assert decoder_states==17
    # Nested return joins, with nonempty residual cycles in the recurrence.
    nested=[]
    for left,right,a,b in [(T,A,0,0),(T,A,5,1),(A,T,1,5)]:
        U=join(left,right,a,b); pred,_=predicted(left,right,a,b); direct=pack(cycles(U)); assert pred==direct
        assert q_poly(pred,len(U),1)==h_subsets(U)[-1]
        assert q_poly(pred,len(U),2)==q_direct(U,2)
        nested.append({'rows':U,'packing':poly_json(pred),'H':q_poly(pred,len(U),1)})
    receipt={'status':'FINITE-EXACT','apex_words_lengths_0_through_5':len(words),'all_return_pairs_word_lengths_0_through_3':trials,'cover_r_values':[1,2,3],'clique_initials':clique,'six_address_gate':{'rows':T,'cycles':C,'only_compatible_pair':nonedges,'packing':poly_json(P),'Q_at_1_2_3':[q_poly(P,6,r) for r in [1,2,3]],'decoded_states':decoder_states},'nested_return_controls':nested}
    Path(__file__).with_name('selector_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__': main()
