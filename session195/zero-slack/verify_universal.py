#!/usr/bin/env python3
from pathlib import Path
from itertools import product
import hashlib,json
import verify as v
HERE=Path(__file__).resolve().parent

def components(a):
    n=len(a);reach=[]
    for i in range(n):
        m=1<<i
        while True:
            nxt=m
            for j in v.bits(m):nxt|=a[j]
            if nxt==m:break
            m=nxt
        reach.append(m)
    rest=(1<<n)-1;parts=[]
    while rest:
        r=next(v.bits(rest));part=sum(1<<j for j in v.bits(rest) if reach[r]>>j&1 and reach[j]>>r&1)
        parts.append(part);rest^=part
    return sorted(parts,key=lambda p:-sum(1 for q in parts if p!=q and a[next(v.bits(p))]>>next(v.bits(q))&1))

def check(a,r):
    H,E,S=v.tables(a);h,e,R,d=v.profile(a,r,(H,E,S));parts=components(a)
    j=next(j for j,m in enumerate(parts) if m>>r&1)
    before=after=1;rootR=rootE=None
    for i,m in enumerate(parts):
        labels=list(v.bits(m));b=v.subgraph(a,labels);tab=v.tables(b);hb=tab[0][-1]
        if i<j:before*=hb
        elif i>j:
            full=(1<<len(b))-1
            after*=sum(tab[0][u]*tab[0][full^u] for u in range(full+1))
        else:
            rr=labels.index(r);rootR=v.profile(b,rr,tab)[2];rootE=tab[1][-1][rr]
    v.need(R==before*rootR*after and e==(before*rootE if j==len(parts)-1 else 0),'SCC response identity')
    kept=0;t=0
    for i,m in enumerate(parts):
        if m.bit_count()>1 or i==j:kept|=m
        elif i>j:t+=1
    labels=list(v.bits(kept));C=v.subgraph(a,labels);rr=labels.index(r);hc,ec,Rc,dc=v.profile(C,rr)
    v.need(hc==h and R==(2**t)*Rc and e==(ec if t==0 else 0),'core transport identity')
    v.need(len(C)<=h+1,'finite core order bound')
    # Original singleton birth versus the canonical core plus terminal padding.
    old,oldr=v.unary(a,r,'A');new=C;mark=rr
    for _ in range(t):new,mark=v.unary(new,mark,'D')
    new,newr=v.unary(new,mark,'A')
    expected=h+(2**t)*Rc-(ec if t==0 else 0)
    v.need(v.profile(old,oldr)[0]==v.profile(new,newr)[0]==expected,'constructive finite-core reconstruction')
    return t,len(a),len(C)

def main():
    tests=graphs=0;maximum_deleted_padding=0
    for n in range(1,6):
        for code in range(1<<(n*(n-1)//2)):
            a=v.tournament(n,code);graphs+=1
            for r in range(n):
                t,_,_=check(a,r);maximum_deleted_padding=max(maximum_deleted_padding,t);tests+=1
    padded=0
    # Arbitrary singleton placements around a nontrivial core, with every mark.
    core=(2,4,1)
    for left,right in product(range(4),repeat=2):
        sizes=(left,3,right);n=sum(sizes);a=[0]*n
        for i in range(n):
            for j in range(i+1,n):
                if left<=i<left+3 and left<=j<left+3:
                    if core[i-left]>>(j-left)&1:a[i]|=1<<j
                    else:a[j]|=1<<i
                else:a[i]|=1<<j
        for r in range(n):check(tuple(a),r);padded+=1
    receipt=dict(status='PASS',labelled_graphs=graphs,rooted_checks=tests,padded_rooted_checks=padded,
                 maximum_deleted_padding_in_exhaustive_universe=maximum_deleted_padding,
                 scope='Exact SCC response, finite-core order bound, and reconstruction controls; no spectrum inference from the finite universe.',
                 script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 helper_sha256=hashlib.sha256((HERE/'verify.py').read_bytes()).hexdigest())
    rendered=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(HERE/'universal-receipt.json').write_text(rendered);print(rendered,end='')
if __name__=='__main__':main()
