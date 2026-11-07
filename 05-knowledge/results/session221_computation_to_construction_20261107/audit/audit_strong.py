#!/usr/bin/env python3
import sys,json,random
from itertools import product
from functools import lru_cache
from audit_physical import hp,graph,swap,f

def strong(out):
    n=len(out)
    if not n:return False
    for v in range(n):
        seen=1<<v;old=-1
        while old!=seen:
            old=seen
            for j in range(n):
                if seen>>j&1:seen|=out[j]
        if seen!=(1<<n)-1:return False
    return True

def main():
    counts={};sranges={};strongtested=0;singlearc=0
    for n in range(4,10):
        for mid in product((0,1),repeat=n-2):
            s=(0,)+mid+(1,)
            for i in range(1,n-2):
                if s[i:i+2]!=(1,0):continue
                t=swap(s,[i]);before,after=graph(s),graph(t)
                z=list(before);z[i]^=1<<(i+1);z[i+1]^=1<<i
                perm=list(range(n+1));perm[i],perm[i+1]=perm[i+1],perm[i]
                relabel=[sum(1<<perm[j] for j in range(n+1) if z[v]>>j&1) for v in range(n+1)]
                relabel=[relabel[perm[v]] for v in range(n+1)]
                if relabel!=after:raise ValueError('singlearc conjugacy')
                singlearc+=1
                if not strong(z) or not strong(before):raise ValueError('strong')
                N=(hp(after)-hp(before))//2;sranges.setdefault((n,i),set()).add(N);strongtested+=1
    for (n,i),vals in sranges.items():
        aa=i-1;bb=n-i-3;lo=1+(2**(aa-1) if aa else 0)+(2**(bb-1) if bb else 0);hi=2**aa+2**bb-1
        if vals!=set(range(lo,hi+1)):raise ValueError(('strong range',n,i,vals,lo,hi))
    for N in range(2,65):
        m=N-1;p=tuple(1-int(c) for c in bin(m)[2:])+(1,)
        A,B=p+(1,0,1),p+(0,1,1)
        if (hp(graph(A)),hp(graph(B)))!=(12*N-9-f(m),14*N-9-f(m)):raise ValueError('strong selector')
    counts.update(strong_physical_moves=strongtested,single_arc_permutation_checks=singlearc,strong_site_ranges=len(sranges),strong_selectors=63)
    res={'status':'PASS','counts':counts,'scope':'Independent literal strong single-arc conjugacies, fixed-site ranges and explicit selectors.'}
    open('strong_receipt.json','w').write(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
if __name__=='__main__':main()
