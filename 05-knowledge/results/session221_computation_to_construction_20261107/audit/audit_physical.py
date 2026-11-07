#!/usr/bin/env python3
"""Independent physical audit: backwards endpoint DP, direct sorting steps,
and arithmetic formulas reconstructed from statements rather than producer code."""
from functools import lru_cache
from itertools import product,combinations
import json,random

def hp(out):
    n=len(out)
    if not n:return 1
    incoming=[sum(1<<u for u in range(n) if out[u]>>v&1) for v in range(n)]
    @lru_cache(None)
    def paths(mask,last):
        rest=mask^(1<<last)
        if not rest:return 1
        prev=rest&incoming[last];ans=0
        while prev:
            b=prev&-prev;prev-=b
            ans+=paths(rest,b.bit_length()-1)
        return ans
    return sum(paths((1<<n)-1,i) for i in range(n))

def graph(bits):
    # Place the apex last, unlike the construction's apex-first labels.
    n=len(bits);out=[sum(1<<j for j in range(i+1,n)) for i in range(n)]+[0]
    for i,b in enumerate(bits):
        if b:out[i]|=1<<n
        else:out[n]|=1<<i
    return out

def w(d):return 1 if d<3 else 2**(d-2)
def h_formula(s):return 1+2*sum(w(j-i) for i in range(len(s)) for j in range(i+1,len(s)) if s[i]==0 and s[j]==1)
def response(s,i):return 1+sum(2**(i-j-2) for j in range(i-1) if s[j]==0)+sum(2**(j-i-3) for j in range(i+3,len(s)) if s[j]==1)
def swap(s,S):
    v=list(s)
    for i in S:v[i],v[i+1]=v[i+1],v[i]
    return tuple(v)
def f(n):
    if n==0:return 0
    return f(n//2)+(n//2+n//2%2 if n%2 else 0)

def main():
    physical=local=parallel=0;ranges={}
    for n in range(10):
        for s in product((0,1),repeat=n):
            h=hp(graph(s));physical+=1
            if h!=h_formula(s):raise ValueError(('physical formula',s,h))
            sites=[i for i in range(n-1) if s[i:i+2]==(1,0)]
            for i in sites:
                d=hp(graph(swap(s,[i])))-h;local+=1
                if d!=2*response(s,i):raise ValueError(('local',s,i,d))
                ranges.setdefault((n,i),set()).add(d//2)
            if sites:
                d=hp(graph(swap(s,sites)))-h;parallel+=1
                expected=2*(sum(response(s,i) for i in sites)+sum(w(j-i-1) for i,j in combinations(sites,2)))
                if d!=expected:raise ValueError(('parallel',s,d,expected))
    for (n,i),vals in ranges.items():
        a=max(i-1,0);b=max(n-i-3,0)
        if vals!=set(range(1,2**a+2**b)):raise ValueError(('range',n,i,vals))
    for n in range(3,10):
        maxgap=max(max(vals) for (nn,i),vals in ranges.items() if nn==n)
        if maxgap!=2**(n-3):raise ValueError(('maximum',n,maxgap))
    examples=[]
    for N in range(1,129):
        m=N-1;digits=() if not m else tuple(int(c) for c in bin(m)[2:]);p=tuple(1-c for c in digits)+(1,)
        a,b=p+(1,0),p+(0,1);ha,hb=hp(graph(a)),hp(graph(b))
        if (ha,hb)!=(4*N-3-f(m),6*N-3-f(m)):raise ValueError(('selector',N,ha,hb))
        edges=sum((x^y).bit_count() for x,y in zip(graph(a),graph(b)))//2
        if edges!=2:raise ValueError(('edges',N,edges))
        if N in (1,2,7,21,47,128):examples.append({'N':N,'H_pair':[ha,hb],'vertices':len(a)+1})
    receipt={'status':'PASS','independent_reverse_endpoint_DP':physical,'local_physical_comparisons':local,'parallel_physical_comparisons':parallel,'fixed_length_site_ranges_checked':len(ranges),'selectors_physical_DP':128,'examples':examples,'boundaries':['Both tournaments share labels and differ in exactly two apex arcs; no isomorphism-distance claim.','Count differences cover positive even integers, not absolute H targets.','Optimal order is proved only among adjacent 10->01 one-apex swaps.','N=1 general recipe uses four vertices; 10->01 gives optimal three.']}
    open('physical_receipt.json','w').write(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
