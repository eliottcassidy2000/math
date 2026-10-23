#!/usr/bin/env python3
"""Independent Stirling-inversion audit using direct substituted tournaments."""
import itertools, math, json
from pathlib import Path
from functools import lru_cache
from fractions import Fraction
from verify_context import tournament, transitive, profile, q_from_profile, h_subsets, check

OUT=Path(__file__).resolve().parent

def substitute(a,sizes):
    starts=[0]
    for m in sizes:starts.append(starts[-1]+m)
    out=[0]*starts[-1]
    for i,m in enumerate(sizes):
        for u in range(m):
            out[starts[i]+u]=sum(1<<(starts[i]+v) for v in range(u+1,m))
            for j,k in enumerate(sizes):
                if a[i]>>j&1: out[starts[i]+u]|=((1<<k)-1)<<starts[j]
    return out

@lru_cache(None)
def s1(n,k):
    if n==0:return int(k==0)
    if k<1 or k>n:return 0
    return s1(n-1,k-1)-(n-1)*s1(n-1,k)

def walk_count(a,content):
    q=len(a)
    @lru_cache(None)
    def ending(c,last):
        if c[last]==0:return 0
        if sum(c)==1:return 1
        reduced=list(c);reduced[last]-=1;reduced=tuple(reduced)
        return sum(ending(reduced,u) for u in range(q) if a[u]>>last&1)
    return sum(ending(tuple(content),last) for last in range(q))

def main():
    stats={}
    for q,bound in [(2,3),(3,3),(4,2)]:
        instances=0
        for mask in range(1<<(q*(q-1)//2)):
            a=tournament(q,mask)
            sizes=list(itertools.product(range(1,bound+1),repeat=q))
            h={m:h_subsets(substitute(a,m))[-1] for m in sizes}
            for c in sizes:
                value=sum(h[m]*math.prod(s1(ci,mi) for ci,mi in zip(c,m)) for m in itertools.product(*(range(1,ci+1) for ci in c)))
                denominator=math.prod(math.factorial(ci) for ci in c)
                expected=walk_count(a,c)
                check(value==denominator*expected,('Stirling inversion',q,mask,c,value,expected))
                instances+=1
        stats[q]={'bound':bound,'quotients':1<<(q*(q-1)//2),'cells':instances}
        print('inversion',q,stats[q],flush=True)
    twins=[]
    for mask in (408,1332):
        a=tournament(6,mask);p=profile(a)
        value=h_subsets(substitute(a,(2,)*6))[-1]
        twins.append({'mask':mask,'profile':p,'Q':q_from_profile(p),'uniform_lift_H':value})
    check(twins[0]['profile']==twins[1]['profile'],'twins profile equality')
    check([x['uniform_lift_H'] for x in twins]==[421425,411513],'twins differ as quotients')
    stats['twins']=twins
    print('twins',twins,flush=True)
    (OUT/'quotient_verification.json').write_text(json.dumps(stats,indent=2)+'\n')
    print('ALL QUOTIENT CHECKS PASSED')

if __name__=='__main__':main()
