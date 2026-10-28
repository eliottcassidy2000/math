#!/usr/bin/env python3
"""Complete even-target fiber counts using even sources and odd squares."""
import json, math
from collections import Counter
from check_aliquot_fibers import need,sieve,fiber,DATA_DIR

def even_fibers(V):
    N=2*(V-1)
    # Exact ordinary sigma on the even-source interval.
    sig=[0]*(N+1)
    for d in range(1,N+1):
        for n in range(d,N+1,d): sig[n]+=d
    F=[[] for _ in range(V//2+1)]
    for n in range(2,N+1,2):
        t=sig[n]-n
        if 2<=t<=V and t%2==0: F[t//2].append(n)
    # Smallest prime factors of roots, NOT a sieve to V^2.
    spf=list(range(V))
    for p in range(2,math.isqrt(V-1)+1):
        if spf[p]==p:
            for n in range(p*p,V,p):
                if spf[n]==n: spf[n]=p
    for m in range(3,V,2):
        z=m; ss=1
        while z>1:
            p=spf[z]; e=0
            while z%p==0: z//=p; e+=1
            ss*= (p**(2*e+1)-1)//(p-1)
        t=ss-m*m
        if 2<=t<=V and t%2==0: F[t//2].append(m*m)
    for row in F: row.sort()
    return F

if __name__=='__main__':
    V=1001000
    # Full source scan via inherited (t-1)^2 bound, independently small.
    small=120
    s=sieve((small-1)**2); F=even_fibers(small)
    for t in range(2,small+1,2):
        need(F[t//2]==[n for n in range(1,(t-1)**2+1) if s[n]==t],f'small {t}')
    F=even_fibers(V)
    with (DATA_DIR / 'results.json').open() as h: targets=json.load(h)
    need(F[1001000//2]==targets['1001000']['fiber'],'large target comparison')
    empty=[t for t in range(2,V+1,2) if not F[t//2]]
    expected=[2,52,88,96,120,124,146,162,188,206,210,216,238,246,248,262,268,276,288,290,292,304,306,322,324,326,336,342,372,406]
    need(empty[:len(expected)]==expected,'published initial even untouchables')
    summary={'status':'FINITE-EXACT','even_target_range':[2,V], 'target_count':V//2,'untouchable_count':len(empty),'initial_even_untouchables':empty[:50], 'fiber_cardinality_histogram':dict(sorted(Counter(map(len,F[1:])).items())), 'source_count':sum(map(len,F)), 'maximum_fiber_cardinality':max(map(len,F)), 'maximal_fiber_targets':[],'complete_target1001000':F[500500], 'small_independent_targets':[2,small]}
    for t in [52,124,560554,640634,895994,1000000]:
        need(F[t//2]==fiber(t)[0],f'independent target DFS {t}')
    summary['independent_DFS_targets']=[52,124,560554,640634,895994,1000000]
    published=[13863,28572,43515,58459,73565,88828,104062,119302,134758,150232]
    for x,c in zip(range(100000,1000001,100000),published):
        need(sum(t<=x for t in empty)+1==c,f'published count {x}')
    summary['published_count_check']='PASS: even count + known odd untouchable 5 matches all ten primary table rows through 10^6'
    summary['even_untouchable_counts_at_100k']={str(x):sum(t<=x for t in empty) for x in range(100000,1000001,100000)}
    maxsize=max(map(len,F)); summary['maximal_fiber_targets']=[2*i for i in range(1,len(F)) if len(F[i])==maxsize]
    with (DATA_DIR / 'even_range_results.json').open('w') as h: json.dump(summary,h,indent=2,sort_keys=True)
    print(json.dumps(summary,indent=2,sort_keys=True))
