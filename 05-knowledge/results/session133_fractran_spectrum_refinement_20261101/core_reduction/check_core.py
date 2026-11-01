#!/usr/bin/env python3
"""Exact, dependency-free checks for the singleton-SCC core floor.
No asserts: checks also execute with python -O.
"""
import argparse, hashlib, itertools, json, math
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def strong_floor(n):
    if n < 3:
        raise ValueError(n)
    q, r = divmod(n, 3)
    return (3 * 5 ** (q-1), 5 ** q, 9 * 5 ** (q-1))[r]


def core_floor(n):
    if n == 0: return 1
    if n in (1,2): return None
    if n == 3: return 3
    if n == 4: return 5
    if n == 5: return 9
    q,r=divmod(n,3)
    return (3 ** q, 5 * 3 ** (q-1), 25 * 3 ** (q-2))[r]


def budget(M, restrict=True):
    if M < 1: raise ValueError(M)
    best=(0,0,0,1)
    b=0; p5=1
    while p5 <= M and (not restrict or b <= 2):
        a=0; cost=p5
        while 3*cost <= M:
            cost*=3; a+=1
        candidate=(3*a+4*b,a,b,cost)
        if candidate > best: best=candidate
        b+=1; p5*=5
    return best


def tournament(n, code):
    rows=[0]*n; bit=0
    for i in range(n):
        for j in range(i+1,n):
            if (code>>bit)&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
            bit+=1
    return tuple(rows)


def path_counts(rows):
    n=len(rows); allmask=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    h=[0]*(1<<n); h[0]=1
    for v in range(n): dp[1<<v][v]=1
    for s in range(1,1<<n):
        for v,value in enumerate(dp[s]):
            if not value: continue
            h[s]+=value
            rest=rows[v]& (allmask^s)
            while rest:
                bit=rest&-rest; rest-=bit
                dp[s|bit][bit.bit_length()-1]+=value
    return h


def literal_h(rows):
    return sum(all(rows[a]>>b&1 for a,b in zip(p,p[1:]))
               for p in itertools.permutations(range(len(rows))))


def core_mask(rows):
    n=len(rows); reach=[rows[v]|(1<<v) for v in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1: reach[i]|=reach[k]
    return sum(1<<v for v in range(n)
               if any(v!=w and reach[v]>>w&1 and reach[w]>>v&1
                      for w in range(n)))


def triangle_support(rows):
    out=0
    for a,b,c in itertools.combinations(range(len(rows)),3):
        if ((rows[a]>>b&1) and (rows[b]>>c&1) and (rows[c]>>a&1)) or \
           ((rows[a]>>c&1) and (rows[c]>>b&1) and (rows[b]>>a&1)):
            out|=(1<<a)|(1<<b)|(1<<c)
    return out


def join(blocks):
    rows=[]; offset=0; total=sum(map(len,blocks))
    for block in blocks:
        right=((1<<total)-1)^((1<<(offset+len(block)))-1)
        rows.extend((row<<offset)|right for row in block)
        offset+=len(block)
    return tuple(rows)


def near_transitive(n):
    rows=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    rows[0]^=1<<(n-1); rows[n-1]|=1
    return tuple(rows)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out', default='normal.json')
    args=parser.parse_args()
    # Independent component-partition dynamic program using the published f(n).
    max_n=400
    dp=[None]*(max_n+1); dp[0]=1
    for n in range(3,max_n+1):
        values=[strong_floor(s)*dp[n-s] for s in range(3,n+1)
                if dp[n-s] is not None]
        dp[n]=min(values) if values else None
    for n in range(max_n+1):
        require(dp[n]==core_floor(n), f'core floor partition DP at {n}')
    for n in range(6,max_n-2):
        require(core_floor(n+3)==3*core_floor(n),f'growth at {n}')
    for M in range(1,100001):
        got=budget(M)
        require(got[0]==budget(M,False)[0], f'budget unrestricted at {M}')
        feasible=[n for n in range(40) if core_floor(n) is not None and core_floor(n)<=M]
        require(got[0]==max(feasible),f'budget inverse floor at {M}')
    # Every labelled tournament through order six; all induced subgraphs share DP.
    census=[]; deletion_tests=0; literal_tests=0
    for n in range(7):
        total=1<<(n*(n-1)//2); full=(1<<n)-1
        counts={}; reduced_min=None; core_hist={}; minimizing_core_count=0
        for code in range(total):
            rows=tournament(n,code); hs=path_counts(rows); H=hs[full]
            if n<=5:
                require(literal_h(rows)==H,f'permutation DP at {(n,code)}')
                literal_tests+=1
            core=core_mask(rows)
            require(core==triangle_support(rows),f'core triangle support {(n,code)}')
            require(hs[core]==H,f'core preserves H {(n,code)}')
            require(core.bit_count()<=budget(H)[0],f'budget violated {(n,code,H)}')
            for s in range(1<<n):
                require((hs[s]==H)==(s&core==core),f'canonical induced core {(n,code,s)}')
                deletion_tests+=1
            counts[H]=counts.get(H,0)+1
            if core==full:
                reduced_min=H if reduced_min is None else min(H,reduced_min)
                if H==core_floor(n): minimizing_core_count+=1
            if H<=21:
                core_hist[core.bit_count()]=core_hist.get(core.bit_count(),0)+1
        require(reduced_min==core_floor(n),f'actual core min at {n}')
        census.append({'n':n,'labelled':total,'H_spectrum':sorted(counts),
                       'min_singleton_free_H':reduced_min,'labelled_minimizing_cores':minimizing_core_count,
                       'low21_core_hist':core_hist})
    require(census[6]['labelled_minimizing_cores']==80,'order6 labelled extremizer count')
    # Independent size/weight word DP checks the extremizer number formulas.
    extremizer_counts=[]
    for N in range(6,101):
        b=N%3; a=(N-4*b)//3
        e=math.comb(a+b,b)
        weights={}
        def visit(n,cost,word):
            if n==N:
                weights[cost]=weights.get(cost,0)+1
                return
            if n+3<=N: visit(n+3,cost*3,word+'3')
            if n+4<=N: visit(n+4,cost*5,word+'4')
        if N<=25:
            visit(0,1,'')
            require(min(weights)==core_floor(N),f'word floor {N}')
            require(weights[core_floor(N)]==e,f'word minimizer count {N}')
        extremizer_counts.append({'N':N,'unlabelled':e,'labelled':e*math.factorial(N)//3**a})
    # Witnesses include all newly sharpened low thresholds and hostile twins.
    C3=near_transitive(3); S4=near_transitive(4)
    witnesses=[]
    for M in [1,3,5,7,9,15,21,25,27,45,75,81,999999]:
        N,a,b,cost=budget(M)
        rows=join([C3]*a+[S4]*b)
        require(len(rows)==N,f'witness order {M}')
        require(core_mask(rows)==(1<<N)-1,f'witness core {M}')
        if N<=15:
            require(path_counts(rows)[-1]==cost,f'witness H {M}')
        witnesses.append({'M':M,'K':N,'triangles':a,'strong4':b,'H':cost,
                          'rows':rows,'direct_DP_checked':N<=15})
    twin_hostiles=[]
    for t in range(1,11):
        rows=near_transitive(t+2)
        H=path_counts(rows)[-1]
        require(H==2**t+1,f'near-transitive formula t={t}')
        # Middle t vertices are a transitive module; contracting them gives C3.
        middle=list(range(1,t+1))
        outside_module=(1<<0)|(1<<(t+1))
        for u,v in itertools.combinations(middle,2):
            require(rows[u]&outside_module==rows[v]&outside_module,f'module {t,u,v}')
        for u,v in zip(middle,middle[1:]):
            others=((1<<(t+2))-1)^((1<<u)|(1<<v))
            require(rows[u]&others==rows[v]&others,f'adjacent twins {t,u,v}')
        twin_hostiles.append({'module_size':t,'H':H,'contracted_H':3})
    result={'status':'PASS','scope':'all-order theorem proved in THEORY.md; finite independent checks here',
            'partition_DP_orders':max_n,'budget_checked_M':100000,
            'labelled_tournaments':sum(x['labelled'] for x in census),
            'induced_equality_tests':deletion_tests,'literal_permutation_tests':literal_tests,
            'census':census,'budget_witnesses':witnesses,'twin_contraction_hostiles':twin_hostiles,
            'core_floor_n0_through24':[core_floor(n) for n in range(25)],
            'extremizer_count_formulas_n6_through100':extremizer_counts}
    Path(args.out).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','out':args.out,'labelled_tournaments':result['labelled_tournaments'],
                      'induced_equality_tests':deletion_tests},sort_keys=True))

if __name__=='__main__':main()
