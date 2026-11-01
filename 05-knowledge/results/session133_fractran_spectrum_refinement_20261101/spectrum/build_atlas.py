#!/usr/bin/env python3
"""Standalone exact tournament class atlas through order six; stdlib only.
Code bit k for pair (u,v), u<v lexicographic: 1 iff u->v.
Canonical representative is minimum code among all vertex permutations.
"""
import json, math
from collections import Counter, defaultdict
from itertools import permutations, combinations
from pathlib import Path

def need(ok,msg):
    if not ok: raise ValueError(msg)

def pairs(n):return list(combinations(range(n),2))
def rows_of(n,code):
    a=[0]*n
    for k,(u,v) in enumerate(pairs(n)):
        if code>>k&1:a[u]|=1<<v
        else:a[v]|=1<<u
    return a

def code_of(a):return sum(((a[u]>>v)&1)<<k for k,(u,v) in enumerate(pairs(len(a))))

def perm_code(a,p):
    return sum(((a[p[u]]>>p[v])&1)<<k for k,(u,v) in enumerate(pairs(len(a))))

def canonical(a):return min(perm_code(a,p) for p in permutations(range(len(a))))

def path_dp(a):
    n=len(a)
    if not n:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if dp[mask][v]:
                remaining=a[v]&~mask
                while remaining:
                    b=remaining&-remaining;remaining-=b
                    dp[mask|b][b.bit_length()-1]+=dp[mask][v]
    return sum(dp[-1])

def direct_paths(a):
    return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(a)-1)) for p in permutations(range(len(a))))

def components(a):
    n=len(a); reach=[]
    for v in range(n):
        r=1<<v;old=-1
        while old!=r:
            old=r
            for w in range(n):
                if r>>w&1:r|=a[w]
        reach.append(r)
    left=set(range(n));out=[]
    while left:
        v=min(left); c=sorted(w for w in left if (reach[v]>>w&1) and (reach[w]>>v&1))
        out.append(c);left-=set(c)
    # Unique dominance order, not label order.
    return sorted(out,key=lambda c:-sum(bool(a[c[0]]>>d[0]&1) for d in out if d!=c))

def odd_cycles(a):
    out=[]
    for k in range(3,len(a)+1,2):
        for s in combinations(range(len(a)),k):
            for tail in permutations(s[1:]):
                p=(s[0],)+tail
                if all(a[p[i]]>>p[(i+1)%k]&1 for i in range(k)):
                    out.append(p)
    return out

def packing(a):
    cs=odd_cycles(a); counts=[1]+[0]*(len(a)//3)
    masks=[sum(1<<v for v in c) for c in cs]
    def visit(start,used,k):
        for i in range(start,len(masks)):
            if not used&masks[i]:
                counts[k+1]+=1;visit(i+1,used|masks[i],k+1)
    visit(0,0,0)
    return counts,dict(sorted(Counter(map(len,cs)).items()))

def main():
    old={0}; all_rows=[];counts=[]
    for n in range(1,7):
        current=set()
        for code in old:
            base=rows_of(n-1,code)
            for mask in range(1<<(n-1)):
                a=base+[0]
                for u in range(n-1):
                    if mask>>u&1:a[u]|=1<<(n-1)
                    else:a[-1]|=1<<u
                current.add(canonical(a))
        labelled=0;strong_hist=Counter();strong_low=[]
        for code in sorted(current):
            a=rows_of(n,code);h=path_dp(a)
            need(h==direct_paths(a),'DP/permutation disagreement')
            alpha,lengths=packing(a)
            need(h==sum((1<<k)*x for k,x in enumerate(alpha)),'OCF disagreement')
            aut=sum(perm_code(a,p)==code for p in permutations(range(n)))
            need(math.factorial(n)%aut==0,'nonintegral orbit')
            labelled+=math.factorial(n)//aut
            cs=components(a); strong=len(cs)==1 and n>=3
            row={'n':n,'code':code,'rows':a,'H':h,'automorphisms':aut,'strong':strong,
                 'self_converse':canonical([((1<<n)-1) ^ (1<<v) ^ a[v] for v in range(n)])==code,
                 'SCC_sizes':[len(c) for c in cs],'alpha':alpha,'cycles_by_length':lengths}
            all_rows.append(row)
            if strong:
                strong_hist[h]+=1
                if h<=21:strong_low.append(row)
        need(labelled==1<<(n*(n-1)//2),'augmentation cover/orbit mass failure')
        counts.append({'n':n,'classes':len(current),'labelled_mass':labelled,
                       'strong_H_multiplicities':dict(sorted(strong_hist.items())),
                       'strong_low_classes':len(strong_low)})
        old=current
    need([r['classes'] for r in counts]==[1,1,2,4,12,56],'known small class count control')
    # Exact all-order input: cited Busch gives strong order <=6 if H<=21.
    need(5**6>21**3,'Busch cutoff arithmetic')
    strong_values={r['H'] for r in all_rows if r['strong']}
    need(7 not in strong_values and 21 not in strong_values,'forbidden strong base failure')
    # Hostile: same H can have different packing profile or SCC count.
    h15=[r for r in all_rows if r['H']==15]
    need(len({tuple(r['alpha']) for r in h15})>1,'lost-packing-coordinate hostile missing')
    payload={'schema':'session133.small-tournament-atlas.v1',
        'scope':'Complete isomorphism classes through order 6; all-order claims require cited Busch cutoff and SCC product.',
        'orientation':'lexicographic upper triangle, bit 1 means lower-labelled vertex points to higher',
        'universe':counts,'classes':all_rows,
        'strong_low_multiplicities':{
            str(h):{str(n):sum(r['n']==n and r['H']==h and r['strong'] for r in all_rows)
                    for n in range(3,7) if any(r['n']==n and r['H']==h and r['strong'] for r in all_rows)}
            for h in sorted(strong_values) if h<=21},
        'strong_self_converse_low_multiplicities':{
            str(h):{str(n):sum(r['n']==n and r['H']==h and r['strong'] and r['self_converse'] for r in all_rows)
                    for n in range(3,7) if any(r['n']==n and r['H']==h and r['strong'] and r['self_converse'] for r in all_rows)}
            for h in sorted(strong_values) if h<=21},
        'controls':{'DP_equals_permutations':'all 76 classes, orders 1..6',
                    'DP_equals_literal_OCF':'all 76 classes, orders 1..6',
                    'labelled_orbit_mass':'all orders 1..6',
                    'H_not_unique_packing_profile':'H=15'}}
    print(json.dumps(payload,indent=2,sort_keys=True))
if __name__=='__main__':main()
