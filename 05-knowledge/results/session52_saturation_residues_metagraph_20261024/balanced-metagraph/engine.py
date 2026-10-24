#!/usr/bin/env python3
"""Exact n-only balanced fixed-path port counts. No tournament classes are input."""
from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def balanced_endpoint_dp(n):
    """Return LL,LH,HL,HH counts and cache states by forgetting finished rows.

    A state records wins against already processed vertices, in the surviving
    path order. The next path arc is mandatory; all other current-row arcs
    are chosen. Last vertex is classified low/high at the terminal step.
    """
    need(n >= 2 and n % 2 == 0, 'positive even order at least two')
    p=n//2
    @lru_cache(None)
    def finish(wins):
        r=len(wins)
        if r == 1:
            return (int(wins[0] == p-1), int(wins[0] == p))
        answer=[0,0]
        # current vertex beats the next one; the next receives no past win.
        free=r-2
        for target in (p-1,p):
            k=target-wins[0]-1
            if not 0 <= k <= free:
                continue
            for chosen in combinations(range(free),k):
                selected=set(chosen)
                new=(wins[1],)+tuple(wins[j+2]+int(j not in selected) for j in range(free))
                # Every future vertex can win at most r-2 still-unseen games.
                if any(a>p or a+(r-2)<p-1 for a in new):
                    continue
                add=finish(new)
                answer[0]+=add[0]; answer[1]+=add[1]
        return tuple(answer)
    result={}
    for label,target in [('L',p-1),('H',p)]:
        total=[0,0]
        k=target-1
        if 0 <= k <= n-2:
            for chosen in combinations(range(n-2),k):
                selected=set(chosen)
                wins=(0,)+tuple(int(j not in selected) for j in range(n-2))
                a,b=finish(wins); total[0]+=a; total[1]+=b
        result[label+'L'],result[label+'H']=total
    return result, finish.cache_info().currsize

def blue_half_polynomial(p):
    """G_p in p degree variables, truncated at degree p per variable."""
    poly={(1,)*p:1}
    factors=[]
    def vec(indices):
        x=[0]*p
        for i in indices:x[i]+=1
        return tuple(x)
    zero=(0,)*p
    for i in range(p):
        for j in range(i+2,p):factors.append((vec((i,)),vec((j,))))
    for i in range(p):
        for j in range(i+1,p):factors.append((zero,vec((i,j))))
    for i in range(p-1):factors.append((zero,vec((i,))))
    peak=len(poly)
    for factor in factors:
        nxt=Counter()
        for d,c in poly.items():
            for monomial in factor:
                e=tuple(a+b for a,b in zip(d,monomial))
                if all(a<=p for a in e):nxt[e]+=c
        poly=nxt;peak=max(peak,len(poly))
    balanced={d:c for d,c in poly.items() if all(a in (p-1,p) for a in d)}
    total=sum(balanced.values())
    internal=sum(c for d,c in balanced.items() if d[0]==p)
    return {'balanced_blue_ports':total,'internal_blue_ports':internal,'boundary_blue_ports':total-internal},peak

def coefficient_polynomial(n):
    """Direct ordinary degree product control, not the row-forgetting DP."""
    p=n//2; poly={(1,)*(n-1)+(0,):1};peak=1
    for i in range(n):
        for j in range(i+2,n):
            nxt=Counter()
            for d,c in poly.items():
                for winner in (i,j):
                    if d[winner]>=p:continue
                    e=list(d);e[winner]+=1;nxt[tuple(e)]+=c
            poly=nxt;peak=max(peak,len(poly))
    counts=Counter()
    for d,c in poly.items():
        if all(a in (p-1,p) for a in d):
            counts[('H' if d[0]==p else 'L')+('H' if d[-1]==p else 'L')]+=c
    return {k:counts[k] for k in ('LL','LH','HL','HH')},peak

def summary(n):
    endpoint,states=balanced_endpoint_dp(n)
    blue,peak=blue_half_polynomial(n//2)
    total=sum(endpoint.values());internal=endpoint['HL']; boundary=total-internal
    r={'n':n,'full_cube_words':1<<comb(n-1,2),'balanced_ports':total,'endpoint_ports':endpoint,
       'internal_ports':internal,'boundary_ports':boundary,**blue,
       'internal_black_ports':internal-blue['internal_blue_ports'],
       'boundary_black_ports':boundary-blue['boundary_blue_ports'],
       'row_dp_states':states,'blue_polynomial_peak_states':peak}
    if n>=4:
        need(internal%2==blue['internal_blue_ports']%2==0,'free C pairing parity')
        r.update(internal_lines=internal//2,internal_blue_lines=blue['internal_blue_ports']//2,
                 internal_black_lines=(internal-blue['internal_blue_ports'])//2,
                 incident_lines=internal//2+boundary)
    return r

if __name__=='__main__':
    import argparse,json
    a=argparse.ArgumentParser();a.add_argument('n',type=int);args=a.parse_args()
    print(json.dumps(summary(args.n),indent=2,sort_keys=True))


def eulerian_anticycle_dp(n):
    """Independent D_n count on K_(n+1) minus the undirected Hamilton cycle."""
    N=n+1; target=n//2-1
    @lru_cache(None)
    def visit(wins):
        r=len(wins);i=N-r
        if r==1:return int(wins[0]==target)
        # Relative future indices j=1,...,r-1; omit cycle neighbors.
        neighbors=tuple(j for j in range(1,r) if j!=1 and not (i==0 and j==r-1))
        k=target-wins[0]
        if not 0<=k<=len(neighbors):return 0
        answer=0;neigh=set(neighbors)
        for selected_tuple in combinations(neighbors,k):
            selected=set(selected_tuple)
            new=tuple(wins[j]+int(j in neigh and j not in selected) for j in range(1,r))
            if any(x>target for x in new):continue
            answer+=visit(new)
        return answer
    count=visit((0,)*N)
    return count,visit.cache_info().currsize
