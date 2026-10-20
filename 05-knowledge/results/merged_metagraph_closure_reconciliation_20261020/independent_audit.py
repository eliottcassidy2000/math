#!/usr/bin/env python3
"""Independent combinatorial audit; stdout is JSON; no checkout writes.
Uses adjacency bitsets and Hamiltonian-path DFS, not audit_incoming.py.
"""
from itertools import combinations, permutations
from fractions import Fraction
from collections import Counter, defaultdict
import json


def require(p, message):
    if not p:
        raise RuntimeError(message)


def tournament(n, mask):
    rows = [0] * n
    for bit, (u, v) in enumerate(combinations(range(n), 2)):
        if mask & (1 << bit):
            rows[u] |= 1 << v
        else:
            rows[v] |= 1 << u
    return tuple(rows)


def arc(A, u, v):
    return bool(A[u] & (1 << v))


def path_counts(A):
    n = len(A)
    full = (1 << n) - 1
    total = closing = 0
    def dfs(first, last, used):
        nonlocal total, closing
        if used == full:
            total += 1
            closing += arc(A, last, first)
            return
        available = A[last] & ~used
        while available:
            bit = available & -available
            available -= bit
            dfs(first, bit.bit_length() - 1, used | bit)
    for start in range(n):
        dfs(start, start, 1 << start)
    return total, closing


def canonical(A, perms):
    pairs = tuple(combinations(range(len(A)), 2))
    codes = [sum(int(arc(A, p[u], p[v])) << bit for bit, (u, v) in enumerate(pairs)) for p in perms]
    original = sum(int(arc(A, u, v)) << bit for bit, (u, v) in enumerate(pairs))
    return min(codes), codes.count(original)


def strong(A):
    all_vertices = (1 << len(A)) - 1
    for s in range(len(A)):
        visited = 1 << s
        while True:
            expanded = visited
            for u in range(len(A)):
                if visited & (1 << u):
                    expanded |= A[u]
            if expanded == visited:
                break
            visited = expanded
        if visited != all_vertices:
            return False
    return True


def closure_census():
    out = []
    for n in range(3, 6):
        perms = list(permutations(range(n)))
        universe = 1 << (n * (n - 1) // 2)
        classes = {}
        closed_labels = 0
        for mask in range(universe):
            A = tournament(n, mask)
            H, C = path_counts(A)
            require(H % 2 == 1, 'Hamiltonian path parity control failed')
            require(C % n == 0, 'cycle rotation count failed')
            if H == C:
                code, aut = canonical(A, perms)
                closed_labels += 1
                classes[code] = {'canonical_code':code, 'H':H, 'HC':C//n, 'Aut':aut, 'scores':sorted(row.bit_count() for row in A)}
        out.append({'n':n, 'labelled_universe':universe, 'fully_closed_labelled':closed_labels, 'classes':list(classes.values())})
    require([len(x['classes']) for x in out] == [1, 0, 1], 'closed classification failed')
    # Explicit Q5 model: singleton a -> C3 -> singleton b -> a.
    Q = [0]*5
    for u, v in [(0,1),(0,2),(0,3),(1,2),(2,3),(3,1),(1,4),(2,4),(3,4),(4,0)]:
        Q[u] |= 1 << v
    require(canonical(Q, list(permutations(range(5))))[0] == out[2]['classes'][0]['canonical_code'], 'Q5 model misidentified')
    # Necessary hostile: regular cyclic C5 is not fully closed.
    C5 = tuple(sum(1 << ((u+d)%5) for d in (1,2)) for u in range(5))
    H, C = path_counts(C5)
    require((H,C)==(15,10), 'C5 hostile failed')
    return out, {'C5_H':H, 'C5_closing_words':C, 'C5_open_words':H-C}


def asymptotic_check():
    T = [1] + [1 << ((n-1)*(n-2)//2) for n in range(1,101)]
    a = [0]*101
    NS = [0]*101
    for n in range(1,101):
        a[n] = T[n] - sum(a[k] * T[n-k] for k in range(1,n))
        NS[n] = T[n] - a[n]
        require(0 <= a[n] <= T[n], 'negative strong count')
        if n>=4:
            scale = 1 << (n-3)
            # Clear denominators; no floating-point tests.
            require(NS[n]*scale >= T[n-1]*(2*scale-(n-2)), 'lower bound failed')
            require(NS[n]*scale <= T[n-1]*(2*scale+(n-3)), 'upper bound failed')
            require(NS[n-1]*scale <= (n-2)*T[n-1], 'previous-order bound failed')
            for k in range(1,n):
                require(Fraction(T[k]*T[n-k],T[n-1]) == Fraction(1,1<<((k-1)*(n-k-1))), 'exponent identity failed')
    return {'bounds_verified_orders':[4,100], 'integer_arithmetic':True, 'controls':[{'n':n,'NS':str(NS[n]),'T_previous':str(T[n-1]),'ratio':float(Fraction(NS[n],T[n-1]))} for n in (4,8,16,40,100)]}


def criterion_hostile():
    n=9
    S={1,3,4,7}
    P=(3,1,2,6,4,5,0,7,8)
    A=tuple(sum(1<<((u+d)%n) for d in S) for u in range(n))
    require(sorted(P)==list(range(n)), 'hostile is not permutation')
    require(all(arc(A,u,v) for u,v in zip(P,P[1:])), 'hostile is not directed path')
    reversal=[None]*n
    for i,u in enumerate(P):
        reversal[u]=P[-1-i]
    checks=[]
    for u in range(n):
        for v in range(n):
            if arc(A,u,v):
                require(arc(A,reversal[v],reversal[u]), 'hostile reversal not anti')
                checks.append([u,v,reversal[v],reversal[u]])
    sums=[(P[i]+P[-1-i])%n for i in range(n)]
    require(len(set(sums))>1, 'hostile sums constant')
    sigma=[(u+3)%n if u%3==0 else u for u in range(n)]
    require(tuple(sigma)==P, 'module automorphism does not produce path')
    require(all(arc(A,sigma[u],sigma[v]) for u in range(n) for v in range(n) if arc(A,u,v)), 'module shift not automorphism')
    return {'n':n,'S':sorted(S),'path':P,'reversal_by_vertex':reversal,'pair_sums':sums,'directed_arcs_checked':len(checks),'arc_images':checks}


def incoming_path_retest():
    result=[]
    for n in range(7,24,2):
        counts=Counter()
        half=(n-1)//2
        for bits in range(1<<(half-1)):
            S={1}|{d if bits & (1<<(d-2)) else n-d for d in range(2,half+1)}
            candidates=[]
            for a in sorted(S-{1}):
                for b in sorted(S-{1}):
                    if a<b and a+b<=n+1:
                        candidates.append(list(range(a-1))+list(range(a+b-2,n))+list(range(a-1,a+b-2)))
            if not (S & set(range(2,half+1))):
                candidates.append([0,1]+[x for j in range(1,half) for x in (n-2*j,n-2*j+1)]+[2])
            require(candidates, 'incoming selection has no candidate')
            counts['connection_sets']+=1
            for P in candidates:
                require(sorted(P)==list(range(n)), 'incoming path not permutation')
                require(all((v-u)%n in S for u,v in zip(P,P[1:])), 'incoming path not directed')
                rev=[0]*n
                for i,u in enumerate(P):
                    rev[u]=P[-1-i]
                # Find an ACTUAL violated directed arc, not a nonconstant sum.
                bad=next(((u,v) for u in range(n) for v in range(n) if (v-u)%n in S and (rev[u]-rev[v])%n not in S),None)
                require(bad is not None, 'incoming candidate actually anti')
                counts['paths_with_violated_arc']+=1
        result.append({'n':n,**counts})
    require(sum(x['paths_with_violated_arc'] for x in result)==41992, 'unexpected incoming path count')
    return result


def incidence_checks():
    result=[]
    for n in range(3,7):
        free=list((u,v) for u,v in combinations(range(n),2) if v>u+1)
        full=(1<<len(free))-1
        perms=list(permutations(range(n)))
        infos=[]
        for mask in range(full+1):
            A=[0]*n
            for u in range(n-1):
                A[u]|=1<<(u+1)
            for bit,(u,v) in enumerate(free):
                if mask & (1<<bit): A[u]|=1<<v
                else: A[v]|=1<<u
            codes=[sum(int(arc(A,p[u],p[v]))<<bit for bit,(u,v) in enumerate(combinations(range(n),2))) for p in perms]
            ordinary=min(codes)
            converse=((1<<(n*(n-1)//2))-1)-max(codes)
            node=min(ordinary,converse)
            blue=all(arc(A,u,v)==arc(A,n-1-v,n-1-u) for u,v in combinations(range(n),2))
            infos.append({'node':node,'closed':arc(A,n-1,0),'blue':blue,'strong':strong(A)})
        nodes=sorted({x['node'] for x in infos})
        strong_nodes={x['node'] for x in infos if x['strong']}
        edges=[]
        for mask in range(full+1):
            other=full^mask
            if mask>other: continue
            x,y=infos[mask],infos[other]
            require(x['closed'] != y['closed'], 'literal complement does not flip closure')
            require(x['blue']==y['blue'], 'literal complement does not preserve colour')
            origin,target=(y,x) if x['closed'] else (x,y)
            require(target['strong'], 'closed target not strong')
            edges.append((origin['node'],target['node'],origin['blue']))
        shores=[set(),set(nodes),strong_nodes,set(nodes)-strong_nodes]+[{u} for u in nodes]
        if n<=5:
            shores=[{nodes[i] for i in range(len(nodes)) if mask&(1<<i)} for mask in range(1<<len(nodes))]
        testcount=0
        for shore in shores:
            for colour in (None,True,False):
                hs=[x for x in infos if x['node'] in shore and (colour is None or x['blue']==colour)]
                C=sum(x['closed'] for x in hs);O=len(hs)-C
                eligible=[e for e in edges if colour is None or e[2]==colour]
                incoming=sum(u not in shore and v in shore for u,v,_ in eligible)
                outgoing=sum(u in shore and v not in shore for u,v,_ in eligible)
                internal=sum(u in shore and v in shore for u,v,_ in eligible)
                d=incoming+outgoing
                require(incoming-outgoing==C-O,'oriented cut difference failed')
                require(d>=abs(C-O),'oriented cut lower bound failed')
                require(d==C+O-2*internal,'oriented incidence identity failed')
                testcount+=1
        if n>=6:
            for u in strong_nodes:
                require(any(v==u for _,v,_ in edges),'strong node has no incoming edge')
                require(any(v==u for v,_,_ in edges),'strong node has no outgoing edge')
        result.append({'n':n,'tilings':full+1,'merged_nodes':len(nodes),'strong_nodes':len(strong_nodes),'lines':len(edges),'shore_colour_tests':testcount,'all_shores':n<=5})
    return result


if __name__=='__main__':
    small, controls=closure_census()
    output={'review_sha256':'637ef8c930f98d37155a48f140f113d4d1831dcf3e6f92d30a004e6005422d7d','small_closure_census':small,'closure_controls':controls,'asymptotic':asymptotic_check(),'criterion_hostile':criterion_hostile(),'incoming_explicit_paths':incoming_path_retest(),'incidence_controls':incidence_checks(),'status':'PASS'}
    print(json.dumps(output,indent=2))
