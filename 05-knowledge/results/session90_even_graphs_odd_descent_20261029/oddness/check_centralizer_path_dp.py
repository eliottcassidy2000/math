#!/usr/bin/env python3
"""Equitable centralizer-orbit compression for the folded Hamilton-path observer.
Imports the companion exact row/codec routines only; stdlib, active tests under -O.
"""
from itertools import permutations,product
from collections import Counter,defaultdict
import json
from check_fixed_orbit_descent import lift,arc,is_aut,comp,folded_paths,literal_paths

CHECKS=0
def need(c,msg='orbit-DP check failed'):
    global CHECKS
    CHECKS+=1
    if not c:raise AssertionError(msg)

def compressed(rows,tau,C):
    root=next(i for i in range(len(rows)) if tau[i]==i)
    pairs=[(i,tau[i]) for i in range(len(rows)) if i<tau[i]]; n=len(pairs)
    if n==0:return {'count':1,'states':0,'orbits':0,'labelled_arcs':0,'quotient_arcs':0,'orbit_sizes':{}}
    vertices=[v for pair in pairs for v in pair];pairindex={v:i for i,pair in enumerate(pairs) for v in pair}
    C=set(C);need(tuple(range(len(rows))) in C)
    for g in C:
        need(g[root]==root and comp(g,tau)==comp(tau,g) and is_aut(rows,g))
        for h in C:need(comp(g,h) in C)
    gp={g:tuple(pairindex[g[pair[0]]] for pair in pairs) for g in C}
    def act(g,state):
        mask,v=state
        return sum(1<<gp[g][i] for i in range(n) if mask>>i&1),g[v]
    states={(mask,v) for mask in range(1,1<<n) for v in vertices if mask>>pairindex[v]&1}
    owner={}; members={}
    while len(owner)<len(states):
        state=min(states-owner.keys());orbit={act(g,state) for g in C};rep=min(orbit)
        need(orbit<=states);need(not orbit&owner.keys());members[rep]=orbit
        for x in orbit:owner[x]=rep
    def successors(state):
        mask,v=state
        return [(mask|(1<<pairindex[w]),w) for w in vertices if not mask>>pairindex[w]&1 and arc(rows,v,w)]
    out={rep:Counter(owner[y] for y in successors(rep)) for rep in members}
    # Verify equitable counts directly for every state, not only representatives.
    for rep,orbit in members.items():
        for state in orbit:need(Counter(owner[y] for y in successors(state))==out[rep])
    P=defaultdict(int);badP=defaultdict(int)
    for rep,orbit in members.items():
        if rep[0].bit_count()==1:P[rep]=len(orbit);badP[rep]=1
    for rep in sorted(members,key=lambda s:(s[0].bit_count(),s)):
        for target,multiplicity in out[rep].items():
            P[target]+=P[rep]*multiplicity;badP[target]+=badP[rep]*multiplicity
    full=(1<<n)-1
    total=sum(P[rep] for rep in members if rep[0]==full and arc(rows,rep[1],root))
    # Compare TOTAL count to explicit labelled state recurrence at each state orbit.
    dp=defaultdict(int)
    for v in vertices:dp[(1<<pairindex[v],v)]=1
    for state in sorted(states,key=lambda s:(s[0].bit_count(),s)):
        for target in successors(state):dp[target]+=dp[state]
    for rep,orbit in members.items():
        need(P[rep]==sum(dp[state] for state in orbit))
        need(len({dp[state] for state in orbit})==1)
    badtotal=sum(badP[rep] for rep in members if rep[0]==full and arc(rows,rep[1],root))
    return {'count':total,'wrong_unit_initialization':badtotal,'states':len(states),'orbits':len(members),
            'labelled_arcs':sum(len(successors(s)) for s in states),
            'quotient_arcs':sum(len(x) for x in out.values()),
            'orbit_sizes':dict(sorted(Counter(len(o) for o in members.values()).items()))}

def normalized_centralizer(rows,n):
    # All centralizing automorphisms preserve p_i -> q_i, so only pair permutations remain.
    result=[]
    for p in permutations(range(n)):
        g=(0,)+tuple(1+p[i] for i in range(n))+tuple(1+n+p[i] for i in range(n))
        if is_aut(rows,g):result.append(g)
    return result

def wreath9():
    rows=tuple(sum(1<<(3*k+l) for k in range(3) for l in range(3)
                   if (k!=i and (k-i)%3==1) or (k==i and (l-j)%3==1))
               for i in range(3) for j in range(3))
    tau=tuple(3*((-i)%3)+(-j)%3 for i in range(3) for j in range(3))
    G=[tuple(3*((i+a)%3)+(j+b[i])%3 for i in range(3) for j in range(3))
       for a in range(3) for b in product(range(3),repeat=3)]
    C=[g for g in G if comp(g,tau)==comp(tau,g)]
    need(len(G)==81 and len(C)==3)
    need(all(is_aut(rows,g) for g in G))
    # Independently exhaustive permutations certify the full automorphism group.
    exact={g for g in permutations(range(9)) if is_aut(rows,g)}
    need(exact==set(G))
    out=compressed(rows,tau,C);need(out['count']==folded_paths(rows,tau))
    paths=literal_paths(rows);fixed=sum(tuple(tau[v] for v in reversed(p))==p for p in paths)
    need(out['count']==fixed);need(out['count']%len(C)==0)
    fixedvertices=[v for v in range(9) if all(g[v]==v for g in C)]
    need(fixedvertices==[0,1,2]);need(set(out['orbit_sizes'])=={1,3});need(out['wrong_unit_initialization']==19)
    out.update({'Aut':len(G),'centralizer':len(C),'fixed_vertices':fixedvertices,
                'blue_orbits':fixed//len(C),'H':len(paths),
                'naive_start_mass_if_divide_states_by_C':'8/3',
                'actual_start_state_orbits':len({min(g[v] for g in C) for v in range(1,9)})})
    return out

def run():
    totals=Counter();cases=0; nontrivial=0; allrows=[]
    for n in range(4):
        for m in range(1<<(n*n)):
            rows,tau=lift(m,n);C=normalized_centralizer(rows,n)
            out=compressed(rows,tau,C);need(out['count']==folded_paths(rows,tau))
            need(out['count']%len(C)==0)
            for key in ['states','orbits','labelled_arcs','quotient_arcs']:totals[key]+=out[key]
            if len(C)>1:nontrivial+=1;allrows.append({'n':n,'matrix':m,'centralizer':len(C),**out})
            cases+=1
    result={'all_lifts':cases,'nontrivial_centralizer_lifts':nontrivial,
            'aggregate':dict(totals),'nontrivial_rows':allrows,'wreath9':wreath9()}
    result['active_checks']=CHECKS
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':run()
