#!/usr/bin/env python3
"""A native extension hostile to the fixed-G25 complement carrier."""
import json
from itertools import combinations
from math import isqrt
def need(b,s):
    if not b:raise RuntimeError(s)
def bank(n):
    E=[(i,j) for i,j in combinations(range(1,n+1),2) if isqrt(i+j)**2==i+j]
    A={i:[] for i in range(1,n+1)}
    for i,j in E:A[i].append(j);A[j].append(i)
    need([v for v in A if len(A[v])==1]==[18],'forced endpoint')
    paths=[]
    def dfs(p,seen):
        if len(p)==n:paths.append(p);return
        for v in A[p[-1]]:
            if v not in seen:dfs(p+(v,),seen|{v})
    dfs((18,),{18})
    return E,sorted(paths+[p[::-1] for p in paths])
def encode(p,E):
    mask=val=0
    for i,j in zip(p,p[1:]):
        k=E.index(tuple(sorted((i,j))));mask|=1<<k
        if i<j:val|=1<<k
    return mask,val
def main():
    E,P=bank(25);EE,PP=bank(26);new=[e for e in EE if e not in E]
    need(len(P)==20 and len(PP)==24 and new==[(10,26),(23,26)],'source banks')
    old=[encode(p,E) for p in P];ext=[encode(p,E+new) for p in PP];full=(1<<len(E))-1
    def legal(x,patterns):return tuple(i for i,(m,v) in enumerate(patterns) if x&m==v)
    def state(x):
        L=legal(x,old)
        return L,tuple(legal(x^(full^old[i][0]),old) for i in L)
    active=set()
    for m,v in old:
        free=[i for i in range(len(E)) if not m>>i&1]
        for bits in range(1<<len(free)):
            active.add(v|sum(1<<k for j,k in enumerate(free) if bits>>j&1))
    states={x:state(x) for x in active};need(len(active)==4864 and len(set(states.values()))==32,'inherited active carrier')
    need(len({(states[x],tuple(states[x^(full^old[i][0])] for i in states[x][0])) for x in active})==32,'old carrier transition closure')
    groups={};scores={}
    for x in sorted(active):
        for a in range(4):
            h=len(legal(x|(a<<32),ext));scores[x,a]=h
            groups.setdefault((states[x],a),{}).setdefault(h,x)
    split=[(s,a,d) for (s,a),d in groups.items() if len(d)>1]
    # Minimal active-source observer for arbitrary old commands followed by
    # one terminal extension26 query, whose two new edge bits are supplied.
    ordered=sorted(active)
    outputs={x:(states[x],tuple(scores[x,a] for a in range(4))) for x in ordered}
    def compress(signatures):
        names={s:i for i,s in enumerate(sorted(set(signatures.values())))}
        return {x:names[s] for x,s in signatures.items()}
    classes=compress(outputs);rounds=[len(set(classes.values()))]
    for _ in range(len(active)):
        sig={x:(classes[x],tuple(classes[x^(full^old[i][0])] for i in states[x][0])) for x in ordered}
        nxt=compress(sig);rounds.append(len(set(nxt.values())))
        if rounds[-1]==rounds[-2]:break
        classes=nxt
    else:raise RuntimeError('refinement did not stabilize')
    # Start with just current legal names and the terminal outputs instead;
    # this must recover the same equivalence without assuming the old quotient.
    alt=compress({x:(states[x][0],tuple(scores[x,a] for a in range(4))) for x in ordered})
    for _ in range(len(active)):
        nxt=compress({x:(alt[x],tuple(alt[x^(full^old[i][0])] for i in states[x][0])) for x in ordered})
        if len(set(nxt.values()))==len(set(alt.values())):break
        alt=nxt
    need(len(set(zip(classes.values(),alt.values())))==len(set(classes.values()))==len(set(alt.values())),'two observer initializations')
    need(rounds==[36,38,38] and len(split)==4 and len({s for s,a,d in split})==4,'declared native refinement census')
    witness=None
    for x in sorted(active):
        for k in range(32):
            y=x^(1<<k)
            if y<=x or y not in states or states[y]!=states[x]:continue
            for a in range(4):
                if scores[x,a]!=scores[y,a]:witness=(x,y,k,a);break
            if witness:break
        if witness:break
    need(witness is not None,'one-edge native hostile')
    x,y,k,a=witness
    # Independent literal orientation/path comparison for the two witnesses.
    def arcs(z):return {(i,j) if z>>k&1 else (j,i) for k,(i,j) in enumerate(E+new)}
    sources=[]
    for z in [x,y]:
        w=z|(a<<32);A=arcs(w)
        good=[p for p in PP if all(e in A for e in zip(p,p[1:]))]
        need(len(good)==scores[z,a],'literal consequence cross-check')
        sources.append({'guard_bits':z,'H26':len(good),'legal_paths26':good})
    return {'status':'exact native extension hostile found','scope':'All4864 active G25 orientations, all four common orientations of the two new edges; the dead G25 fibre is not enumerated.',
            'edge_order25':E,'new_edges':new,'path_counts_directed':[len(P),len(PP)],
            'paths26_with_new_endpoint':sum(p[0]==26 or p[-1]==26 for p in PP),
            'paths26_with_new_interior':sum(p[0]!=26 and p[-1]!=26 for p in PP),
            'active_carrier_states':32,'fixed_attachment_groups':len(groups),'split_groups':len(split),
            'split_old_states':len({s for s,a,d in split}),
            'terminal_extension_refinement_counts':rounds,'minimal_active_terminal_extension_states':rounds[-1],
            'witness':{'flipped_edge':E[k],'new_bits':a,'new_arcs':[e if a>>i&1 else e[::-1] for i,e in enumerate(new)],
                       'same_old_state':states[x],'old_legal_paths':[P[i] for i in states[x][0]],'sources':sources}}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
