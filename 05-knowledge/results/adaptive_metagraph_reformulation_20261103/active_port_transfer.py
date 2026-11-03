#!/usr/bin/env python3
"""Exact temporal active-port path-cover transfer; no repository imports.

The degree control is finite. Integer cover weights remain unbounded.
Only small history controls are compared with actual subset DP / literal covers.
"""
from collections import Counter, deque
from functools import lru_cache
from itertools import permutations
import json


def need(test, why):
    if not test: raise RuntimeError(why)


def birth(T,v):
    n=len(T)
    return tuple(row|((1<<n) if i==v else 0) for i,row in enumerate(T))+(((1<<n)-1)^(1<<v),)


def active(T,d):
    return tuple(v for v,row in enumerate(T) if len(T)-1-row.bit_count()<=d)


def degree_control(T,d):
    return tuple(len(T)-1-T[v].bit_count() for v in active(T,d))


def control_step(degrees,root,d):
    need(0<=root<len(degrees),'legal active root rank')
    raised=tuple(a+(i!=root) for i,a in enumerate(degrees))
    keep=tuple(i for i,a in enumerate(raised) if a<=d)
    return tuple(raised[i] for i in keep)+(1,),keep


def controls(d):
    T=birth((2,4,1),0);seed=degree_control(T,d)
    seen={seed};todo=deque([seed]);arcs=[]
    while todo:
        state=todo.popleft()
        for root in range(len(state)):
            dest,_=control_step(state,root,d)
            need(len(dest)<=2*d,'active bound')
            arcs.append((state,root,dest))
            if dest not in seen:seen.add(dest);todo.append(dest)
    return seed,tuple(sorted(seen)),tuple(sorted(arcs))


@lru_cache(None)
def endpoint_step(state,degrees,root,d):
    """State is an unordered tuple of paired (start,end) ports; -1 forgotten."""
    dest_degrees,keep=control_step(degrees,root,d)
    r=len(degrees);renumber={old:new for new,old in enumerate(keep)}
    renumber[r]=len(keep)
    incoming=[None]+[i for i,(_,e) in enumerate(state) if e==root]
    outgoing=[None]+[i for i,(s,_) in enumerate(state) if s!=root]
    out=Counter()
    for i in incoming:
        for j in outgoing:
            if i is not None and i==j:continue
            parts=[p for k,p in enumerate(state) if k!=i and k!=j]
            parts.append((r if i is None else state[i][0],r if j is None else state[j][1]))
            parts=[(renumber.get(s,-1),renumber.get(e,-1)) for s,e in parts]
            if sum(e==-1 for _,e in parts)>1:continue
            out[tuple(sorted(parts))]+=1
    return dest_degrees,tuple(sorted(out.items()))


def advance(vec,degrees,root,d):
    out=Counter();dest=None
    for state,weight in vec.items():
        dest,terms=endpoint_step(state,degrees,root,d)
        for target,multiplicity in terms:out[target]+=weight*multiplicity
    need(dest is not None,'nonempty cover vector')
    return dest,out


def literal_cover_vector(T,ports):
    """Independent full permutations and cuts, deduplicated actual path covers."""
    n=len(T);covers=set();port={v:i for i,v in enumerate(ports)}
    for word in permutations(range(n)):
        for cuts in range(1<<(n-1)):
            paths=[];start=0;valid=True
            for i in range(n-1):
                if cuts>>i&1:
                    paths.append(word[start:i+1]);start=i+1
                elif not (T[word[i]]>>word[i+1]&1):valid=False;break
            if valid:
                paths.append(word[start:]);covers.add(tuple(sorted(paths)))
    vec=Counter()
    for cover in covers:
        state=tuple(sorted((port.get(p[0],-1),port.get(p[-1],-1)) for p in cover))
        if sum(e==-1 for _,e in state)<=1:vec[state]+=1
    return vec


def H(T):
    n=len(T);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,value in enumerate(dp[mask]):
            choices=T[v]&~mask
            while value and choices:
                bit=choices&-choices;choices-=bit
                dp[mask|bit][bit.bit_length()-1]+=value
    return sum(dp[-1])


def output(vec):return sum(weight for state,weight in vec.items() if len(state)==1)


def endpoint_bounds(k):
    # g: no forgotten end permitted; f: at most one permitted.
    g=[1];f=[2]
    for r in range(1,k+1):
        g.append(3*g[-1]+(2*(r-1)*g[r-2] if r>1 else 0))
        f.append(3*f[-1]+g[r-1]+(2*(r-1)*f[r-2] if r>1 else 0))
    return f


def carrier(d,seed,vec):
    seen={(seed,state) for state in vec};todo=deque(seen)
    while todo:
        control,state=todo.popleft()
        for root in range(len(control)):
            dest,terms=endpoint_step(state,control,root,d)
            for target,_ in terms:
                pair=(dest,target)
                if pair not in seen:seen.add(pair);todo.append(pair)
    per=Counter(control for control,_ in seen)
    return {'reachable_control_endpoint_pairs':len(seen),
            'endpoint_pairs_by_control':[{ 'degrees':c,'count':per[c]} for c in sorted(per)]}


def verify_small(d,maxn):
    T=birth((2,4,1),0);degrees=degree_control(T,d)
    vec=literal_cover_vector(T,active(T,d));todo=[(T,degrees,vec)]
    checks=0;cover_checks=1
    while todo:
        T,degrees,vec=todo.pop();n=len(T)
        need(output(vec)==H(T),'actual subset DP')
        need(degrees==degree_control(T,d),'actual active-degree control')
        checks+=1
        if n>=maxn:continue
        ports=active(T,d)
        for root,v in enumerate(ports):
            U=birth(T,v);new_degrees,new_vec=advance(vec,degrees,root,d)
            if len(U)<=6:
                need(new_vec==literal_cover_vector(U,active(U,d)),'full literal endpoint distribution')
                cover_checks+=1
            todo.append((U,new_degrees,new_vec))
    return {'actual_named_graphs':checks,'literal_endpoint_distributions':cover_checks}


def main():
    receipt={'scope':'all-order finite carrier theorem; only small named controls, no new target-support census',
             'endpoint_signature_upper_bounds_k0_to6':endpoint_bounds(6),'controls':{}}
    for d in (1,2,3):
        seed,states,arcs=controls(d)
        entry={'seed':seed,'number':len(states),'states':states,
               'arcs':[{'source':s,'root_rank':r,'target':t} for s,r,t in arcs]}
        T=birth((2,4,1),0);vec=literal_cover_vector(T,active(T,d))
        entry.update(carrier(d,seed,vec))
        entry['small_actual_controls']=verify_small(d,7 if d<=2 else 6)
        receipt['controls'][str(d)]=entry
    # Exact cap-three realization and forbidden final cap-two root in the parent witness.
    word=(0,0,3,3,0);T=(2,4,1);used=[]
    for v in word:
        used.append(len(T)-1-T[v].bit_count());T=birth(T,v)
    need(used==[1,1,2,2,3] and H(T)==127,'actual degree-three 127 witness')
    receipt['H127']={'word':word,'root_indegrees':used,'rows':T,'H':127}
    # Unbounded static feedback cost inside cap one: literal disjoint triangles.
    T=birth((2,4,1),0)
    for n in range(4,30):T=birth(T,n-1)
    triples=[(3*i,3*i+1,3*i+2) for i in range(10)]
    need(all(sum((T[u]>>v)&1 for v in triple if v!=u)==1 for triple in triples for u in triple),'disjoint triangle family')
    receipt['static_feedback_control']={'order':len(T),'disjoint_C3':len(triples),'degree_cap':1}
    return receipt


if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
