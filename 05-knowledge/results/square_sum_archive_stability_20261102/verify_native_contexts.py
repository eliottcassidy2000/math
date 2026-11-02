#!/usr/bin/env python3
"""Exact induced square-sum contexts; source labels and selected edges retained."""
from itertools import combinations,permutations
from math import isqrt
import json

def need(b,s):
    if not b: raise RuntimeError(s)
def edge(a,b):return tuple(sorted((a,b)))
def issquare(n):return n>0 and isqrt(n)**2==n
def induced(xs):return {edge(a,b) for a,b in combinations(xs,2) if issquare(a+b)}
def components(xs,es):
    adj={v:set() for v in xs}
    for a,b in es:adj[a].add(b);adj[b].add(a)
    seen=set();out=[]
    for v in xs:
        if v in seen:continue
        todo=[v];c=[];seen.add(v)
        while todo:
            x=todo.pop();c.append(x)
            for y in adj[x]:
                if y not in seen:seen.add(y);todo.append(y)
        out.append(sorted(c))
    return sorted(out)
def compile_context(source,connectors):
    xs=list(source);base=induced(xs);added=set();paths=[]
    for A,B in connectors:
        H=max(xs);r=2*H+5
        if (r-A-B-1)%2:r+=1
        s=(A+B+1-r*r)//2;v=r*r-B;w=s*s-A
        need(H<v<w,'range and injectivity')
        path=[A,w,v,B]
        added.update(edge(a,b) for a,b in zip(path,path[1:]))
        xs.extend([w,v]);paths.append(path)
        need(induced(xs)==base|added,'unintended induced edge')
    return xs,added,paths
def main():
    rows=[];cases=0
    for scale in [1,4,9,25]:
        S=[scale*x for x in [1,3,35,46]]
        F={edge(S[0],S[1]),edge(S[2],S[3])}
        Q={edge(S[0],S[2]),edge(S[1],S[3])}
        need(induced(S)==F|Q,'source four-cycle')
        for accept,reject in [(F,Q),(Q,F)]:
            for order in permutations(sorted(reject)):
                for flags in range(4):
                    ctx=[pair if not flags>>i&1 else pair[::-1] for i,pair in enumerate(order)]
                    xs,added,paths=compile_context(S,ctx)
                    ca=components(xs,accept|added);cr=components(xs,reject|added)
                    for es in [accept|added,reject|added]:
                        need(all(sum(x in e for e in es)==2 for x in xs),'degree ledger')
                    need(len(ca)==1 and sorted(map(len,cr))==[4,4],'separator')
                    cases+=1
                    if scale==1 and len(rows)<2:
                        rows.append({'source':S,'accepted_edges':sorted(accept),'rejected_edges':sorted(reject),'connector_paths':paths,'induced_edges':sorted(induced(xs)),'accepted_components':ca,'rejected_components':cr})
    # Eight shared ports: the accepted state must join all four fragments,
    # including those unaffected by the pairing difference.
    S=[1,3,35,46,100,300,3500,4600]
    F={edge(1,3),edge(35,46),edge(100,300),edge(3500,4600)}
    Q={edge(1,35),edge(3,46),edge(100,300),edge(3500,4600)}
    need(F|Q<=induced(S),'eight-port source edges')
    N=[(3,100),(300,3500),(4600,46),(35,1)]
    xs,added,paths=compile_context(S,N)
    ca=components(xs,F|added);cr=components(xs,Q|added)
    need(len(ca)==1 and len(cr)>1,'eight-port separator')
    need(all(sum(x in e for e in F|added)==2 for x in xs),'eight-port accepted degrees')
    need(all(sum(x in e for e in Q|added)==2 for x in xs),'eight-port rejected degrees')
    cases+=1
    extended={'ports':len(S),'added_vertices':len(xs)-len(S),'max_label_bits':max(xs).bit_length(),'accepted_component_sizes':list(map(len,ca)),'rejected_component_sizes':list(map(len,cr))}
    print(json.dumps({'status':'passed','cases':cases,'source_sets':5,'witnesses':rows,'eight_port_control':extended,'scope':'finite exact implementation controls; all-size proof is NATIVE_CONTEXTS.md'},sort_keys=True,indent=2))
if __name__=='__main__':main()
