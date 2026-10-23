"""Literal cut-budget certificates and hostile controls; no asserts are used."""
from pathlib import Path
from collections import Counter
from itertools import combinations,permutations
from math import isqrt
import argparse,json

HERE=Path(__file__).resolve().parent
checks=0
def need(condition,why):
    global checks
    checks+=1
    if not condition:raise RuntimeError(why)

def edges(n):
    return [(a,b) for a in range(1,n+1) for b in range(a+1,n+1) if isqrt(a+b)**2==a+b]

def components(vertices,edge_list):
    remaining=set(vertices);result=[]
    adjacent={v:set() for v in remaining}
    for a,b in edge_list:
        if a in remaining and b in remaining:adjacent[a].add(b);adjacent[b].add(a)
    while remaining:
        stack=[min(remaining)];found=set(stack);remaining-=found
        while stack:
            v=stack.pop()
            new=adjacent[v]&remaining;remaining-=new;found|=new;stack.extend(new)
        result.append(sorted(found))
    return sorted(result)

def cycle_edge_witnesses(cycle,E):
    """Every edge of a finite simple Hamiltonian graph is used and omitted.

    A chord uv is included by deleting u--succ(u) and v--succ(v) from
    the oriented cycle and then adding uv. Keep the full spanning path.
    """
    n=len(cycle);vertices=set(cycle)
    need(n>=3 and len(vertices)==n,'simple cycle has at least three vertices')
    CE={tuple(sorted(e)) for e in zip(cycle,cycle[1:]+cycle[:1])}
    need(len(CE)==n and CE<=E,'cycle edges present')
    successor={cycle[i]:cycle[(i+1)%n] for i in range(n)}
    def path_from_edges(F):
        adjacent={v:[] for v in cycle}
        for a,b in F:adjacent[a].append(b);adjacent[b].append(a)
        ends=sorted(v for v in cycle if len(adjacent[v])==1)
        need(len(F)==n-1 and len(ends)==2,'path degrees after cycle surgery')
        need(all(len(adjacent[v]) in (1,2) for v in cycle),'no branching after cycle surgery')
        p=[ends[0]]
        while len(p)<n:
            choices=[v for v in adjacent[p[-1]] if len(p)==1 or v!=p[-2]]
            need(len(choices)==1 and choices[0] not in p,'connected acyclic surgery')
            p.append(choices[0])
        need(set(p)==vertices and p[-1]==ends[1],'spanning endpoint witness')
        return p
    included={};omitted={}
    for e in sorted(E):
        cut=e if e in CE else min(CE)
        no=path_from_edges(CE-{cut})
        if e in CE:
            deleted={next(f for f in sorted(CE) if f!=e)}
            yes=path_from_edges(CE-deleted)
        else:
            u,v=e;deleted={tuple(sorted((u,successor[u]))),tuple(sorted((v,successor[v])))}
            need(len(deleted)==2,'two distinct cuts for a chord')
            yes=path_from_edges((CE-deleted)|{e})
        need(e in {tuple(sorted(x)) for x in zip(yes,yes[1:])},('edge inclusion',e))
        need(e not in {tuple(sorted(x)) for x in zip(no,no[1:])},('edge omission',e))
        included[str(e)]={'path':yes,'deleted_cycle_edges':[list(x) for x in sorted(deleted)]}
        omitted[str(e)]=no
    return included,omitted

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    bank={row['n']:row['paths'] for row in json.loads((HERE/'paths.json').read_text())['orders']}
    # Explicit certificates, not an assertion that all forbidden edges have this form.
    cases=[(15,[1,3]),(16,[1,3]),(17,[1,3]),(23,[1,4,6,7,13,14,15]),
           (23,[3,4,6,7,14,15]),(23,[3,4,7,13,14]),(23,[3,4,5,7,14]),
           (26,[1,3,5,7,10,12,14,23]),
           (29,[7,9,20]),(30,[7,9,20])]
    out={'status':'PROVED cut-budget lemma; FINITE-EXACT explicit square-sum instances and hostile controls','tight_cuts':[],'obstructions':[]}
    for n,S in cases:
        outside=set(range(1,n+1))-set(S);cc=components(outside,edges(n));budget=len(S)+1-len(cc)
        need(bank[n] and budget==0,('nonvacuous tight cut',n,S))
        inside=[list(e) for e in edges(n) if set(e)<=set(S)]
        for p in bank[n]:
            path_edges=list(zip(p,p[1:]));ends_in=sum(v in S for v in (p[0],p[-1]));internal=sum(a in S and b in S for a,b in path_edges)
            actual=len(components(outside,path_edges))
            need(actual==len(S)+1-ends_in-internal,('exact path deletion identity',n,S,p))
            need(len(cc)<=actual and ends_in+internal<=budget,('cut budget',n,S,p))
        out['tight_cuts'].append({'n':n,'cut':S,'outside_components':cc,'budget':budget,'forbidden_edges_inside_cut':inside,'path_count':len(bank[n])})
    for n,S in [(18,[7]),(19,[7]),(20,[5,7]),(21,[5,7]),(22,[5,7,14])]:
        cc=components(set(range(1,n+1))-set(S),edges(n))
        need(len(cc)>len(S)+1 and not bank[n],('strict cut obstruction',n,S))
        out['obstructions'].append({'n':n,'cut':S,'outside_components':cc,'negative_budget':len(S)+1-len(cc)})
    # Hostile: degree two does not force both edges when that vertex is an endpoint.
    witness=next(p for p in bank[23] if 22 in (p[0],p[-1]) and not any({a,b}=={14,22} for a,b in zip(p,p[1:])))
    need(sum(22 in e for e in edges(23))==2,'degree-two hostile vertex')
    out['degree_two_hostile']={'n':23,'vertex':22,'missing_incident_edge':[14,22],'path':witness}
    # Correction: a Hamiltonian cycle implies BOTH inclusion and omission for
    # every graph edge, including chords; the two-cut construction is essential.
    cycle=[1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15]
    E=set(edges(32));cycle_edges={tuple(sorted(e)) for e in zip(cycle,cycle[1:]+cycle[:1])}
    need(set(cycle)==set(range(1,33)) and len(cycle_edges)==32 and cycle_edges<=E,'literal Hamiltonian cycle')
    included,omitted=cycle_edge_witnesses(cycle,E)
    out['G32_cycle']=cycle
    out['G32_edge_inclusion_witnesses']=included
    out['G32_no_forced_edge_witnesses']=omitted
    # Exhaust all graphs containing the fixed labelled cycle on3..6 vertices.
    controls=[]
    for n in range(3,7):
        C=list(range(1,n+1));CE={tuple(sorted(e)) for e in zip(C,C[1:]+C[:1])}
        chords=[e for e in combinations(C,2) if e not in CE];edge_witnesses=0
        for mask in range(1<<len(chords)):
            F=CE|{e for i,e in enumerate(chords) if mask>>i&1}
            yes,no=cycle_edge_witnesses(C,F);need(len(yes)==len(no)==len(F),'complete chord control');edge_witnesses+=len(F)
        controls.append({'n':n,'graphs':1<<len(chords),'edges_with_both_witnesses':edge_witnesses})
    out['general_cycle_controls']=controls
    positions=[]
    for n in range(3,65):
        C=list(range(1,n+1));complete=set(combinations(C,2))
        yes,no=cycle_edge_witnesses(C,complete)
        need(len(yes)==len(no)==n*(n-1)//2,'all cycle/chord positions')
        positions.append({'n':n,'positions_with_both_witnesses':len(complete)})
    out['cycle_chord_position_controls']=positions
    # Hostile to weakening the hypothesis from a spanning cycle to a path:
    # a triangle with leaves attached to two vertices has a unique spanning path.
    H={(1,2),(1,3),(2,3),(1,4),(2,5)}
    HP=[p for p in permutations(range(1,6)) if p[0]<p[-1] and all(tuple(sorted(e)) in H for e in zip(p,p[1:]))]
    need(HP==[(4,1,3,2,5)],'traceable graph can have a forbidden edge')
    out['traceable_is_not_enough']={'edges':[list(e) for e in sorted(H)],'unique_path':list(HP[0]),'never_edge':[1,2]}
    C4={(1,2),(2,3),(3,4),(1,4)}
    pinned=[p for p in permutations(range(1,5)) if p[0]==1 and p[-1]==2 and all(tuple(sorted(e)) in C4 for e in zip(p,p[1:]))]
    need(pinned==[(1,4,3,2)],'prescribed endpoints do not inherit optionality')
    out['prescribed_endpoint_hostile']={'cycle':[1,2,3,4],'endpoints':[1,2],'unique_path':list(pinned[0]),'never_edge':[1,2]}
    out['scope']='Tight cuts shown are sufficient, not complete. Every edge of every Hamiltonian finite simple graph is optional among its spanning paths. The cited square-sum Hamiltonian-cycle theorem therefore makes every G_n edge optional for all n>=32; exact usage frequencies are enumerated only through32 here.'
    out['checks']=checks
    if args.write:(HERE/'cuts.json').write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads((HERE/'cuts.json').read_text())==out,'frozen cut receipt')
    print('cut checks',checks)

if __name__=='__main__':main()
