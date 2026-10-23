"""Hybrid cut/degree/cycle proof that no G25 Hamiltonian path uses1--15.

This consumer checks explicit cuts and three deterministic propagation traces;
it does not search for Hamiltonian paths. Completeness comes from the mandatory
leaf endpoint18 and the24 possible other endpoints.
"""
from pathlib import Path
import argparse,json
from cuts import components,edges

HERE=Path(__file__).resolve().parent;checks=0
def need(condition,why):
    global checks
    checks+=1
    if not condition:raise RuntimeError(why)

def propagate(other):
    E=set(edges(25));selected={(1,15)};excluded=set();steps=[]
    while True:
        old=(set(selected),set(excluded))
        for v in range(1,26):
            incident={e for e in E if v in e};present=incident&selected;unknown=incident-selected-excluded
            target=1 if v in (18,other) else 2;needed=target-len(present)
            if needed<0 or needed>len(unknown):
                return {'other_endpoint':other,'steps':steps,'terminal':{'kind':'degree_deficit_or_excess','vertex':v,'target_degree':target,'selected_incident':sorted(present),'available_incident':sorted(unknown)}}
            if needed==0 and unknown:
                excluded|=unknown;steps.append({'rule':'exclude','at_vertex':v,'edges':sorted(unknown),'target_degree':target})
            elif needed==len(unknown) and unknown:
                selected|=unknown;steps.append({'rule':'select','at_vertex':v,'edges':sorted(unknown),'target_degree':target})
        need(not selected&excluded,'inconsistent propagation state')
        if old==(selected,excluded):break
    cycle=[2,14,22,3,13,23]
    cycle_edges={tuple(sorted(e)) for e in zip(cycle,cycle[1:]+cycle[:1])}
    need(cycle_edges<=selected and len(cycle)<25,('forced proper cycle',other))
    return {'other_endpoint':other,'steps':steps,'terminal':{'kind':'forced_proper_cycle','cycle':cycle,'selected_edges':sorted(selected)}}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    receipt=json.loads((HERE/'marked_cut_search.json').read_text());E=edges(25)
    need(sum(18 in e for e in E)==1,'mandatory endpoint is a leaf')
    need(sorted(row['other_endpoint'] for row in receipt['endpoints'])==[v for v in range(1,26) if v!=18],'complete endpoint universe')
    cuts=[];residual=[]
    for row in receipt['endpoints']:
        t=row['other_endpoint']
        if not row['ruled_out']:residual.append(t);continue
        S=set(row['minimal_violating_cut']);cc=components(set(range(1,26))-S,E);budget=len(S)+1-len(cc)
        obligation=sum(v in S for v in (18,t))+int({1,15}<=S)
        need(obligation>budget,('endpoint-and-edge cut violation',t))
        cuts.append({'other_endpoint':t,'cut':sorted(S),'outside_components':cc,'budget':budget,'required_spend':obligation})
    need(residual==[8,10,21],'exact cut survivors')
    # Independent validation of every positive unmarked tight-cut certificate.
    unmarked=json.loads((HERE/'tight_cut_search.json').read_text());positive=[]
    for row in unmarked['queries']:
        S=set(row['minimal_tight_cut']);n=row['n']
        if not S:continue
        need(set(row['edge'])<=S,'queried edge lies inside cut')
        cc=components(set(range(1,n+1))-S,edges(n));need(len(cc)==len(S)+1,'tight-cut components')
        positive.append({'n':n,'edge':row['edge'],'cut':sorted(S),'outside_components':cc})
    out={'status':'PROVED by explicit finite certificate, independently replayed: G25 has no Hamiltonian path containing edge1--15','cut_eliminated_endpoints':cuts,'residual_traces':[propagate(t) for t in residual],'positive_unmarked_cut_certificates':positive,'checks':checks}
    out=json.loads(json.dumps(out))
    if args.write:(HERE/'guarded_forbidden.json').write_text(json.dumps(out,indent=2)+'\n')
    else:need(json.loads((HERE/'guarded_forbidden.json').read_text())==out,'frozen guarded certificate')
    print('guarded forbidden checks',checks)

if __name__=='__main__':main()
