#!/usr/bin/env python3
"""Refine synchronized bits until the union graph has exactly the decoded cube."""
import json,itertools
from pathlib import Path
from edge_atlas import Atlas,edge,edges,need,OPS,arity
HERE=Path(__file__).resolve().parent

def leaf_index(n,word):
    for s in word:
        q=OPS[s][0];n=(n-(q-1)//2)//q
    return n
def refined_schedule(n,word):
    m=leaf_index(n,word)
    while m>=2812:
        word+='A';m=(m-12)//25
    need(112<=m<=2811,'fully refined finite leaf range')
    return word
def union_graph(cube):return cube.base|{e for removed,added in cube.moves.values() for e in added}

def undirected_hamilton_cycles(n,es):
    adj=[set() for _ in range(n)]
    for a,b in es:adj[a].add(b);adj[b].add(a)
    dp={(1,0):1}
    for used in range(1,n):
        new={}
        for (mask,v),ways in dp.items():
            for z in adj[v]:
                if not(mask>>z&1):
                    key=(mask|(1<<z),z);new[key]=new.get(key,0)+ways
        dp=new
    return sum(ways for (mask,v),ways in dp.items() if 0 in adj[v])//2

def main():
    # Enumerate only local variable edges; vertices on internal common paths
    # force one common edge and hence exactly one variable edge at their ports.
    local=[]
    for kind,var,forced,vertices in [
        ('reverse',[(0,1),(1,3),(3,2),(2,0)],{1,2},4),
        ('exchange',[(0,1),(1,4),(4,5),(5,2),(2,3),(3,0)],{1,2,3,4},6)]:
        cases=[]
        for mask in range(1<<len(var)):
            chosen=[e for j,e in enumerate(var) if mask>>j&1];deg=[0]*vertices
            for u,v in chosen:deg[u]+=1;deg[v]+=1
            if not all(1<=x<=2 for x in deg) or not all(deg[v]==1 for v in forced):continue
            cases.append({'edges':chosen,'degrees':deg})
        need(len(cases)==(2 if kind=='reverse' else 3),'complete local degree possibilities')
        local.append({'kind':kind,'feasible_degree_cases':cases})
    reverse={(0,1),(1,3),(3,2),(2,0),(1,4),(4,2),(0,5),(5,3)}
    exchange={(0,1),(1,4),(4,5),(5,2),(2,3),(3,0),(1,6),(6,7),(7,2),(3,8),(8,9),(9,4),(5,10),(10,0)}
    need(undirected_hamilton_cycles(6,reverse)==2,'reversal gadget exact Hamilton count')
    need(undirected_hamilton_cycles(11,exchange)==2,'exchange gadget excludes its isolated internal cycle')
    bank=json.loads((HERE.parent/'square_sum_spectrum_overnight_20261031e/phase_switches.json').read_text())
    types={};minimum_exchange=10**9
    for n,moves in bank['moves'].items():
        for move in moves:
            types[move[0]]=types.get(move[0],0)+1
            if move[0]=='reverse':need(move[2]-move[1]>=2,'nonempty reversal interior')
            else:
                minimum_exchange=min(minimum_exchange,move[2]-move[1],move[3]-move[2])
                need(min(move[2]-move[1],move[3]-move[2])>=4,'both exchange interiors have forced degree-two vertices')
    a=Atlas();coarse=a.cycle(70312,'A');word=refined_schedule(70312,'A');fine=a.cycle(70312,word)
    need(word=='AA' and coarse.bits==25 and fine.bits==625,'saturation sizes')
    need(coarse.base==fine.base and union_graph(coarse)==union_graph(fine),'same union graph after splitting ties')
    q=[(min(fine.moves[0][0]),True),(min(fine.moves[1][1]),True)]
    need(coarse.compile(q) is None,'coarse synchronized bit rejects the local query')
    compiled=fine.compile(q);need(compiled is not None,'refined graph permits independent local changes')
    es=edges(fine.walk(compiled['one_code']),True)
    need(es<=union_graph(coarse) and all((e in es)==wanted for e,wanted in q),'native graph hostile witness')
    print(json.dumps({'status':'passed','local_exhaustive_cases':local,'exact_small_gadget_cycles':[2,2],
                      'complete_leaf_move_types':types,'minimum_exchange_block_length':minimum_exchange,
                      'tie_erasure_hostile':{'N':70312,'coarse_schedule':'A','refined_schedule':word,'coarse_bits':25,'refined_bits':625,
                        'same_union_graph':True,'queries':q,'coarse_accepts':False,'refined_accepts':True,'refined_free_bits':compiled['free_bits']},
                      'scope':'Fully refined union graph has exactly the independent-switch Hamiltonian cycles. Coarse synchronized codes are not closed under its local edge queries.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
