#!/usr/bin/env python3
"""Exhaust the blue slices of all six failed affine cells; replay literal repairs."""
import importlib.util
import itertools
import json
from pathlib import Path
from search_repairs import arc, directed, blue, edge, move

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('low_orders',HERE.parent/'port_residual_integral_20261025/low_orders.py')
low=importlib.util.module_from_spec(spec);spec.loader.exec_module(low)
bank=low.bank

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def reframe(T,q):
    return tuple(sum(arc(T,q[u],q[v])<<v for v in range(len(T)) if u!=v) for u in range(len(T)))

def replay(T,qs):
    n=len(T);P=tuple(range(n));target=move(T,P);xor=frozenset()
    for q in qs:
        need(sorted(q)==list(P),'permutation')
        need(directed(T,q),'directed path')
        need(not blue(T,q),'black path')
        xor^=edge(q);T=move(T,q)
    need(T==target and xor==edge(P),'labelled target and XOR')

def main():
    inherited=json.loads((HERE.parent/'port_residual_integral_20261025/low_orders.normal.json').read_text())
    results=[];total=0
    for row in inherited['failed_ledger_hostiles']:
        n=row['n'];I=set(row['ports']);opened=row['middle_open'];chord=row['x_to_c']
        _,qs=bank.words(n,I,opened,chord)
        A=bank.symbolic_source(n,I,[(1,n-2,int(opened)),(2,n-2,int(chord))])
        current=A
        for q in qs[:row['blue_step']-1]:current=bank.k_move(current,q)
        basis=bank.eliminate(bank.blue_equations(current,qs[row['blue_step']-1]))
        need(basis is not None,'blue slice feasible')
        R=[[bank.reduce_form(f,basis) for f in a] for a in A]
        mask=0
        for a in R:
            for f in a:mask|=f&~1
        bits=[i for i in range(1,mask.bit_length()) if mask>>i&1]
        need(len(bits)==(3 if n==12 else 0),'exact blue-slice dimension')
        sources=set();certs=[]
        for values in itertools.product((0,1),repeat=len(bits)):
            assignment=sum(v<<i for v,i in zip(values,bits))
            T=tuple(sum(((f&1)^((f&assignment).bit_count()&1))<<v for v,f in enumerate(a) if u!=v) for u,a in enumerate(R))
            need(T not in sources,'distinct free assignments');sources.add(T)
            low.literal_check(T,row['ports'],opened,chord,qs,row['literal_blue_flags'])
            if n==12:
                frame=qs[0];U=reframe(T,frame)
                newI={v for v in range(1,n-1) if arc(U,v,0)}
                need(arc(U,n-1,0) and U[0].bit_count()==n-5 and U[-1].bit_count()==4,'closed same endpoint layer')
                result=bank.words(n,newI,bool(arc(U,1,n-2)),bool(arc(U,2,n-2)))
                need(result is not None,'new observer has base repair')
                name,local=result
                repair=[[frame[v] for v in q] for q in local]+qs[1:]
                need(len(repair)==5,'five black stages')
                extra={'reframed_ports':sorted(newI),'branch':name}
            else:
                sat=json.loads((HERE/f"sat_repair_{len([r for r in results if r['n']==14])}.json").read_text())
                need(sat['solver_result']=='sat' and tuple(sat['source_rows'])==T,'SAT witness for exact unique source')
                repair=sat['words'];extra={}
            replay(T,repair)
            # Converse is transported with path reversal; no port is identified.
            converse=tuple(((1<<n)-1)^(1<<u)^T[u] for u in range(n))
            P=tuple(range(n));rev=P[::-1]
            # Relabel by reversed original observer so replay again targets its P move.
            Tc=reframe(converse,rev)
            rc=[[n-1-v for v in q[::-1]] for q in repair]
            replay(Tc,rc)
            certs.append({'source_rows':T,'words':repair,'converse_replay':True,**extra});total+=1
        need(tuple(row['source_rows']) in sources,'inherited hostile included')
        results.append({'n':n,'ports':row['ports'],'blue_step':row['blue_step'],
                        'free_arc_bits':bits,'blue_slice_size':len(sources),'certificates':certs})
    need(total==20,'complete exceptional universe')
    print(json.dumps({'status':'FINITE-EXACT PASS','exceptional_sources':total,
        'universe':'all completions in the offending blue slice of each of the six inherited failed cells',
        'method':'exact GF2 elimination followed by independent literal arc replay; no solver trust in verification',
        'labelled_targets_preserved':True,'converse_checked':True,'results':results},indent=2,sort_keys=True))

if __name__=='__main__':main()
