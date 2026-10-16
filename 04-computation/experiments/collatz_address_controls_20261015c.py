#!/usr/bin/env python3
"""Exact controls for address folds, signed graph quotients, and Pell switching."""
import itertools
import json
from collections import Counter
from pathlib import Path
import sympy as s
from collatz_gram_slack_20261015b import close_graph, gram, need, transport_slack

ROOT=Path(__file__).resolve().parents[2]


def main():
    bits=list(itertools.product((0,1),repeat=4))
    fold=lambda v:(v[0]+v[1],v[2],v[3])
    fibres=Counter(map(fold,bits))
    need(len(bits)==16 and len(fibres)==12,'address quotient cardinalities')
    need([sum(1 for v in bits if v[0]+v[1]==j) for j in range(3)]==[4,8,4],'unfolded selected-pair sectors')
    proper=[v for v in bits if 0<sum(v)<4]
    need([sum(sum(v)==j for v in proper) for j in (1,2,3)]==[4,6,4],'Boolean proper rank sectors')
    qproper=[u for u in fibres if 0<sum(u)<4]
    need([sum(sum(u)==j for u in qproper) for j in (1,2,3)]==[3,4,3],'divisor rank quotient')
    need([sum(u[0]==j for u in qproper) for j in range(3)]==[3,4,3],'different repeated-prime exponent sectors')
    need(all(fold(tuple(1-x for x in v))==(2-fold(v)[0],1-v[2],1-v[3]) for v in bits),'complement intertwines')
    need(all(fibres[u]==(2 if u[0]==1 else 1) for u in fibres),'retained multiplicity sidecar')
    need(sum(fibres[u] for u in qproper)==14,'weighted quotient recovers proper source count')

    # Same quotient and occupancy, but opposite-sign and same-sign row weights.
    fold_test=[]
    for row in ({0:1,1:-1},{0:1,1:1}):
        G=gram([row],[0,1])
        _,_,loss=transport_slack(G,[0,1],{0:(0,1),1:(0,1)})
        fold_test.append(loss)
    need([z['new_norm'] for z in fold_test]==[0,4],'occupancy is not paid cancellation')

    # All 125 maps on three nonroot sources, and all six orders of their edges.
    sources=[1,2,3,4]
    rows=[{1:1,2:2,3:-1},{1:-1,2:1,3:1,4:-2}]
    G=gram(rows,sources)
    graph_orders=0
    for targets in itertools.product(range(1,6),repeat=3):
        full=dict(zip((2,3,4),targets))
        cc,_=close_graph(full,sources)
        _,_,optimal=transport_slack(G,sources,{e:cc[e][:2] for e in sources})
        for order in itertools.permutations((2,3,4)):
            edges={}
            last=sum(abs(v) for r in G for v in r)
            for e in order:
                edges[e]=full[e]
                cache,_=close_graph(edges,sources)
                _,_,loss=transport_slack(G,sources,{v:cache[v][:2] for v in sources})
                need(loss['new_norm']<=last,'every finite edge order is monotone')
                last=loss['new_norm']
            need(last==optimal['new_norm'],'same final quotient in every order')
            graph_orders+=1

    # The new source refinement uses only this inherited changing-generator control.
    u=s.symbols('u',integer=True)
    B=s.Matrix([[0,-1],[1,1]])
    U=s.Matrix([[1,u],[0,1]])
    C=U.inv()*B.inv()*U
    need(B**3==-s.eye(2) and C**3==-s.eye(2),'finite-order individual generators')
    product=s.simplify(B*C)
    need(product==s.Matrix([[1,u],[u,u*u+1]]),'changing-generator product')
    V=s.Matrix([[1,-1],[0,1]])
    need(V.inv()*product.subs(u,2)*V==s.Matrix([[3,4],[2,3]]),'exact Q-Pell conjugacy')
    K=s.Matrix([[2,1],[1,2]])
    need(B.T*K*B==K,'first metric')
    KU=U.T*K*U
    need(s.simplify(C.T*KU*C-KU)==s.zeros(2),'second metric')

    A=s.Matrix([[1,-2],[0,1]])
    r=s.Matrix([[1,2]])
    need(r*A==s.Matrix([[1,0]]),'unimodular change hides coefficient')
    need(max((r*s.Matrix(z))[0]**2 for z in itertools.product((-1,1),repeat=2))==9,'old cube energy')
    need(max(((r*A)*(A.inv()*s.Matrix(z)))[0]**2 for z in itertools.product((-1,1),repeat=2))==9,'transported cube retains energy')

    # The halving relation alone retains a linear-drift hostile.
    chi=lambda n:(0,1,-1)[n%3]
    H=[0]
    for n in range(1,16385):
        need(chi(2*n)==-chi(n),'hostile preserves every halving edge')
        H.append(H[-1]+(2 if n%3==2 else -1)*chi(n))
        need(H[-1]==-n-(n%3==2),'exact hostile prefix discrepancy')
    for n in range(1,4097):
        need(H[4*n]-H[n]==-3*n,'halving-only annular drift')

    out={'status':'FINITE-EXACT algebra and hostile controls; no Collatz asymptotic claim',
         'address_fold':{'unfolded':16,'selected_pair_layers':[4,8,4],'proper_rank_layers':[4,6,4],
                         'quotient':12,'proper_quotient_rank_layers':[3,4,3],
                         'middle_coordinate_fibre_weight':2},
         'fold_weight_controls':fold_test,'toy_graph_orders':graph_orders,
         'Pell_product':str(product),'basis_cube_hostile':{'original':9,'false_fresh_cube':1,'transported_cube':9},
         'halving_only_hostile':{'checked_annuli':4096,'Delta':'-3N','first_failed_odd_source':3}}
    (ROOT/'05-knowledge/results/collatz_address_controls_20261015c.json').write_text(json.dumps(out,indent=2)+'\n')
    print('address, Gram, 750 graph-order, Pell, and cube checks passed')


if __name__=='__main__':main()
