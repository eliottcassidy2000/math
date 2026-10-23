#!/usr/bin/env python3
"""Independent exact Smith checks of the two SYNTHETIC order-five switches."""
import json
from pathlib import Path
from collections import Counter
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from check_metagraph_independent import laplacian,binary_kernel,determinant,check
HERE=Path(__file__).resolve().parent

def edges(table):
    out=Counter()
    for a,b,w in table:out[tuple(sorted((a,b)))]+=w
    return out

def degrees(n,E):return [sum(w*((i==a)+(i==b)) for (a,b),w in E.items()) for i in range(n)]

def inspect(n,B,K):
    E=B.copy()
    for e,w in K.items():E[e]+=2*w
    M=laplacian(n,E);A=[r[:-1] for r in M[:-1]]
    D=smith_normal_form(S.Matrix(A),domain=S.ZZ);d=[abs(int(D[i,i])) for i in range(n-1)]
    check(all(d),'graph disconnected');v=[(x&-x).bit_length()-1 for x in d if x%2==0]
    _,basis=binary_kernel(A);V=S.Matrix(basis);raw=V*S.Matrix(A)*V.T
    Q=raw.applyfunc(lambda x:(int(x)//2)%2);rank,_=binary_kernel(Q.tolist())
    return dict(tree_count=determinant(A),positive_two_valuations=v,two_rank=len(basis),half_form_rank=rank,half_form_matrix=[[int(x) for x in row] for row in Q.tolist()]),M

def run():
    source=json.loads((HERE.parent/'metagraph/dyadic_network_receipt.json').read_text())[2]
    controls=json.loads((HERE.parent/'metagraph/bockstein_receipt.json').read_text());n=source['nodes'];B=edges(source['blue_edges']);K=edges(source['paired_black_edges'])
    actual,M=inspect(n,B,K)
    blue_control=controls['degree_preserving_rewiring'];BB=edges(blue_control['rewired_blue_edges'])
    check(degrees(n,B)==degrees(n,BB),'blue degrees changed');blue_switch,MB=inspect(n,BB,K)
    check(blue_switch['tree_count']==37184 and blue_switch['two_rank']==4,'blue switch claim failed')
    black_control=controls['blue_fixed_black_rewiring'];KK=edges(black_control['rewired_paired_black_edges'])
    check(degrees(n,K)==degrees(n,KK),'paired-black degrees changed');black_switch,MK=inspect(n,B,KK)
    check(all((M[i][j]-MK[i][j])%2==0 for i in range(n) for j in range(n)),'mod2 graph changed')
    check(black_switch['tree_count']==37952,'black switch determinant failed')
    check(black_switch['positive_two_valuations']==[1,1,4] and black_switch['half_form_rank']==2,'black switch torsion failed')
    check([actual['half_form_matrix'][i][i] for i in range(3)]==[black_switch['half_form_matrix'][i][i] for i in range(3)]==[1,1,1],'half-form diagonals disagree')
    check(actual['tree_count']==41096 and actual['positive_two_valuations']==[1,1,1],'actual graph mismatch')
    out=dict(status='PASS',scope='Synthetic colored multigraph controls, not alternative tournament realizations',actual=actual,blue_switch=blue_switch,black_switch=black_switch)
    (HERE/'switch_controls_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':run()
