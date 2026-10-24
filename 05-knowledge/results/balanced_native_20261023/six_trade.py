#!/usr/bin/env python3
"""Exact small-move lattice certificate. No symbolic algebra dependency in verification."""
from pathlib import Path
from itertools import combinations
from collections import defaultdict,Counter,deque
from fractions import Fraction as F
import json,subprocess
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(m)
def determinant(A):
    A=[list(row) for row in A];n=len(A);previous=1;sign=1
    for k in range(n-1):
        p=next((i for i in range(k,n) if A[i][k]),None)
        if p is None:return 0
        if p!=k:A[k],A[p]=A[p],A[k];sign=-sign
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                z=A[i][j]*pivot-A[i][k]*A[k][j]
                need(z%previous==0,'fraction-free determinant division is exact');A[i][j]=z//previous
            A[i][k]=0
        previous=pivot
    return sign*A[-1][-1]
def main():
    cert=json.loads((HERE/'six_trade_basis.json').read_text());ids=cert['sources']
    path='0f9dec9c3e1d:05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json'
    inv=json.loads(subprocess.check_output(['git','show',path],cwd=ROOT));opts={r['source']:r for r in inv['records']}
    need(ids==sorted(opts),'complete original28-source supply')
    lab=[(sum(b<<j for j,b in enumerate(opts[s]['source_label'][:4])),opts[s]['source_label'][4]) for s in ids]
    def total(I):
        p=x=0
        for i in I:p^=lab[i][0];x+=lab[i][1]
        return p,x%128
    moves=set()
    for k in (2,3):
        bank=defaultdict(list)
        for I in combinations(range(28),k):bank[total(I)].append(I)
        for bucket in bank.values():
            for a,b in combinations(bucket,2):
                v=tuple(int(i in a)-int(i in b) for i in range(28))
                v=v if next(x for x in v if x)>0 else tuple(-x for x in v);moves.add(v)
    need(sorted(moves)==list(map(tuple,cert['trades'])),'complete four/six trade catalogue')
    for v in moves:
        need(sum(v)==0 and set(v)<= {-1,0,1},'binary native endpoint trade, exact target count')
        need(total(i for i,x in enumerate(v) if x==1)==total(i for i,x in enumerate(v) if x==-1),'every trade has the full paired charge')
    columns=cert['basis'];need(len(columns)==27 and all(len(c)==27 for c in columns),'rank27 basis dimensions')
    for c,expr in zip(columns,cert['basis_trade_combinations']):
        rebuilt=[sum(int(coef)*cert['trades'][int(j)][i] for j,coef in expr.items()) for i in range(28)]
        need(rebuilt[:27]==c and sum(rebuilt)==0,'each basis vector has an exact signed native-trade expression')
    B=list(map(list,zip(*columns)));det=abs(determinant(B));need(det==cert['determinant_index']==2048,'literal determinant of actual trade basis')
    generators=[(p^lab[-1][0],(x-lab[-1][1])%128) for p,x in lab[:-1]]
    subgroup={(0,0)};todo=deque(subgroup)
    while todo:
        p,x=todo.popleft()
        for q,y in generators:
            z=(p^q,(x+y)%128)
            if z not in subgroup:subgroup.add(z);todo.append(z)
    need(len(subgroup)==det,'charge kernel index equals the generated native-trade lattice index')
    # Exact formal coordinates of the marked primitive18 difference; not a legal route order.
    plus={176,1394,1457,1477,1532,1923,2087,2176,3369};minus={1425,1549,1841,1921,2185,2553,2835,2905,3000}
    delta=[int(s in plus)-int(s in minus) for s in ids];M=[[F(x) for x in row]+[F(y)] for row,y in zip(B,delta)]
    for k in range(27):
        p=next(i for i in range(k,27) if M[i][k]);M[k],M[p]=M[p],M[k];z=M[k][k];M[k]=[x/z for x in M[k]]
        for i in range(27):
            if i!=k:
                z=M[i][k];M[i]=[x-z*y for x,y in zip(M[i],M[k])]
    coords=[r[-1] for r in M];need(all(x.denominator==1 for x in coords),'integer formal charge decomposition')
    need([sum(B[i][j]*coords[j] for j in range(27)) for i in range(27)]==delta[:27],'exact marked formal endpoint')
    return {'status':'FINITE-EXACT lattice equality for original28-source supply; feasible route inference excluded',
      'moves':len(moves),'supports':dict(Counter(sum(x!=0 for x in v) for v in moves)),
      'basis_rank':27,'basis_determinant':det,'charge_image_order':len(subgroup),
      'primitive18_formal_coordinates':[int(x) for x in coords],'formal_coefficient_L1':sum(abs(int(x)) for x in coords),'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
