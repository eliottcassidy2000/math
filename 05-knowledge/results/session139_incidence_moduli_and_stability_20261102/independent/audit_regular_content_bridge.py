#!/usr/bin/env python3
"""Exact bounded check of the R11/Paley11 coarse-spectrum/content bridge."""
from itertools import combinations
from pathlib import Path
import hashlib,json

def check(ok,label):
    if not ok:raise RuntimeError(label)

def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def transpose(a):return [list(col) for col in zip(*a)]
def graph(S):return [[int((j-i)%11 in S) for j in range(11)] for i in range(11)]

def main():
    records=[]
    for name,S in [('R11',set(range(1,6))),('P11',{1,3,4,5,9})]:
        A=graph(S);At=transpose(A)
        check(all(A[i][i]==0 for i in range(11)) and all(A[i][j]+A[j][i]==1 for i in range(11) for j in range(11) if i!=j),('tournament',name))
        check(all(sum(row)==5 for row in A),('regular',name))
        check(all(sum(A[i][j]-At[i][j] for j in range(11))==0 for i in range(11)),('skew fixes ones',name))
        check(mm(A,At)==mm(At,A),('normality',name))
        triangles=sum(A[i][j]*A[j][k]*A[k][i]+A[i][k]*A[k][j]*A[j][i] for i,j,k in combinations(range(11),3))
        check(triangles==55,('triangles',name))
        A2=mm(A,A);A4=mm(A2,A2);trace4=sum(A4[i][i] for i in range(11))
        records.append({'name':name,'order':11,'outdegree':5,'cyclic_triangles':triangles,'trace_A4':trace4,'exact_normality':True,'skew_annihilates_ones':True})
    check([r['trace_A4'] for r in records]==[880,660],'noncospectral hostile')
    P=graph({1,3,4,5,9});G={tuple((a*x+b)%11 for x in range(11)) for a in {1,3,4,5,9} for b in range(11)}
    check(len(G)==55,'affine distinctness')
    check(all(P[i][j]==P[g[i]][g[j]] for g in G for i in range(11) for j in range(11)),'affine automorphisms')
    check(all(tuple(g[h[x]] for x in range(11)) in G for g in G for h in G),'affine closure')
    receipt=Path(__file__).resolve().parents[1]/'all_moduli/content.normal.json'
    raw=receipt.read_bytes();expected='5ee1af479cff6ecef4a34874a94768e3d86432d24da7bed4e0c8d0dc28c18605'
    check(hashlib.sha256(raw).hexdigest()==expected,'audited content receipt pin')
    gR=json.loads(raw)['R11']['context_content'];check(gR==33 and gR%5!=0,'R11 content separation')
    print(json.dumps({'status':'PASS','graphs':records,'R11_context_content':gR,'Paley11_certified_context_content_divisor':55,'Paley11_exact_content_claimed':False,'cospectral':False,'content_receipt_sha256':expected},sort_keys=True,indent=2))
if __name__=='__main__':main()
