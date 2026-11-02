#!/usr/bin/env python3
"""Independent exact check of the order-eight naming correction, not the census."""
import json

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def pfaffian(a):
    if not a:
        return 1
    return sum((-1)**(j+1)*a[0][j]*pfaffian([
        [a[k][l] for l in range(len(a)) if l not in (0,j)]
        for k in range(len(a)) if k not in (0,j)])
        for j in range(1,len(a)))

qr={1,2,4}
s=[[0 if i==j else (1 if i==0 or (j!=0 and (j-i)%7 in qr) else -1)
    for j in range(8)] for i in range(8)]
require(all(s[i][j]==-s[j][i] for i in range(8) for j in range(8)), 'skew')
gram=[[sum(s[i][k]*s[j][k] for k in range(8)) for j in range(8)] for i in range(8)]
require(gram==[[7*(i==j) for j in range(8)] for i in range(8)], 'SS^T=7I')
pf=pfaffian(s)
require(pf==-49, 'Pfaffian orientation')
scores=[row.count(1) for row in s]
require(scores==[7]+[3]*7, 'source plus regular Paley7 scores')
require(len(set(scores))!=1, 'order8 tournament is not regular')
require(pf*pf==7**4, 'Hadamard equality determinant')
print(json.dumps({'status':'PASS','construction':'Paley7 QR={1,2,4}, augmented by a source',
    'matrix':s,'gram_diagonal':7,'pfaffian':pf,'scores':scores,
    'scope':'Naming/witness correction only; original finite census unchanged.'},indent=2,sort_keys=True))
