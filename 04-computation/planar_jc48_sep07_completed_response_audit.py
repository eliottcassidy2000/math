#!/usr/bin/env python3
"""Independent Fraction reconstruction of the complete source-image matrix.

Producer output is read only as a proposed kernel certificate and is checked
against this independently assembled matrix. No producer is imported.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
from hashlib import sha256
import json

GATES = 0
def check(x, label):
    global GATES
    GATES += 1
    if not x:
        raise RuntimeError(label)

def add(*polys):
    r = {}
    for f in polys:
        for m,c in f.items(): r[m] = r.get(m,F(0))+c
    return {m:c for m,c in r.items() if c}

def mul(f,g):
    r = {}
    for (a,b),c in f.items():
        for (i,j),d in g.items():
            m=(a+i,b+j); r[m]=r.get(m,F(0))+c*d
    return {m:c for m,c in r.items() if c}

def scale(f,c): return {m:v*c for m,v in f.items() if v*c}
def deriv(f,k):
    r={}
    for m,c in f.items():
        if m[k]:
            n=list(m); n[k]-=1; r[tuple(n)]=c*m[k]
    return r

def literal(f):
    r={}
    for (a,b),c in f.items():
        for j in range(a+b+1):
            n=a+2*b+j
            if n<=15:
                key=(n,b+2*j); r[key]=r.get(key,F(0))+c*comb(a+b,j)
    return {m:c for m,c in r.items() if c}

def rref(a):
    a=[row[:] for row in a]; piv=[]; row=0
    for col in range(len(a[0])):
        j=next((j for j in range(row,len(a)) if a[j][col]),None)
        if j is None: continue
        a[row],a[j]=a[j],a[row]
        z=a[row][col]; a[row]=[v/z for v in a[row]]
        for j in range(len(a)):
            if j!=row and a[j][col]:
                z=a[j][col]; a[j]=[v-z*w for v,w in zip(a[j],a[row])]
        piv.append(col); row+=1
        if row==len(a): break
    return a,piv

D={(3,0):F(1),(0,2):F(-1)}
H={(1,0):F(-3),(2,0):F(8,3),(3,0):F(-1376,135),
   (4,0):F(896,15),(0,2):F(-32,5),(5,0):F(-731648,2025),
   (1,2):F(512,75)}
RM=[(a,b) for a in range(12) for b in range(6) if 6<=a+2*b<=11]
TM=[(a,b) for a in range(16) for b in range(8)
    if 11<=a+2*b<=15 and 2*a+3*b<=22]+[(3,6)]
RM.sort(key=lambda m:(m[0]+2*m[1],m[1]))
TM.sort(key=lambda m:(m[0]+2*m[1],m[1]))
check(len(RM)==30 and len(TM)==17,'independently generated universes')
cols=[]
for a,b in RM:
    R={(a,b):F(1)}
    Ar=add(mul(add(scale(D,2),{(3,0):F(3)}),R),
           mul(mul({(1,0):F(1)},D),deriv(R,0)))
    Br=mul({(1,0):F(1)},add(mul({(0,1):F(-2)},R),mul(D,deriv(R,1))))
    psi=add({(a+3,b+1):F(-(10+2*a+3*b),2)},
            mul(D,add(mul(Ar,deriv(H,1)),scale(mul(Br,deriv(H,0)),-1))))
    cols.append(literal(psi))
cols += [scale(literal({m:F(1)}),-1) for m in TM]
keys=sorted(set().union(*(set(col) for col in cols)))
matrix=[[col.get(k,F(0)) for col in cols] for k in keys]
red,piv=rref(matrix)
check(len(piv)==36,'independent rational elimination rank36')
free=[j for j in range(47) if j not in piv]
check(free==list(range(24,30))+[31,37,42,44,46],'independent complete free columns')
check(all(not any(cols[j].values()) for j in range(24,30)), 'all six terminal source columns vanish')
check(len(rref([row[:30] for row in matrix])[1])==24,'source action kernel exactly six')

root=Path(__file__).resolve().parents[1]
report=root/'05-knowledge/results/planar_jc48_sep07_completed_response.out'
lines=report.read_text().splitlines()
def bank(name):
    raw=next(line[len(name)+2:] for line in lines if line.startswith(name+': '))
    return [[F(v) for v in row] for row in json.loads(raw)]
SB,LB=bank('SOURCE_BASIS'),bank('CARRIER_LIFT_BASIS')
check(len(SB)==17 and len(LB)==30, 'typed certificate shapes')
for j in range(5):
    vector=[row[j] for row in LB+SB]
    check(all(sum(a*b for a,b in zip(row,vector))==0 for row in matrix),
          'independent full substitution of each source and lift')
for i in range(5):
    check(SB[[1,7,12,14,16][i]]==[F(int(i==j)) for j in range(5)],
          'source identity coordinates')
det=LB[1][3]*LB[3][4]-LB[1][4]*LB[3][3]
check(det==F(1751787,413887398592), 'independent row10 obstruction determinant')
Rstar={(0,3):F(12015,15962432),(1,3):F(28467,3990608),
       (2,3):F(2499,498826),(3,3):F(-4806,1247065),
       (4,3):F(639368,11223585),(0,5):F(532468,6235325)}
Tstar={(7,2):F(108135,15962432),(3,4):F(-12015,1680256),
       (8,2):F(208143,3990608),(4,4):F(-741987,7981216),
       (5,4):F(-196239,997652),(1,6):F(180225,15962432),
       (2,6):F(346905,3990608),(3,6):F(-1)}
v=[Rstar.get(m,F(0)) for m in RM]+[Tstar.get(m,F(0)) for m in TM]
check(all(sum(a*b for a,b in zip(row,v))==0 for row in matrix),
      'independent exact six-term supplier')
old={(5,4):F(-27945,235202),(1,6):F(39123,470404),
     (2,6):F(52578,117601),(3,6):F(-1)}
relations=[row[30:] for row in red if not any(row[:30]) and any(row[30:])]
check(len(relations)==12,'independent twelve source conditions')
check(any(sum(a*old.get(m,F(0)) for a,m in zip(row,TM)) for row in relations),
      'old packet independent obstruction')
blob=json.dumps([[str(v) for v in row] for row in red],separators=(',',':')).encode()
print('Independent completed-response Fraction image audit PASS')
print('Nonzero coefficient rows:',len(keys))
print('Rank: 36; source dimension: 5; coordinate kernel: 6')
print('Always-active gates:',GATES)
print('RREF SHA256:',sha256(blob).hexdigest())
