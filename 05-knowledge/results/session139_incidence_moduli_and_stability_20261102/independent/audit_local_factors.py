#!/usr/bin/env python3
"""Independent symbolic exact audit of conditional local-factor hostiles."""
from fractions import Fraction as F
import json
from audit_bridges import P,det,check
x=P({(1,0):1});K=[[P() for _ in range(5)] for _ in range(5)];K[0][1]=x;K[1][0]=-x
U=[[F(1,2),F(1,2),F(0)],[F(-1,4),F(1,4),F(1,2)]]
for i in range(2):
    for j in range(3):K[i][j+2]=P(U[i][j]);K[j+2][i]=-P(U[i][j])
for i in range(2,5):
    for j in range(i+1,5):K[i][j]=P(1);K[j][i]=P(-1)
A=[[(int(i!=j)+K[i][j])*F(1,2) for j in range(5)] for i in range(5)]
u=det([[A[i][j]+int(i==j) for j in range(5)] for i in range(5)]);da=det(A);f=u-da
check(f==(x*x-x*F(3,8)+F(81,64))*F(1,4),'symbolic n5 margin')
z=F(3,4);p=(1+z)/2;q=p*(1-p);val=f.at(z,0);lo=f.at(-1,0);hi=f.at(1,0);L=(1-p)*lo+p*hi
check(all(da.at(k,0)>0 for k in (-1,z,1)),'minimum margin branch')
check(val-(1-2*q)*L==F(-7,8192),'n5 conditional factor hostile')
# Block determinant certificate: pair det=-q, det(I+pair)=1-q;
# each directed triangle det=1,det(I+C3)=2. Product gives full ordinal matrix.
q8=F(1,4);u8=4*(1-q8);da8=-q8;minimum=u8-abs(da8);L8=F(4)
check(minimum-(1-q8)*L8==-q8,'n8 conditional factor hostile')
check(minimum>1-q8,'absolute matching floor not refuted')
print(json.dumps({'status':'PASS','n5_symbolic_fplus':str(f),'n5_residual':str(val-(1-2*q)*L),'n8_residual':str(minimum-(1-q8)*L8),'scope':'conditional factors only; absolute floor survives'},sort_keys=True,indent=2))
