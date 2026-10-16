#!/usr/bin/env python3
"""Exact scoped audit of incoming THM-4658: clocks, energy seam and branch cross term."""
from fractions import Fraction as F
from pathlib import Path
import json

def need(p,msg):
    if not p:raise AssertionError(msg)

def T(n):return (3*n+1)//2 if n%2 else n//2

N=2**16
t=[0]*(N+1);t[1]=1
for n in range(2,N+1):
    x=n;steps=0
    while x>=n:
        x=T(x);steps+=1
        need(steps<10000,'independent finite descent cap')
    t[n]=(-1 if steps%2 else 1)*t[x]
C3=[[0]*(N+1) for _ in range(3)]
C9=[[0]*(N+1) for _ in range(9)]
for n in range(1,N+1):
    for r in range(3):C3[r][n]=C3[r][n-1]+(t[n] if n%3==r else 0)
    for r in range(9):C9[r][n]=C9[r][n-1]+(t[n] if n%9==r else 0)
H=[3*C3[2][n]-sum(C3[r][n] for r in range(3)) for n in range(N+1)]
def channel2(r,n):return 3*C9[r][n]-C3[r%3][n]
for n in range(N//6+1):
    need(H[4*n]==H[n]+channel2(2,3*n)-channel2(8,6*n),'exact two-octave residue identity')
energies=[]
for k in range(14):
    X=2**k
    E=sum(v*v for v in H[X:2*X])
    E4=sum(v*v for v in H[4*X:8*X])
    delta=[H[4*n]-H[n] for n in range(X,2*X)]
    D=sum(v*v for v in delta)
    C=sum(H[n]*delta[n-X] for n in range(X,2*X))
    R=E4-4*E-8*C-4*D
    need(abs(R)<=48*X or (abs(R)-48*X)**2<=576*X*E4,'integer audit of remainder bound')
    energies.append({'k':k,'E':E,'E4':E4,'D':D,'C':C,'R':R})
f=[-1,-2,2,-1,-2,2,-1,1,2]
need(sum(f)==0,'mean-zero function')
u=[f[(2*y)%9] for y in range(27)]
v=[f[((2*y-1)//3)%9] if y%3==2 else 0 for y in range(27)]
norm=F(sum(f[y%9]**2 for y in range(27)),27)
norm_u=F(sum(z*z for z in u),27)
norm_v=F(sum(z*z for z in v),27)
cross=F(sum(a*b for a,b in zip(u,v)),27)
need((norm,norm_u,norm_v,cross)==(F(8,3),F(8,3),F(8,9),F(2,3)),'exact branch cross-term witness')
values=[]
for mask in range(4):
    edges={n:e for bit,n,e in [(0,3,5),(1,5,8)] if mask>>bit&1}
    row={}
    for n,w in [(3,1),(8,-1)]:
        while n in edges:n,w=edges[n],-w
        row[n]=row.get(n,0)+w
    values.append(sum(abs(v) for v in row.values())**2)
need(values==[4,4,4,0],'non-submodular native graph payment')
report={'status':'FINITE-EXACT scoped audit, no all-scale energy bound',
        'source_universe_max':N,'two_octave_identity_N_count':N//6+1,'energy_cases':energies,
        'mean_zero_cylinder_function_mod9':f,
        'L2_values':[str(x) for x in (norm,norm_u,norm_v,cross)],
        'operator_norm_ratio_at_real_three_quarters':'(1+sqrt(3))/(2sqrt(2))+(1/2)(3/4)^(3/4) > 1',
        'branch_square_sum_at_three_quarters':'(1+sqrt(3))/(2sqrt(2)) < 1',
        'native_complementarity_norms_by_mask':values}
p=Path(__file__).resolve().parents[2]/'05-knowledge/results/collatz_incoming_audit_20261015b.json'
p.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='energy_cases'},indent=2))
