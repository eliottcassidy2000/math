#!/usr/bin/env python3
"""Independent q=2 rank0/1/2 checks using source/target colors directly."""
from fractions import Fraction as F
from itertools import product
import json

def check(x,msg):
    if not x:raise RuntimeError(msg)
def mul(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def sub(A,B):return [[a-b for a,b in zip(r,s)] for r,s in zip(A,B)]
def scale(A,c):return [[a*c for a in r] for r in A]
def ident(n):return [[F(i==j) for j in range(n)] for i in range(n)]
colors=range(2); arrows=list(product(colors,repeat=2))
# D deletes either source or target. U attaches a fresh vertex in either orientation.
D1=[[F(1),F(1)]];E0=[[F(1,2)],[F(1,2)]]
D2=[[F((a==c)+(b==c)) for a,b in arrows] for c in colors]
U1=[[F((a==c)+(b==c)) for c in colors] for a,b in arrows]
E1=scale(U1,F(1,4))
check(sub(mul(D2,E1),mul(E0,D1))==ident(2),'Weyl rank1')
P1=mul(E0,D1);P2=scale(mul(E1,D2),F(1,2))
check(mul(P1,P1)==P1,'rank1 refresh eigenvalues 1,0')
# Spectrum rank2:1,1/2,0,0. Polynomial plus traces establish multiplicities exactly.
P2mI=sub(P2,ident(4)); P2mhalf=sub(P2,scale(ident(4),F(1,2)))
check(not any(x for row in mul(mul(P2,P2mI),P2mhalf) for x in row),'rank2 spectral polynomial')
check(sum(P2[i][i] for i in range(4))==F(3,2),'rank2 trace')
check(all(sum(row)==1 for row in P2),'stochastic symmetry')
check(all(P2[i][j]==P2[j][i] for i in range(4) for j in range(4)),'uniform rank2 detailed balance')
# All two-vertex color classes contain a unique source; source-free count0 vs kerD dimension2.
check(D2[0][0]*D2[1][3]-D2[0][3]*D2[1][0]!=0,'D2 full row rank')
# Forgetting source color creates a kernel contrast; distinguished-color deletion avoids it.
delta1=[[F(1),F(1)]]; deltac1=[[F(1),F(0)]]; U0=[[F(1)],[F(1)]]
deltac2=[[F(a==0 and b==c) for a,b in arrows] for c in colors]
check(mul(delta1,[[F(1)],[F(-1)]])==[[F(0)]],'source-color contrast in all-source deletion kernel')
check(sub(mul(deltac2,U1),mul(U0,deltac1))==ident(2),'distinguished-source Weyl relation')
check(sum(all(deltac2[i][j]==0 for i in range(2)) for j in range(4))==2,'distinguished-source zero basis count')
print(json.dumps({'status':'PASS','q':2,'a_0_to_2':[1,2,4],'rank1_spectrum':[1,0],'rank2_spectrum':[1,'1/2',0,0],'rank2_source_free':0,'rank2_kernel_D_dimension':2},indent=2))
