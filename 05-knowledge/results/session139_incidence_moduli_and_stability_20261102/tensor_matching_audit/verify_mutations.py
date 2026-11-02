#!/usr/bin/env python3
"""Exact adversarial mutation controls; intended false strengthenings fail."""
from pathlib import Path
from fractions import Fraction as F
import json
# Load only independent helper definitions, not the full positive-control loop.
src=Path(__file__).with_name('verify_exact.py').read_text()
exec(src.split('# Full low-dimensional lattice controls')[0])

def matrix(n,r,bits,x):
    k=[[F(0) for _ in range(n)] for _ in range(n)]
    pairs={(2*i,2*i+1):x[i] for i in range(r)};h=0
    for i in range(n):
        for j in range(i+1,n):
            value=pairs.get((i,j))
            if value is None:value=F(1 if bits>>h&1 else -1);h+=1
            k[i][j]=value;k[j][i]=-value
    return k

out={}
def negative(name,gap,**data):
    require(gap<0,(name,'mutation not killed',gap))
    out[name]={'false_claim_gap':str(gap),**data}

negative('missing_coordinate_domain',tensor([1,-2],[F(2)])-(1+F(2)**2),b=[1,-2],x='2')
negative('missing_even_tail',tensor([1,-1],[F(1)])-2,b=[1,-1],x='1')
negative('missing_odd_constant',tensor([0,0],[F(0)])-1,b=[0,0],x='0')
# Full-pair determinant floors need the parity correction.
k=matrix(2,1,0,[F(0)]);u,v=uv(k_to_a(k));q=F(1,4)
negative('omit_no_isolate_correction',u-abs(v)-(1-q),u=str(u),v=str(v))
negative('wrong_signed_parity_n2',u+v-((1-q)+q),u=str(u),v=str(v))
# A true matching case where individual Pfaffian >= monomial is false.
found=None
for bits in range(1<<14):
    k=matrix(6,1,bits,[F(0)]);pf=pf_factory(k)
    if pf((2,3,4,5))==3 and pf(tuple(range(6)))==-2:
        found=(bits,k);break
require(found is not None,'3x-2 Pfaffian witness missing')
bits,k=found;k[0][1]=F(2,3);k[1][0]=F(-2,3);pf=pf_factory(k)
negative('termwise_pfaffian_floor',pf(tuple(range(6)))**2-k[0][1]**2,binary_cross_bits=bits,skew_pair='2/3',polynomial='3x-2',skew_matrix=[[str(v) for v in row] for row in k])
# Conditional (1-q) times endpoint margin is false even for matching uncertainty.
# Ordinal sum: fair pair, 3-cycle, 3-cycle.
n=8;k=[[F(0) if i==j else F(1 if i<j else -1) for j in range(n)] for i in range(n)]
k[0][1]=k[1][0]=0
for a,b,c in [(2,3,4),(5,6,7)]:k[a][c]=-1;k[c][a]=1
u,v=uv(k_to_a(k));q=F(1,4);endpoint=F(4)
negative('conditional_product_factor',u-abs(v)-(1-q)*endpoint,u=str(u),v=str(v),true_absolute_floor=str(1-q))
# Other fractional cross entries invalidate the matching-specific absolute floor.
C=[[0,1,1],[-1,0,1],[-1,-1,0]];U=[[F(1,2),F(1,2),0],[F(-1,4),F(1,4),F(1,2)]]
k=[[F(0) for _ in range(5)] for _ in range(5)];k[0][1]=F(3,4);k[1][0]=F(-3,4)
for i in range(2):
    for j in range(3):k[i][j+2]=U[i][j];k[j+2][i]=-U[i][j]
for i in range(3):
    for j in range(3):k[i+2][j+2]=C[i][j]
u,v=uv(k_to_a(k));q=F(7,64)
negative('ignore_other_fractional_edges',u-abs(v)-(1-q),u=str(u),v=str(v))
print(json.dumps({'status':'all eight false strengthenings rejected exactly','mutations':out},indent=2,sort_keys=True))
