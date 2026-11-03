"""Fixed processing-order alphabets and actual-source constructors."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
PRIOR=HERE.parent/'return_channels_20261103'
sys.path.insert(0,str(PRIOR))
from channels import SEEDS,T,LOOPS
from common import ATLAS,MATS as OLD_MATS,COSTS as OLD_COSTS,need,step,actual_target,literal_frontier
MATS=dict(OLD_MATS);MATS['s78']=tuple(map(tuple,ATLAS['s78']['matrix']))
COSTS=dict(OLD_COSTS);COSTS['s78']=6
MACROS={'b':('b',),'q':('q2',),'s':('s78',),**LOOPS}
D=LOOPS['L2']+('b',)+LOOPS['L2']
REPAIRS=((287,206,239),(289,191,239),(293,210,281),(295,216,260))

def vector(word):
    v=(1,0,0)
    for l in word:v=step(v,MATS[l])
    return v

def matrix(word):
    p=((1,0,0),(0,1,0),(0,0,1))
    for l in word:
        m=MATS[l];p=tuple(tuple(sum(p[i][k]*m[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    return p

def bulk(v):return (v[2],v[0]+2*v[1]-v[2])

def quotient(p):
    z=[(r[2],r[0]+2*r[1]-r[2]) for r in p]
    q=(tuple(z[0][j]+z[2][j] for j in range(2)),z[0])
    need(z[1]==tuple(2*x for x in q[1]),'error-line quotient identity')
    return q

def det2(q):return q[0][0]*q[1][1]-q[0][1]*q[1][0]

def inverse_row(v,q):
    from fractions import Fraction as F
    a,b=q[0];c,d=q[1];det=det2(q)
    return (F(v[0]*d-v[1]*c,det),F(-v[0]*b+v[1]*a,det))
