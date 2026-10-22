import json,sympy as s
from pathlib import Path; x=json.load(open(Path(__file__).parent/'transfer-results.json'));C={};bad=[]
for t,m in x['matrices'].items():
 rows=[[sum(m[i][:3]),m[i][0],m[i][1]+sum(m[i][3:])] for i in (2,0,1,3)]
 C[t]=s.Matrix.hstack(s.Matrix(rows[0]),s.Matrix([rows[1][j]-rows[0][j] for j in range(3)]),s.Matrix(rows[3]))
 # condition is E coeff rows2-rows0 vs Q coeff rows3
 if s.Matrix(rows[2])-s.Matrix(rows[0])!=s.Matrix(rows[3]):bad.append(('EvsQ',t))
if bad: raise RuntimeError(str(bad))
print('COMMON THREE-STATE COMPRESSION PASS all18')
# mirror map computed without heavy imports
import itertools as it
E=list(it.combinations(range(5),2))
def arc(m,i,j):return (m>>E.index((i,j)))&1 if i<j else 1-arc(m,j,i)
def tau(m):return sum((1-arc(m,4-i,4-j))<<b for b,(i,j) in enumerate(E))
G=s.diag(1,-1,1)
if any(C[str(tau(int(t)))]!=G*M.T*G for t,M in C.items()): raise RuntimeError('adjoint')
print('PER LETTER MIRROR ADJOINT PASS all18')
if any(M[:,0]!=s.Matrix([sum(x['initial'][t][:3]),x['initial'][t][0],x['initial'][t][1]+sum(x['initial'][t][3:])]) for t,M in C.items()): raise RuntimeError('initial')
print('HOMOGENEOUS INITIALIZATION PASS all18')
open(Path(__file__).parent/'compressed-matrices.json','w').write(json.dumps({t:M.tolist() for t,M in C.items()},default=int))
