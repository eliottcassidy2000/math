from sympy import Matrix,ones
import json

def basis(cols):return Matrix.hstack(*cols).columnspace()
def closure(A,cols,joint):
 B=basis([ones(A.rows,1)]+cols)
 while True:
  trial=B+[A*v for v in B]
  if joint:trial += [v.multiply_elementwise(w) for i,v in enumerate(B) for w in B[i:]]
  new=basis(trial)
  if len(new)==len(B):return B
  B=new
A=Matrix.zeros(6)
for i,j in ((0,2),(0,5),(1,3),(1,4)):A[i,j]=1
lv=Matrix([1,1,0,0,0,0]);a=Matrix([0,0,0,0,1,1]);b=Matrix([0,0,0,1,0,1])
ranks=[len(closure(A,[lv,a,b],joint)) for joint in (False,True)]
if ranks!=[4,6]:raise RuntimeError('six-state ranks')
if list(A*a.multiply_elementwise(b))[:2]!=[1,0]:raise RuntimeError('six-state separator')
A=Matrix.zeros(7)
for i,j,w in ((0,2,1),(0,3,1),(1,4,2),(2,5,2),(3,6,2),(4,5,1),(4,6,1)):A[i,j]=w
f=Matrix([0,0,0,0,0,1,0]);lv=[Matrix([int(i in S) for i in range(7)]) for S in ((0,1),(2,3,4),(5,6))]
B=closure(A,[f]+lv,False)
if not all(v[0]==v[1] for v in B):raise RuntimeError('linear blindness')
g=A*f;sep=list(A*g.multiply_elementwise(g))[:2]
if sep!=[4,2]:raise RuntimeError('variance')
print(json.dumps({'six_state_linear_joint_dimensions':ranks,'seven_state_separator':list(map(int,sep))}))
