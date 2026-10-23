import itertools, sympy as s
z=s.symbols('z')
def require(condition):
 if not condition: raise RuntimeError('Exact identity check failed')
def D(M):
 n=M.rows; I=s.eye(n)
 return M.row_join(M+I).col_join((M-I).row_join(-M))
def E(M):
 n=M.rows; o=s.ones(1,n)
 return s.Matrix([[0]]).row_join(o).row_join(s.ones(1,1)).col_join((-o.T).row_join(M).row_join(o.T)).col_join(s.Matrix([[-1]]).row_join(-o).row_join(s.zeros(1,1)))
count=0
for n in range(1,5):
 for bits in itertools.product([-1,1],repeat=n*(n-1)//2):
  M=s.zeros(n)
  for (i,j),a in zip(itertools.combinations(range(n),2),bits): M[i,j]=a;M[j,i]=-a
  Q=s.eye(n)-M*M; N=D(M)
  require(s.eye(2*n)-N*N==s.diag(2*Q,2*Q))
  B=z*s.eye(n)-M
  require(s.expand(E(M).charpoly(z).as_expr()-((z*z+1)*B.det()+2*z*(s.ones(1,n)*B.adjugate()*s.ones(n,1))[0]))==0)
  count+=1
print('Exact checks: all',count,'labeled tournaments n=1..4, Q doubling and E determinant pass.')
M=s.zeros(1)
for name,N in [('D(E(seed))',D(E(M))),('E(E(D(seed)))',E(E(D(M))))]:
 print(name,'scores',sorted([int((sum(N.row(i))+N.rows-1)/2) for i in range(N.rows)]),'charpoly',s.factor(N.charpoly(z).as_expr()))
# isospectral same-size tournaments need resolvent sidecar? exhaustive n=4 group.
groups={}
for bits in itertools.product([-1,1],repeat=6):
 M=s.zeros(4)
 for (i,j),a in zip(itertools.combinations(range(4),2),bits):M[i,j]=a;M[j,i]=-a
 p=str(M.charpoly(z).as_expr()); pe=str(E(M).charpoly(z).as_expr())
 groups.setdefault(p,set()).add(pe)
print('n=4 spectra extension fibers:',{k: sorted(v) for k,v in sorted(groups.items())})
seen={}
for bits in itertools.product([-1,1],repeat=6):
 M=s.zeros(4)
 for (i,j),a in zip(itertools.combinations(range(4),2),bits):M[i,j]=a;M[j,i]=-a
 p=str(M.charpoly(z).as_expr()); pe=str(E(M).charpoly(z).as_expr())
 key=(p,pe)
 if key not in seen:
  seen[key]=bits
  print('sidecar representative',bits,'base',p,'extended',pe,'numerator',s.factor((s.ones(1,4)*(z*s.eye(4)-M).adjugate()*s.ones(4,1))[0]))
for bits in itertools.product([-1,1],repeat=3):
 M=s.zeros(3)
 for (i,j),a in zip(itertools.combinations(range(3),2),bits):M[i,j]=a;M[j,i]=-a
 def F(N):return (s.ones(1,N.rows)*(z*s.eye(N.rows)-N).inv()*s.ones(N.rows,1))[0]
 f=s.factor(F(M)); fe=s.factor(F(E(M)))
 require(s.cancel(fe-((z*z+1)*f+2*z)/(z*z+1+2*z*f))==0)
 R=((z*z+1)*s.eye(3)-2*M*M).inv()
 require(s.cancel(F(D(M))-2*z*(s.ones(1,3)*R*s.ones(3,1))[0])==0)
print('Eight n=3 tournaments: closed marked resolvent E and D recursions verified exactly.')
