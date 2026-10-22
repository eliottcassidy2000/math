"""Exact finite controls for the separate general proofs in SCORE_SCALE_AND_RANDOM_CONTROLS.md."""
import itertools,json
from fractions import Fraction as F
out=[]
for n in range(2,6):
 vals=[];trace4=[]
 for bits in itertools.product([-1,1],repeat=n*(n-1)//2):
  s=[[0]*n for _ in range(n)]
  for v,(i,j) in zip(bits,itertools.combinations(range(n),2)):s[i][j]=v;s[j][i]=-v
  vals.append(sum(sum(r)**2 for r in s))
  s2=[[sum(s[i][k]*s[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
  trace4.append(sum(v*v for r in s2 for v in r))
 mu=F(sum(vals),len(vals));var=F(sum((v-mu)**2 for v in vals),len(vals));t=F(sum(trace4),len(trace4))
 if mu!=n*(n-1) or var!=2*n*(n-1)*(n-2) or t!=n*(n-1)*(2*n-3):raise RuntimeError('failure')
 out.append({'n':n,'count':len(vals),'E_W':str(mu),'Var_W':str(var),'E_trS4':str(t)})
print(json.dumps(out,indent=2))
