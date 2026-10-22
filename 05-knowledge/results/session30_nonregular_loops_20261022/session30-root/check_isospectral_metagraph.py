import itertools,json
from fractions import Fraction
from pathlib import Path
import sympy as sp
BITS=['000110001000000','000110011000000']
out=[];n=6
for bits in BITS:
 s=[[0]*n for _ in range(n)]
 for bit,(i,j) in zip(bits,itertools.combinations(range(n),2)):
  s[i][j]=1 if bit=='1' else -1;s[j][i]=-s[i][j]
 H=aut=anti=blue=0
 for p in itertools.permutations(range(n)):
  aut+=all(s[p[i]][p[j]]==s[i][j] for i in range(n) for j in range(i+1,n))
  anti+=all(s[p[i]][p[j]]==-s[i][j] for i in range(n) for j in range(i+1,n))
  if all(s[p[i]][p[i+1]]==1 for i in range(n-1)):
   H+=1
   blue+=all(s[p[n-1-i]][p[n-1-j]]==-s[p[i]][p[j]] for i in range(n) for j in range(i+1,n))
 S=sp.Matrix(s);rho=(n*sp.eye(n)-S).det()/abs((sp.I*n*sp.eye(n)-S).det())
 if sorted(map(sum,s))!=[-1]*3+[1]*3:raise RuntimeError('scores')
 if H%aut or blue%aut or anti!=aut:raise RuntimeError('orbit failure')
 out.append({'bits':bits,'H':H,'Aut':aut,'anti':anti,'h':H//aut,'b':blue//aut,'deficit':(H-blue)//aut,'fundamental_black_degree':(H-blue)//(2*aut),'charpoly':list(map(int,S.charpoly().all_coeffs())),'rho':str(rho)})
if out[0]['rho']!=out[1]['rho'] or out[0]['charpoly']!=out[1]['charpoly']:raise RuntimeError('isospectral failure')
if [(a['H'],a['Aut'],a['h'],a['b']) for a in out]!=[(41,1,41,9),(45,3,15,9)]:raise RuntimeError('unexpected counts')
print(json.dumps({'status':'PASS','universe':'two explicit six-vertex witnesses; all720 permutations for each','witnesses':out},indent=2))
