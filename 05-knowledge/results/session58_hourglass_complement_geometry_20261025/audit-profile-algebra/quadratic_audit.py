import json
from pathlib import Path
from fractions import Fraction as F
D=Path(__file__).parent
q=json.loads((Path(__file__).resolve().parent.parent/'observer-profile-algebra/quadratic_certificate.json').read_text());a=json.loads((D/'independent_results.json').read_text())['census'][-1]
def req(v,m):
 if not v:raise RuntimeError(m)
features=q['selected_features'];cols=a['cols'];A=a['matrix']
M=[[row[cols.index(f[1])]**(1 if f[0]=='linear' else 2) for f in features] for row in A]
req(M==q['selected_integer_minor'],'independent reconstruction of minor')
B=[[F(x) for x in row] for row in M];det=F(1)
for j in range(len(B)):
 k=next((i for i in range(j,len(B)) if B[i][j]),None);req(k is not None,'full rank')
 if k!=j:B[k],B[j]=B[j],B[k];det=-det
 pivot=B[j][j];det*=pivot
 for i in range(j+1,len(B)):
  ratio=B[i][j]/pivot
  B[i]=[x-ratio*y for x,y in zip(B[i],B[j])]
req(str(det)==q['determinant'],'exact determinant')
print(json.dumps({'rank':34,'determinant':str(det),'square_columns':[f[1] for f in features if f[0]=='square']},indent=2))
